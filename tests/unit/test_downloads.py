"""Tests for the PubChem SDF mirror client (HTTP mocked with respx)."""

import hashlib
import os
import zlib

import httpx
import pytest
import respx

from pubchem_finder.core.exceptions import DownloadError, InvalidInput
from pubchem_finder.data.downloads import DEFAULT_SDF_URL, SdfMirror, file_md5

NAME = "Compound_000000001_000500000.sdf.gz"
PAYLOAD = os.urandom(3 * 1024 * 1024 + 123)  # a few blocks, not block-aligned
MD5 = hashlib.md5(PAYLOAD).hexdigest()
URL = DEFAULT_SDF_URL
LISTING = f"""<html><body>
<a href="../">../</a>
<a href="Compound_000500001_001000000.sdf.gz">Compound_000500001_001000000.sdf.gz</a>
<a href="{NAME}">{NAME}</a>
<a href="{NAME}.md5">{NAME}.md5</a>
<a href="{NAME}">{NAME}</a>
<a href="README">README</a>
</body></html>"""


@pytest.fixture
def mirror():
    # The suite must stay instant: the retries are real, the pauses are not.
    with SdfMirror(retry_delay=0) as m:
        yield m


@pytest.fixture
def api():
    with respx.mock(assert_all_called=False) as router:
        router.get(URL).respond(200, text=LISTING)
        router.get(URL + NAME + ".md5").respond(200, text=f"{MD5}  {NAME}\n")
        yield router


def test_list_chunks(api, mirror):
    assert mirror.list_chunks() == [NAME, "Compound_000500001_001000000.sdf.gz"]


def test_mirror_keeps_no_cookies():
    """One client downloads every chunk of a build; NCBI's session cookie must not
    go back with the next request."""
    cookie = {"set-cookie": "ncbi_sid=ABC123; domain=.nih.gov; path=/"}
    with respx.mock() as router, SdfMirror() as mirror:
        router.get(DEFAULT_SDF_URL).respond(
            200, text=f'<a href="{NAME}">{NAME}</a>', headers=cookie
        )
        md5 = router.get(DEFAULT_SDF_URL + NAME + ".md5").respond(
            200, text=f"{'0' * 32}  {NAME}\n", headers=cookie
        )
        assert mirror.list_chunks() == [NAME]
        mirror.expected_md5(NAME)
        mirror.expected_md5(NAME)
        assert [call.request.headers.get("cookie") for call in md5.calls] == [
            None,
            None,
        ]
        assert len(mirror._client.cookies.jar) == 0


def test_list_chunks_errors(mirror):
    with respx.mock() as router:
        router.get(URL).respond(200, text="<html>nothing here</html>")
        with pytest.raises(DownloadError, match="No Compound"):
            mirror.list_chunks()
    with respx.mock() as router:
        router.get(URL).respond(503)
        with pytest.raises(DownloadError, match="Could not list"):
            mirror.list_chunks()


def test_list_chunks_retries_a_transient_failure(mirror):
    """The listing is the build's first request, so a blip there must not end it."""
    with respx.mock() as router:
        route = router.get(URL).mock(
            side_effect=[httpx.Response(503), httpx.Response(200, text=LISTING)]
        )
        assert mirror.list_chunks() == [NAME, "Compound_000500001_001000000.sdf.gz"]
    assert route.call_count == 2


def test_list_chunks_does_not_retry_a_permanent_failure(mirror):
    with respx.mock() as router:
        route = router.get(URL).respond(404)
        with pytest.raises(DownloadError, match="Could not list"):
            mirror.list_chunks()
    assert route.call_count == 1


def test_download_verifies_checksum(api, mirror, tmp_path):
    api.get(URL + NAME).respond(200, content=PAYLOAD)
    seen: list[int] = []
    path = mirror.download(
        NAME, tmp_path, progress=lambda n, done, total: seen.append(done)
    )
    assert path == tmp_path / NAME
    assert path.read_bytes() == PAYLOAD
    assert not (tmp_path / f"{NAME}.part").exists()
    assert seen[-1] == len(PAYLOAD)


def test_checksum_mismatch(api, mirror, tmp_path):
    # Flip the last byte: replacing it with a fixed byte is a no-op whenever the
    # random payload already ends in that byte.
    corrupt = PAYLOAD[:-1] + bytes([PAYLOAD[-1] ^ 0xFF])
    api.get(URL + NAME).respond(200, content=corrupt)
    with pytest.raises(DownloadError, match="Checksum mismatch"):
        mirror.download(NAME, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_existing_verified_file_is_reused(api, mirror, tmp_path):
    (tmp_path / NAME).write_bytes(PAYLOAD)
    route = api.get(URL + NAME).respond(200, content=b"should not be fetched")
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD
    assert not route.called


def test_existing_corrupt_file_is_replaced(api, mirror, tmp_path):
    (tmp_path / NAME).write_bytes(b"corrupt")
    api.get(URL + NAME).respond(200, content=PAYLOAD)
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD


def test_resume_with_range(api, mirror, tmp_path):
    half = len(PAYLOAD) // 2
    (tmp_path / f"{NAME}.part").write_bytes(PAYLOAD[:half])

    def respond(request: httpx.Request) -> httpx.Response:
        assert request.headers["Range"] == f"bytes={half}-"
        return httpx.Response(206, content=PAYLOAD[half:])

    api.get(URL + NAME).mock(side_effect=respond)
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD


def test_server_ignoring_range_restarts(api, mirror, tmp_path):
    (tmp_path / f"{NAME}.part").write_bytes(b"stale partial bytes")
    api.get(URL + NAME).respond(200, content=PAYLOAD)
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD


def test_unsatisfiable_range_restarts_once(api, mirror, tmp_path):
    (tmp_path / f"{NAME}.part").write_bytes(PAYLOAD + b"extra")
    api.get(URL + NAME).mock(
        side_effect=[httpx.Response(416), httpx.Response(200, content=PAYLOAD)]
    )
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD


def test_unsatisfiable_range_twice_fails(api, mirror, tmp_path):
    (tmp_path / f"{NAME}.part").write_bytes(b"x")
    api.get(URL + NAME).mock(side_effect=[httpx.Response(416), httpx.Response(416)])
    with pytest.raises(DownloadError, match="refused to resume"):
        mirror.download(NAME, tmp_path)


def test_http_and_transport_errors(api, mirror, tmp_path):
    api.get(URL + NAME).respond(500)
    with pytest.raises(DownloadError, match="HTTP 500"):
        mirror.download(NAME, tmp_path)
    api.get(URL + NAME).mock(side_effect=httpx.ReadTimeout("slow"))
    with pytest.raises(DownloadError, match="failed"):
        mirror.download(NAME, tmp_path)


def test_transient_status_is_retried(api, mirror, tmp_path):
    route = api.get(URL + NAME).mock(
        side_effect=[
            httpx.Response(503),
            httpx.Response(429),
            httpx.Response(200, content=PAYLOAD),
        ]
    )
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD
    assert route.call_count == 3


def test_transport_error_is_retried(api, mirror, tmp_path):
    route = api.get(URL + NAME).mock(
        side_effect=[
            httpx.ReadError("connection reset"),
            httpx.Response(200, content=PAYLOAD),
        ]
    )
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD
    assert route.call_count == 2


def test_retry_resumes_where_the_first_attempt_stopped(api, mirror, tmp_path):
    half = len(PAYLOAD) // 2
    (tmp_path / f"{NAME}.part").write_bytes(PAYLOAD[:half])
    ranges: list[str | None] = []

    def second(request: httpx.Request) -> httpx.Response:
        ranges.append(request.headers.get("Range"))
        return httpx.Response(206, content=PAYLOAD[half:])

    api.get(URL + NAME).mock(
        side_effect=[httpx.ReadTimeout("slow"), second],
    )
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD
    assert ranges == [f"bytes={half}-"]


def test_transient_errors_give_up_after_all_attempts(api, mirror, tmp_path):
    route = api.get(URL + NAME).respond(503)
    with pytest.raises(DownloadError, match="HTTP 503; gave up after 4 attempts"):
        mirror.download(NAME, tmp_path)
    assert route.call_count == 4


def test_permanent_status_is_not_retried(api, mirror, tmp_path):
    route = api.get(URL + NAME).respond(404)
    with pytest.raises(DownloadError, match="HTTP 404"):
        mirror.download(NAME, tmp_path)
    assert route.call_count == 1


def test_transient_checksum_failure_is_retried(api, mirror, tmp_path):
    api.get(URL + NAME + ".md5").mock(
        side_effect=[httpx.Response(503), httpx.Response(200, text=f"{MD5}  {NAME}\n")]
    )
    api.get(URL + NAME).respond(200, content=PAYLOAD)
    assert mirror.download(NAME, tmp_path).read_bytes() == PAYLOAD


def test_pauses_double_between_attempts(api, tmp_path):
    pauses: list[float] = []
    api.get(URL + NAME).respond(503)
    with (
        SdfMirror(attempts=4, retry_delay=2.0, sleep=pauses.append) as m,
        pytest.raises(DownloadError, match="gave up after 4 attempts"),
    ):
        m.download(NAME, tmp_path)
    assert pauses == [2.0, 4.0, 8.0]


def test_default_pauses_are_five_ten_and_twenty_seconds(api, tmp_path):
    """The shipped default: four attempts, paused 5, 10 and 20 seconds apart."""
    pauses: list[float] = []
    api.get(URL + NAME).respond(503)
    with (
        SdfMirror(sleep=pauses.append) as m,
        pytest.raises(DownloadError, match="gave up after 4 attempts"),
    ):
        m.download(NAME, tmp_path)
    assert pauses == [5.0, 10.0, 20.0]


def test_attempts_below_one_is_refused():
    with pytest.raises(InvalidInput, match="attempts must be at least 1"):
        SdfMirror(attempts=0)


def test_bad_checksum_file(mirror, tmp_path):
    with respx.mock() as router:
        router.get(URL + NAME + ".md5").respond(200, text="not a checksum")
        with pytest.raises(DownloadError, match="Unrecognized checksum"):
            mirror.download(NAME, tmp_path)
    with respx.mock() as router:
        router.get(URL + NAME + ".md5").respond(404)
        with pytest.raises(DownloadError, match="Could not fetch checksum"):
            mirror.download(NAME, tmp_path)


def test_injected_client_and_base_url_normalization(tmp_path):
    client = httpx.Client()
    mirror = SdfMirror(base_url="https://example.org/sdf", client=client)
    assert mirror.base_url == "https://example.org/sdf/"
    mirror.close()
    assert not client.is_closed
    client.close()


def test_file_md5(tmp_path):
    path = tmp_path / "f"
    path.write_bytes(PAYLOAD)
    assert file_md5(path) == MD5


@pytest.mark.live
@pytest.mark.skipif(not os.environ.get("PUBCHEM_LIVE"), reason="set PUBCHEM_LIVE=1")
def test_live_listing_and_first_record_tags():
    """The mirror lists chunks and the first record carries the tags we read."""
    with SdfMirror() as mirror:
        names = mirror.list_chunks()
        assert names[0].startswith("Compound_000000001_")
        assert len(mirror.expected_md5(names[0])) == 32
        with httpx.Client(follow_redirects=True, timeout=60) as client:
            head = client.get(
                DEFAULT_SDF_URL + names[0], headers={"Range": "bytes=0-262143"}
            ).content
    text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(head).decode("latin-1")
    first = text.split("$$$$")[0]
    for tag in (
        "PUBCHEM_COMPOUND_CID",
        "PUBCHEM_IUPAC_INCHIKEY",
        "PUBCHEM_MOLECULAR_WEIGHT",
    ):
        assert f"<{tag}>" in first, tag
    assert "<PUBCHEM_SMILES>" in first or "<PUBCHEM_OPENEYE_ISO_SMILES>" in first
