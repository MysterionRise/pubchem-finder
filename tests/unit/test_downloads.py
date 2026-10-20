"""Tests for the PubChem SDF mirror client (HTTP mocked with respx)."""

import hashlib
import os
import zlib

import httpx
import pytest
import respx

from pubchem_finder.core.exceptions import DownloadError
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
    with SdfMirror() as m:
        yield m


@pytest.fixture
def api():
    with respx.mock(assert_all_called=False) as router:
        router.get(URL).respond(200, text=LISTING)
        router.get(URL + NAME + ".md5").respond(200, text=f"{MD5}  {NAME}\n")
        yield router


def test_list_chunks(api, mirror):
    assert mirror.list_chunks() == [NAME, "Compound_000500001_001000000.sdf.gz"]


def test_list_chunks_errors(mirror):
    with respx.mock() as router:
        router.get(URL).respond(200, text="<html>nothing here</html>")
        with pytest.raises(DownloadError, match="No Compound"):
            mirror.list_chunks()
    with respx.mock() as router:
        router.get(URL).respond(503)
        with pytest.raises(DownloadError, match="Could not list"):
            mirror.list_chunks()


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
    api.get(URL + NAME).respond(200, content=PAYLOAD[:-1] + b"x")
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
