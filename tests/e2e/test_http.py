"""The real MCP server process over streamable HTTP, the way a URL-based client or
a shared team server uses it (tests/e2e/test_stdio.py covers stdio).

Starts ``pubchem-finder mcp --transport http`` on a free port, waits for
``GET /health``, then runs the same calls as the stdio test through the official
client, plus several sessions at once.
"""

import os
import socket
import subprocess
import sys
import time
from collections.abc import Iterator
from pathlib import Path

import anyio
import httpx
import pytest
from mcp.client.session import ClientSession
from mcp.client.streamable_http import streamable_http_client

from pubchem_finder import __version__

pytestmark = pytest.mark.anyio

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
# Importing RDKit, PyTables and DuckDB takes a while on Windows runners.
STARTUP_TIMEOUT = 120.0


def _free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def _tail(log_path: Path) -> str:
    return log_path.read_text(encoding="utf-8", errors="replace")[-2000:]


def _wait_until_healthy(proc: subprocess.Popen, port: int, log_path: Path) -> None:
    deadline = time.monotonic() + STARTUP_TIMEOUT
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            pytest.fail(f"server exited with {proc.returncode}:\n{_tail(log_path)}")
        try:
            response = httpx.get(f"http://127.0.0.1:{port}/health", timeout=1.0)
        except httpx.HTTPError:
            time.sleep(0.2)
            continue
        if response.status_code == 200 and response.json()["version"] == __version__:
            return
        time.sleep(0.2)
    pytest.fail(f"server not ready after {STARTUP_TIMEOUT:.0f}s:\n{_tail(log_path)}")


@pytest.fixture(scope="module")
def http_server(substructure_home, tmp_path_factory) -> Iterator[tuple[str, Path]]:
    """Start the server on a free port; yield (MCP URL, path of its stderr log)."""
    port = _free_port()
    log_path = tmp_path_factory.mktemp("http") / "server-stderr.log"
    with log_path.open("w", encoding="utf-8") as errlog:
        proc = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "pubchem_finder",
                "mcp",
                "--transport",
                "http",
                "--port",
                str(port),
            ],
            env={
                **os.environ,
                "PUBCHEM_FINDER_HOME": str(substructure_home),
                "PUBCHEM_FINDER_CACHE": "0",
                "PYTHONIOENCODING": "utf-8",
            },
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=errlog,
        )
        try:
            _wait_until_healthy(proc, port, log_path)
            yield f"http://127.0.0.1:{port}/mcp", log_path
        finally:
            proc.terminate()
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()


async def test_http_server_end_to_end(http_server):
    url, log_path = http_server
    async with (
        streamable_http_client(url) as (read, write),
        ClientSession(read, write) as session,
    ):
        init = await session.initialize()
        tools = await session.list_tools()
        props = await session.call_tool("compound_properties", {"smiles": ASPIRIN})
        similar = await session.call_tool(
            "similarity_search", {"smiles": ASPIRIN, "threshold": 0.3}
        )
        fragment = await session.call_tool(
            "substructure_search", {"query": "c1ccccc1C(=O)O"}
        )
        image = await session.call_tool(
            "draw_molecules",
            {
                "smiles": [h["smiles"] for h in fragment.structured_content["hits"]],
                "highlight": "smarts",
                "smarts": "c1ccccc1C(=O)O",
            },
        )
        bad = await session.call_tool("compound_properties", {"smiles": "C1CC"})

    assert init.server_info.name == "pubchem-finder"
    assert init.server_info.version == __version__
    assert len(tools.tools) == 17
    assert props.structured_content["canonical_smiles"] == ASPIRIN
    assert similar.structured_content["hits"][0]["cid"] == 2244
    assert [h["cid"] for h in fragment.structured_content["hits"]] == [
        243,
        338,
        2244,
        4133,
    ]
    assert image.is_error is False
    assert image.content[0].mime_type == "image/png"
    assert bad.is_error is True
    assert url in _tail(log_path)  # the server announced its address on stderr


async def test_http_server_serves_concurrent_sessions(http_server):
    url, _ = http_server
    top_hits: dict[int, int] = {}

    async def one_session(n: int) -> None:
        async with (
            streamable_http_client(url) as (read, write),
            ClientSession(read, write) as session,
        ):
            await session.initialize()
            props = await session.call_tool("compound_properties", {"smiles": ASPIRIN})
            similar = await session.call_tool(
                "similarity_search", {"smiles": ASPIRIN, "threshold": 0.3}
            )
        assert props.structured_content["canonical_smiles"] == ASPIRIN
        top_hits[n] = similar.structured_content["hits"][0]["cid"]

    async with anyio.create_task_group() as tg:
        for n in range(3):
            tg.start_soon(one_session, n)

    assert top_hits == {0: 2244, 1: 2244, 2: 2244}
