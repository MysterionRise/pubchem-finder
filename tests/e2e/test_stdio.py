"""The real MCP server process over stdio, the way Claude Desktop and Claude Code run it.

The in-memory tests in tests/unit exercise the tools; this one proves the
installed entry point starts, speaks the protocol on stdout and nothing else
writes there.
"""

import os
import sys

import pytest
from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

from pubchem_finder import __version__

pytestmark = pytest.mark.anyio

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


async def test_stdio_server_end_to_end(substructure_home, tmp_path):
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "pubchem_finder", "mcp"],
        env={
            **os.environ,
            "PUBCHEM_FINDER_HOME": str(substructure_home),
            "PUBCHEM_FINDER_CACHE": "0",
        },
    )
    with (tmp_path / "server-stderr.log").open("w", encoding="utf-8") as errlog:
        async with (
            stdio_client(params, errlog=errlog) as (read, write),
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
                    "smiles": [
                        h["smiles"] for h in fragment.structured_content["hits"]
                    ],
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
