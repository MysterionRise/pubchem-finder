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
            prompts = await session.list_prompts()
            templates = await session.list_resource_templates()
            status = await session.read_resource("pubchem://index/status")
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
    assert init.server_info.title == "PubChem Finder"
    assert init.server_info.icons[0].src.startswith("https://")
    assert len(tools.tools) == 17
    assert len(prompts.prompts) == 3
    assert templates.resource_templates[0].uri_template == "pubchem://compound/{cid}"
    assert '"state": "ready"' in status.contents[0].text
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


async def test_first_search_builds_the_demo_index_over_stdio(tmp_path):
    """No index: the first search builds the demo index inside the server, and the
    build writes nothing to stdout, the transport."""
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "pubchem_finder", "mcp"],
        env={
            **os.environ,
            "PUBCHEM_FINDER_HOME": str(tmp_path / "home"),
            "PUBCHEM_FINDER_CACHE": "0",
            "PUBCHEM_FINDER_DEMO": "1",
        },
    )
    with (tmp_path / "server-stderr.log").open("w", encoding="utf-8") as errlog:
        async with (
            stdio_client(params, errlog=errlog) as (read, write),
            ClientSession(read, write) as session,
        ):
            await session.initialize()
            status = await session.call_tool("index_status", {})
            similar = await session.call_tool(
                "similarity_search", {"smiles": ASPIRIN, "threshold": 0.4}
            )
            fragment = await session.call_tool(
                "substructure_search", {"query": "c1ccc2[nH]ccc2c1", "max_hits": 5}
            )

    assert status.structured_content["search_scope"] == "demo"
    assert similar.structured_content["search_scope"] == "demo"
    assert similar.structured_content["hits"][0]["cid"] == 2244
    assert len(fragment.structured_content["hits"]) == 5
    assert (tmp_path / "home" / "demo").is_dir()
