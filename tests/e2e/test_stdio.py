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


@pytest.mark.parametrize("source", ["option", "environment", "module"])
@pytest.mark.parametrize("warning_filter", ["default", "error", "ignore"])
async def test_relative_home_warns_without_corrupting_stdio(
    tmp_path, source, warning_filter
):
    # Run the real entry point, tools and HTTP client, replacing only the network.
    # Unexpected requests fail rather than reaching PubChem.
    mocked_server = """
import runpy
import sys
import respx

module = sys.argv.pop(1)
with respx.mock as router:
    router.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/JSON").respond(
        200, json={"InformationList": {"Information": [{"CID": 2244, "Synonym": ["aspirin"]}]}}
    )
    runpy.run_module(module, run_name="__main__")
"""
    args = ["pubchem_finder", "mcp"]
    env = {
        **os.environ,
        "PUBCHEM_FINDER_HOME": "relative-data",
        "PUBCHEM_FINDER_CACHE": "0",
        "PYTHONWARNINGS": warning_filter,
    }
    if source == "option":
        args += ["--home", "relative-data"]
        env["PUBCHEM_FINDER_HOME"] = str(tmp_path / "ignored")
    elif source == "module":
        args = ["pubchem_finder.mcp_server"]
    params = StdioServerParameters(
        command=sys.executable, args=["-c", mocked_server, *args], env=env, cwd=tmp_path
    )
    log = tmp_path / "server-stderr.log"
    with log.open("w", encoding="utf-8") as errlog:
        async with (
            stdio_client(params, errlog=errlog) as (read, write),
            ClientSession(read, write) as session,
        ):
            await session.initialize()
            for _ in range(2):
                status = await session.call_tool("index_status", {})
                assert status.is_error is False
                assert status.structured_content["home"] == str(
                    tmp_path / "relative-data"
                )
                synonyms = await session.call_tool("pubchem_synonyms", {"cid": 2244})
                assert synonyms.is_error is False
                assert synonyms.structured_content["synonyms"] == ["aspirin"]
    stderr = log.read_text()
    assert stderr.startswith("warning: The relative data directory")
    assert str(tmp_path / "relative-data") in stderr
    assert stderr.count("relative data directory") == 1
    assert stderr.count("warning: ") == 1
    assert "UserWarning" not in stderr
    assert not (tmp_path / "relative-data").exists()
