"""Tests for the MCP server using the SDK's in-memory client."""

import pytest
from mcp.client.client import Client

from pubchem_finder.mcp_server import mcp

pytestmark = pytest.mark.anyio

EXPECTED_TOOLS = {
    "pubchem_lookup",
    "compound_properties",
    "standardize",
    "compare",
    "similarity_search",
    "index_status",
}


async def test_tools_are_listed_with_schemas():
    async with Client(mcp) as client:
        result = await client.list_tools()
    tools = {t.name: t for t in result.tools}
    assert set(tools) == EXPECTED_TOOLS
    for tool in tools.values():
        assert tool.description, tool.name
        assert tool.output_schema is not None, tool.name
    assert "smiles" in tools["compound_properties"].input_schema["properties"]
    assert tools["pubchem_lookup"].input_schema["properties"]["by"]["enum"] == [
        "name",
        "cid",
        "smiles",
        "inchikey",
    ]


async def test_compound_properties_structured(aspirin_smiles):
    async with Client(mcp) as client:
        result = await client.call_tool(
            "compound_properties", {"smiles": aspirin_smiles}
        )
    assert result.is_error is False
    assert result.structured_content["canonical_smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert result.structured_content["lipinski_ok"] is True


async def test_invalid_smiles_is_tool_error(invalid_smiles):
    async with Client(mcp) as client:
        result = await client.call_tool(
            "compound_properties", {"smiles": invalid_smiles}
        )
    assert result.is_error is True
    assert "Invalid SMILES" in result.content[0].text


async def test_standardize_and_compare():
    async with Client(mcp) as client:
        std = await client.call_tool("standardize", {"smiles": "CC(=O)[O-].[Na+]"})
        cmp = await client.call_tool("compare", {"smiles_a": "CCO", "smiles_b": "CCO"})
    assert std.structured_content["smiles"] == "CC(=O)O"
    assert cmp.structured_content["ecfp4_tanimoto"] == 1.0


async def test_pubchem_lookup(pubchem_api, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    async with Client(mcp) as client:
        ok = await client.call_tool("pubchem_lookup", {"query": "caffeine"})
        missing = await client.call_tool("pubchem_lookup", {"query": "notacompound"})
    assert ok.structured_content["compounds"][0]["cid"] == 2519
    assert ok.structured_content["compounds"][0]["pubchem_url"].endswith("/2519")
    assert missing.is_error is True
    assert "Check the spelling" in missing.content[0].text


async def test_similarity_search_and_status(built_home, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(built_home))
    async with Client(mcp) as client:
        status = await client.call_tool("index_status", {})
        result = await client.call_tool(
            "similarity_search", {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "threshold": 0.3}
        )
    assert status.structured_content["state"] == "ready"
    assert result.is_error is False
    top = result.structured_content["hits"][0]
    assert (top["cid"], top["similarity"]) == (2244, 1.0)
    assert top["pubchem_url"].endswith("/2244")


async def test_similarity_search_without_index(tmp_path, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(tmp_path))
    async with Client(mcp) as client:
        status = await client.call_tool("index_status", {})
        result = await client.call_tool("similarity_search", {"smiles": "CCO"})
    assert status.structured_content["state"] == "not_built"
    assert result.is_error is True
    assert "pubchem-finder index build" in result.content[0].text
