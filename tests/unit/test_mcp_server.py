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
    "substructure_search",
    "index_status",
    "pubchem_similarity_search",
    "pubchem_substructure_search",
    "pubchem_exact_match",
    "pubchem_synonyms",
    "pubchem_xrefs",
    "pubchem_description",
    "pubchem_annotations",
    "pubchem_bioactivity",
    "draw_molecules",
    "medchem_profile",
}


async def test_tools_are_listed_with_schemas():
    async with Client(mcp) as client:
        result = await client.list_tools()
    tools = {t.name: t for t in result.tools}
    assert set(tools) == EXPECTED_TOOLS
    for tool in tools.values():
        assert tool.description, tool.name
        if tool.name != "draw_molecules":  # returns image content, not JSON
            assert tool.output_schema is not None, tool.name
    assert "smiles" in tools["compound_properties"].input_schema["properties"]
    assert tools["pubchem_lookup"].input_schema["properties"]["by"]["enum"] == [
        "name",
        "cid",
        "smiles",
        "inchikey",
    ]


async def test_every_tool_is_annotated_read_only_with_a_title():
    """Anthropic's Software Directory Policy requires readOnlyHint, destructiveHint
    and title on every tool; openWorldHint marks the tools that call PubChem."""
    async with Client(mcp) as client:
        result = await client.list_tools()
    titles = []
    for tool in result.tools:
        notes = tool.annotations
        assert notes is not None, tool.name
        assert tool.title, tool.name
        assert tool.title == notes.title, tool.name
        assert notes.read_only_hint is True, tool.name
        assert notes.destructive_hint is False, tool.name
        assert notes.idempotent_hint is True, tool.name
        assert notes.open_world_hint is tool.name.startswith("pubchem_"), tool.name
        wire = tool.model_dump(by_alias=True, exclude_none=True)["annotations"]
        assert {"title", "readOnlyHint", "destructiveHint"} <= set(wire), tool.name
        titles.append(tool.title)
    assert len(titles) == len(set(titles)) == len(EXPECTED_TOOLS)


async def test_unexpected_errors_name_their_cause(monkeypatch, aspirin_smiles):
    def full_disk(smiles):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr("pubchem_finder.service.properties", full_disk)
    async with Client(mcp) as client:
        result = await client.call_tool(
            "compound_properties", {"smiles": aspirin_smiles}
        )
    assert result.is_error is True
    text = result.content[0].text
    assert "OSError" in text
    assert "No space left" in text
    assert "/issues" in text


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
    assert "index build --limit 1" in result.content[0].text


async def test_substructure_search_tool(substructure_home, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(substructure_home))
    async with Client(mcp) as client:
        status = await client.call_tool("index_status", {})
        smiles = await client.call_tool(
            "substructure_search", {"query": "c1ccccc1C(=O)O"}
        )
        smarts = await client.call_tool(
            "substructure_search",
            {"query": "[CX3](=O)[OX2H1]", "query_type": "smarts", "max_hits": 2},
        )
        bad = await client.call_tool(
            "substructure_search", {"query": "[C", "query_type": "smarts"}
        )
    assert status.structured_content["substructure"] is True
    assert smiles.is_error is False
    assert [h["cid"] for h in smiles.structured_content["hits"]] == [
        243,
        338,
        2244,
        4133,
    ]
    assert smiles.structured_content["exhaustive"] is True
    assert [h["cid"] for h in smarts.structured_content["hits"]] == [243, 338]
    assert smarts.structured_content["exhaustive"] is False
    assert bad.is_error is True
    assert "Invalid SMARTS" in bad.content[0].text


async def test_substructure_search_tool_without_pattern_index(built_home, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(built_home))
    async with Client(mcp) as client:
        result = await client.call_tool("substructure_search", {"query": "c1ccccc1"})
    assert result.is_error is True
    assert "index build --with-substructure" in result.content[0].text


async def test_pubchem_knowledge_tools(pubchem_knowledge_api):
    async with Client(mcp) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}
        sim = await client.call_tool(
            "pubchem_similarity_search", {"smiles": "CC(=O)Oc1ccccc1C(=O)O"}
        )
        syn = await client.call_tool("pubchem_synonyms", {"cid": 2244})
        ghs = await client.call_tool(
            "pubchem_annotations", {"cid": 2244, "heading": "GHS Classification"}
        )
        bio = await client.call_tool("pubchem_bioactivity", {"cid": 2244})
        xrefs = await client.call_tool("pubchem_xrefs", {"cid": 2244})
        desc = await client.call_tool("pubchem_description", {"cid": 2244})
        exact = await client.call_tool(
            "pubchem_exact_match", {"smiles": "CC(=O)Oc1ccccc1C(=O)O"}
        )
        sub = await client.call_tool(
            "pubchem_substructure_search", {"smiles": "c1ccccc1O"}
        )
    headings = tools["pubchem_annotations"].input_schema["properties"]["heading"]
    assert "GHS Classification" in str(headings)
    assert sim.structured_content["hits"][0]["cid"] == 2244
    assert syn.structured_content["cas_numbers"] == ["50-78-2"]
    assert "European Chemicals Agency (ECHA)" in ghs.structured_content["sources"]
    assert bio.structured_content["n_assays"] == 6
    assert xrefs.structured_content["counts"]["PubMedID"] == 2
    assert desc.structured_content["title"] == "Aspirin"
    assert exact.structured_content["hits"][0]["cid"] == 2244
    assert len(sub.structured_content["hits"]) == 2


async def test_annotation_and_timeout_errors_are_actionable(pubchem_knowledge_api):
    import httpx

    pubchem_knowledge_api["substructure"].mock(
        return_value=httpx.Response(504, text="Gateway Timeout")
    )
    async with Client(mcp) as client:
        missing = await client.call_tool(
            "pubchem_annotations", {"cid": 2244, "heading": "Toxicity"}
        )
        slow = await client.call_tool("pubchem_substructure_search", {"smiles": "C"})
    assert missing.is_error is True
    assert "GHS Classification" in missing.content[0].text  # lists what exists
    assert slow.is_error is True
    assert "more specific query" in slow.content[0].text


async def test_lookup_timeouts_suggest_a_retry_not_a_narrower_query(
    pubchem_knowledge_api,
):
    import httpx

    pubchem_knowledge_api["xrefs_ids"].mock(
        return_value=httpx.Response(504, text="Gateway Timeout")
    )
    async with Client(mcp) as client:
        slow = await client.call_tool("pubchem_xrefs", {"cid": 2244})
    assert slow.is_error is True
    assert "Retry in a moment" in slow.content[0].text
    assert "max_records" not in slow.content[0].text


async def test_xrefs_tool_reports_patents_it_could_not_get(pubchem_knowledge_api):
    import httpx

    pubchem_knowledge_api["xrefs_patents"].mock(
        return_value=httpx.Response(504, text="Gateway Timeout")
    )
    async with Client(mcp) as client:
        result = await client.call_tool("pubchem_xrefs", {"cid": 2244})
    assert result.is_error is False
    assert result.structured_content["unavailable"] == ["PatentID"]
    assert "PatentID" not in result.structured_content["counts"]


async def test_draw_molecules_returns_png_image():
    import base64

    async with Client(mcp) as client:
        result = await client.call_tool(
            "draw_molecules",
            {
                "smiles": ["CC(=O)Oc1ccccc1C(=O)O", "COC(=O)c1ccccc1O"],
                "legends": ["aspirin", "methyl salicylate"],
                "highlight": "diff",
                "reference_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            },
        )
        bad = await client.call_tool(
            "draw_molecules", {"smiles": ["CCO"], "highlight": "diff"}
        )
    assert result.is_error is False
    (content,) = result.content
    assert content.type == "image"
    assert content.mime_type == "image/png"
    assert base64.b64decode(content.data).startswith(b"\x89PNG")
    assert bad.is_error is True
    assert "reference" in bad.content[0].text


async def test_medchem_profile_and_scaffold_grouping(built_home, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(built_home))
    async with Client(mcp) as client:
        triage = await client.call_tool(
            "medchem_profile", {"smiles": "O=C1C(=Cc2ccccc2)SC(=S)N1"}
        )
        grouped = await client.call_tool(
            "similarity_search",
            {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "threshold": 0.2,
                "group_by_scaffold": True,
            },
        )
    assert triage.structured_content["pains_free"] is False
    scaffolds = grouped.structured_content["scaffolds"]
    assert scaffolds[0]["scaffold"] == "c1ccccc1"
    assert 2244 in scaffolds[0]["cids"]


# --- Prompts ------------------------------------------------------------------


async def test_prompts_are_listed_with_titles_and_single_token_arguments():
    async with Client(mcp) as client:
        prompts = {p.name: p for p in (await client.list_prompts()).prompts}
    assert set(prompts) == {"compound_dossier", "find_analogues", "triage_structures"}
    assert prompts["compound_dossier"].title == "Compound dossier"
    args = {p: [(a.name, a.required) for a in prompts[p].arguments] for p in prompts}
    assert args == {
        "compound_dossier": [("compound", True)],
        "find_analogues": [("compound", True), ("threshold", False)],
        "triage_structures": [("structures", False)],
    }


async def text_of(name, arguments):
    async with Client(mcp) as client:
        result = await client.get_prompt(name, arguments)
    assert [m.role for m in result.messages] == ["user"]
    return result.messages[0].content.text


async def test_compound_dossier_prompt():
    text = await text_of("compound_dossier", {"compound": "aspirin"})
    assert text.startswith("Write a concise dossier on aspirin")
    for tool in ("pubchem_lookup", "pubchem_annotations", "pubchem_bioactivity"):
        assert tool in text
    assert '"GHS Classification"' in text


async def test_find_analogues_prompt():
    text = await text_of(
        "find_analogues", {"compound": "ibuprofen", "threshold": "0.70"}
    )
    assert "similarity_search with threshold=0.7 and group_by_scaffold=true" in text
    assert 'search_scope is "demo"' in text
    # PubChem's search sends the structure there: offered, not run unasked.
    assert "offer to run pubchem_similarity_search" in text
    assert "ask first" in text
    assert "also run" not in text
    default = await text_of("find_analogues", {"compound": "ibuprofen"})
    assert "threshold=0.6" in default


@pytest.mark.parametrize("threshold", ["abc", "0", "1.5"])
async def test_find_analogues_rejects_a_bad_threshold(threshold):
    from mcp.shared.exceptions import MCPError

    async with Client(mcp) as client:
        with pytest.raises(MCPError, match="threshold must be a number"):
            await client.get_prompt(
                "find_analogues", {"compound": "aspirin", "threshold": threshold}
            )


async def test_triage_structures_prompt():
    listed = await text_of("triage_structures", {"structures": "CCO,c1ccccc1,"})
    assert "these structures:\n- CCO\n- c1ccccc1\n" in listed
    assert "medchem_profile" in listed
    empty = await text_of("triage_structures", {})
    assert "ask me for them as SMILES" in empty
    # The assistant has to pass on what standardizing changed.
    assert "say what it removed" in listed
    assert "pass its note on" in listed


# --- Resources ----------------------------------------------------------------


async def test_resources_are_listed():
    async with Client(mcp) as client:
        templates = (await client.list_resource_templates()).resource_templates
        resources = (await client.list_resources()).resources
    assert [(t.uri_template, t.mime_type) for t in templates] == [
        ("pubchem://compound/{cid}", "application/json")
    ]
    assert [str(r.uri) for r in resources] == ["pubchem://index/status"]


async def test_compound_resource(pubchem_api):
    import json

    async with Client(mcp) as client:
        result = await client.read_resource("pubchem://compound/2244")
    record = json.loads(result.contents[0].text)
    assert (record["cid"], record["title"]) == (2244, "Aspirin")
    assert result.contents[0].mime_type == "application/json"


async def test_compound_resource_reads_a_zero_padded_cid(pubchem_api):
    """The length limit used to count the padding: 11 or more digits were refused."""
    import json

    async with Client(mcp) as client:
        result = await client.read_resource("pubchem://compound/00000002244")
    assert json.loads(result.contents[0].text)["cid"] == 2244


@pytest.mark.parametrize(
    "cid", ["999999999", "abc", "0", "%C2%B2", "%D9%A3", "9" * 400]
)  # the two escapes are superscript two and an Arabic-Indic three
async def test_missing_compound_resource_is_not_found(pubchem_api, cid):
    from mcp.shared.exceptions import MCPError

    from tests.conftest import fault

    pubchem_api.get(path__regex=r"/compound/cid/999999999/property/.*").respond(
        404, json=fault("PUGREST.NotFound", "No CID found")
    )
    async with Client(mcp) as client:
        with pytest.raises(MCPError) as err:
            await client.read_resource(f"pubchem://compound/{cid}")
    assert err.value.error.code == -32602
    assert "No CID found" in str(err.value) or "Not a PubChem CID" in str(err.value)


async def test_index_status_resource(tmp_path, monkeypatch):
    import json

    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(tmp_path))
    async with Client(mcp) as client:
        result = await client.read_resource("pubchem://index/status")
    status = json.loads(result.contents[0].text)
    assert (status["state"], status["home"]) == ("not_built", str(tmp_path))


async def test_resource_failures_name_the_cause(monkeypatch):
    from mcp.shared.exceptions import MCPError

    from pubchem_finder import service

    def broken(*args, **kwargs):
        raise OSError(5, "Input/output error")

    monkeypatch.setattr(service, "index_status", broken)
    async with Client(mcp) as client:
        with pytest.raises(MCPError, match=r"Unexpected error \(OSError"):
            await client.read_resource("pubchem://index/status")


async def test_large_result_tools_raise_claude_codes_size_limit():
    async with Client(mcp) as client:
        tools = (await client.list_tools()).tools
    raised = {t.name: t.meta["anthropic/maxResultSizeChars"] for t in tools if t.meta}
    assert raised == {
        "similarity_search": 100_000,
        "substructure_search": 100_000,
        "pubchem_annotations": 100_000,
    }


async def test_oversized_smiles_get_a_tool_error_not_a_crash():
    # Hardening audit C1: 'C' * 80,000 overflowed RDKit's recursive SMILES writer
    # and killed the server process, dropping every session.
    async with Client(mcp) as client:
        result = await client.call_tool("compound_properties", {"smiles": "C" * 80_000})
        after = await client.call_tool("compound_properties", {"smiles": "CCO"})
    assert result.is_error is True
    assert "the limit is 10,000" in result.content[0].text
    assert after.is_error is False
