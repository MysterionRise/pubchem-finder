# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=2.3,<3"]
# ///
"""Call pubchem-finder's MCP tools from plain Python, with no agent framework.

    uv run examples/mcp_client.py      # starts `uvx --python 3.14 pubchem-finder mcp`
    PUBCHEM_FINDER_COMMAND=.venv/bin/pubchem-finder uv run examples/mcp_client.py

Lists the tools, computes aspirin's properties (offline) and finds its closest
analogues: in the local index once one is built (`uvx --python 3.14 pubchem-finder index
build --limit 1 --with-substructure`), until then in the built-in demo set of
FDA-approved drugs. No API key
needed. Agent frameworks (OpenAI Agents SDK, LangChain, PydanticAI...) connect to the
same server with the same command; see
https://github.com/MysterionRise/pubchem-finder/issues/34.
"""

import asyncio
import os
import shlex

from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


async def main() -> None:
    command = shlex.split(
        os.environ.get("PUBCHEM_FINDER_COMMAND", "uvx --python 3.14 pubchem-finder"),
        posix=os.name != "nt",
    )
    params = StdioServerParameters(
        command=command[0], args=[*command[1:], "mcp"], env=dict(os.environ)
    )
    async with (
        stdio_client(params) as (read, write),
        ClientSession(read, write) as session,
    ):
        await session.initialize()
        tools = await session.list_tools()
        print(f"{len(tools.tools)} tools: " + ", ".join(t.name for t in tools.tools))

        props = (
            await session.call_tool("compound_properties", {"smiles": ASPIRIN})
        ).structured_content
        print(
            f"aspirin: MW {props['molecular_weight']}, logP {props['logp']}, "
            f"Lipinski ok: {props['lipinski_ok']}"
        )

        status = (await session.call_tool("index_status", {})).structured_content
        # Releases before 2.2 have no demo set and no search_scope.
        ready = "index" if status["state"] == "ready" else "none"
        if status.get("search_scope", ready) == "none":
            print(f"no local index: {status['hint']}")
            return
        found = (
            await session.call_tool(
                "similarity_search", {"smiles": ASPIRIN, "threshold": 0.3, "top_k": 5}
            )
        ).structured_content
        if found.get("note"):
            print(found["note"])
        print(f"closest analogues in {found['n_indexed']:,} indexed compounds:")
        for hit in found["hits"]:
            print(
                f"  CID {hit['cid']:>9}  {hit['similarity']:.3f}  {hit['iupac_name']}"
            )


asyncio.run(main())
