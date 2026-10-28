# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=2.3,<3"]
# ///
"""Call pubchem-finder's MCP tools from plain Python, with no agent framework.

    uv run examples/mcp_client.py      # starts `uvx --python 3.14 pubchem-finder mcp`
    PUBCHEM_FINDER_COMMAND=.venv/bin/pubchem-finder uv run examples/mcp_client.py

Lists the tools, computes aspirin's properties (offline) and, if a local index exists
(`uvx --python 3.14 pubchem-finder index build --limit 1`), finds its closest analogues. No API key
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
        if status["state"] != "ready":
            print(f"no local index: {status['hint']}")
            return
        hits = await session.call_tool(
            "similarity_search", {"smiles": ASPIRIN, "threshold": 0.6, "top_k": 5}
        )
        print(f"closest analogues in {status['n_compounds']:,} indexed compounds:")
        for hit in hits.structured_content["hits"]:
            print(
                f"  CID {hit['cid']:>9}  {hit['similarity']:.3f}  {hit['iupac_name']}"
            )


asyncio.run(main())
