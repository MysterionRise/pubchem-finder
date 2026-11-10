# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "langchain-mcp-adapters>=0.3,<0.4",
#     "langgraph>=1,<2",
#     "langchain-anthropic>=1,<2",
# ]
# ///
"""Give a LangChain / LangGraph agent pubchem-finder's tools through MCP.

    ANTHROPIC_API_KEY=sk-ant-... uv run examples/langchain_agent.py
    PUBCHEM_FINDER_COMMAND=.venv/bin/pubchem-finder uv run examples/langchain_agent.py

Needs an ANTHROPIC_API_KEY (the agent runs on Claude). The server itself needs no
key. `langchain-mcp-adapters` starts `uvx --python 3.14 pubchem-finder mcp` over stdio
and turns its 17 tools into LangChain tools; `create_react_agent` lets the model call
them. The question works on a fresh install: until a local index is built
(`uvx --python 3.14 pubchem-finder index build --limit 1`), `similarity_search` answers
from the built-in set of FDA-approved drugs and says so in its `note`.

The plain MCP client without a framework is examples/mcp_client.py.
"""

import asyncio
import os
import shlex

from langchain_anthropic import ChatAnthropic
from langchain_mcp_adapters.client import MultiServerMCPClient
from langgraph.prebuilt import create_react_agent

MODEL = os.environ.get("ANTHROPIC_MODEL", "claude-sonnet-4-5")
QUESTION = (
    "Find close analogues of aspirin in the local index and summarise what changed."
)


async def main() -> None:
    command = shlex.split(
        os.environ.get("PUBCHEM_FINDER_COMMAND", "uvx --python 3.14 pubchem-finder"),
        posix=os.name != "nt",
    )
    client = MultiServerMCPClient(
        {
            "pubchem-finder": {
                "command": command[0],
                "args": [*command[1:], "mcp"],
                "transport": "stdio",
                "env": dict(os.environ),
            }
        }
    )
    tools = await client.get_tools()
    print(f"{len(tools)} tools: " + ", ".join(t.name for t in tools))

    agent = create_react_agent(ChatAnthropic(model=MODEL), tools)
    result = await agent.ainvoke({"messages": [("user", QUESTION)]})
    print(result["messages"][-1].content)


asyncio.run(main())
