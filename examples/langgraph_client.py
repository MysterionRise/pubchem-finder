# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "langchain-mcp-adapters==0.3.2",
#     "langchain-openai==1.8.5",
#     "langgraph==1.2.24",
# ]
# ///
"""Use pubchem-finder's MCP tools from a LangGraph agent.

Requires OPENAI_API_KEY. Run with:
    OPENAI_API_KEY=... uv run examples/langgraph_client.py

The MCP server runs through:
    uvx --python 3.14 pubchem-finder mcp
"""

import asyncio
import os
from typing import Any

from langchain_mcp_adapters.client import MultiServerMCPClient
from langchain_openai import ChatOpenAI
from langgraph.graph import START, MessagesState, StateGraph
from langgraph.prebuilt import ToolNode, tools_condition


async def main() -> None:
    if not os.environ.get("OPENAI_API_KEY"):
        raise SystemExit("Set OPENAI_API_KEY before running this example.")

    client = MultiServerMCPClient(
        {
            "pubchem-finder": {
                "command": "uvx",
                "args": ["--python", "3.14", "pubchem-finder", "mcp"],
                "transport": "stdio",
                "env": dict(os.environ),
            }
        }
    )
    tools = await client.get_tools()
    print(f"Loaded {len(tools)} pubchem-finder tools.")

    model = ChatOpenAI(model="gpt-4.1-mini", temperature=0).bind_tools(tools)

    async def call_model(state: MessagesState) -> dict[str, list[Any]]:
        response = await model.ainvoke(state["messages"])
        return {"messages": [response]}

    builder = StateGraph(MessagesState)
    builder.add_node("call_model", call_model)
    builder.add_node("tools", ToolNode(tools))
    builder.add_edge(START, "call_model")
    builder.add_conditional_edges("call_model", tools_condition)
    builder.add_edge("tools", "call_model")
    graph = builder.compile()

    result = await graph.ainvoke(
        {
            "messages": [
                {
                    "role": "user",
                    "content": (
                        "Find close analogues of aspirin using the local "
                        "similarity search. Summarize how their structures "
                        "differ from aspirin and include their similarity "
                        "scores. Explain any limitations of the search."
                    ),
                }
            ]
        }
    )
    print("\n" + str(result["messages"][-1].content))


if __name__ == "__main__":
    asyncio.run(main())
