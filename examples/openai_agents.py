# /// script
# requires-python = ">=3.10"
# dependencies = ["openai-agents"]
# ///
"""Give an OpenAI Agents SDK agent pubchem-finder's tools over MCP.

Requires OPENAI_API_KEY in the environment.

    uv run examples/openai_agents.py
"""

import asyncio

from agents import Agent, Runner
from agents.mcp import MCPServerStdio

QUESTION = "Find close analogues of aspirin and summarise what changed."


async def main() -> None:
    async with MCPServerStdio(
        name="pubchem-finder",
        params={"command": "uvx", "args": ["--python", "3.14", "pubchem-finder", "mcp"]},
        # The first uvx run downloads packages; the default 5 s timeout is too short.
        client_session_timeout_seconds=120,
    ) as server:
        agent = Agent(
            name="Chemist",
            instructions=(
                "Answer with pubchem-finder's tools. If a tool result includes a "
                "note (e.g. results came from the built-in approved-drug set), "
                "mention it."
            ),
            mcp_servers=[server],
        )
        result = await Runner.run(agent, QUESTION)
        print(result.final_output)


if __name__ == "__main__":
    asyncio.run(main())