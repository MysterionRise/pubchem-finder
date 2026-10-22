# Using pubchem-finder with MCP clients

`pubchem-finder mcp` is a standard [Model Context Protocol](https://modelcontextprotocol.io)
server over stdio: the client starts it as a local process and talks to it through stdin and
stdout. Any client that can launch a local (stdio) MCP server can use it, not just Claude. CI
checks the server with the official MCP Python client on Linux, macOS and Windows.

This page collects tested setups. Claude Code and Claude Desktop are in the
[README](https://github.com/MysterionRise/pubchem-finder#use-it-from-claude); for the
others, [help is welcome](#add-your-client).

## What every client needs

| Setting | Value |
|---|---|
| Transport | stdio (a local command) |
| Command | `uvx` (install [uv](https://docs.astral.sh/uv/getting-started/installation/) first) |
| Arguments | `pubchem-finder`, `mcp` |
| Environment (optional) | `PUBCHEM_FINDER_HOME=/path` to keep the local index somewhere other than the default data directory |

Before adding it to a client, run `uvx pubchem-finder --version` once in a terminal. The
first run downloads Python packages (RDKit, PyArrow, DuckDB...) and can take a minute,
longer than some clients wait for a server to start; later starts take about a second.

If a desktop app can't find `uvx`, use its full path as the command: `which uvx` on macOS
and Linux (usually `~/.local/bin/uvx`), `where uvx` on Windows. Apps started from the dock
or Start menu often don't see your shell's `PATH`. If you'd rather not use uv,
`pip install pubchem-finder` and use `pubchem-finder` as the command with `mcp` as the only
argument.

Most clients take the JSON shape Claude Desktop uses:

```json
{
  "mcpServers": {
    "pubchem-finder": {
      "command": "uvx",
      "args": ["pubchem-finder", "mcp"]
    }
  }
}
```

A few use their own shape: VS Code's `mcp.json` has `servers` with `"type": "stdio"`,
Codex CLI uses TOML (`[mcp_servers.pubchem-finder]`), and Zed uses `context_servers` in its
settings. The sections below give the exact config for each client someone has tested.

## Check that it works

1. The client shows the `pubchem-finder` server with **17 tools**.
2. Ask *"What is the molecular weight of caffeine?"* It should call `pubchem_lookup` and
   answer 194.19 g/mol (CID 2519).
3. Ask *"Draw aspirin and ibuprofen side by side."* It should call `draw_molecules`, which
   returns a PNG. Note whether the client shows the picture: not every client renders
   images from tools.
4. With a local index (`pubchem-finder index build --limit 1`, a few minutes), ask
   *"Find close analogues of aspirin in the local index."* It should call
   `similarity_search` and list CID 2244 first.

## Clients

| Client | Status | Setup |
|---|---|---|
| Claude Code | Documented | [README](https://github.com/MysterionRise/pubchem-finder#use-it-from-claude) |
| Claude Desktop | Documented | [README](https://github.com/MysterionRise/pubchem-finder#use-it-from-claude) |
| Cursor | Not tried yet | Help wanted: [#24](https://github.com/MysterionRise/pubchem-finder/issues/24) |
| VS Code (GitHub Copilot agent mode) | Tested (VS Code 1.142.0, Windows 11) | [docs/clients.md](#vs-code-github-copilot)  |
| Windsurf | Not tried yet | Help wanted: [#26](https://github.com/MysterionRise/pubchem-finder/issues/26) |
| Zed | Not tried yet | Help wanted: [#27](https://github.com/MysterionRise/pubchem-finder/issues/27) |
| JetBrains IDEs (AI Assistant) | Not tried yet | Help wanted: [#28](https://github.com/MysterionRise/pubchem-finder/issues/28) |
| Google Antigravity / Gemini CLI | Not tried yet | Help wanted: [#29](https://github.com/MysterionRise/pubchem-finder/issues/29) |
| OpenAI Codex CLI | Not tried yet | Help wanted: [#30](https://github.com/MysterionRise/pubchem-finder/issues/30) |
| Cline | Not tried yet | Help wanted: [#31](https://github.com/MysterionRise/pubchem-finder/issues/31) |
| LM Studio (local models) | Not tried yet | Help wanted: [#32](https://github.com/MysterionRise/pubchem-finder/issues/32) |

Using a client that isn't listed? Open an issue or a pull request in the same format.

Web-based clients that only connect to remote servers (for example ChatGPT connectors)
need an HTTP transport, which pubchem-finder doesn't offer yet ([#33](https://github.com/MysterionRise/pubchem-finder/issues/33)).
All of this is tracked in [#23](https://github.com/MysterionRise/pubchem-finder/issues/23).

## Add your client

1. Set it up with the settings above and go through [Check that it works](#check-that-it-works).
2. Add a section below, in this format:

   ````markdown
   ### <Client name>

   Tested with <client> <version> on <OS and version>, pubchem-finder <version>.

   Config file: <where it lives, or the menu that opens it>.

   ```json
   <the exact config you used>
   ```

   Pictures from `draw_molecules`: shown inline / shown as an attachment / not shown.

   Notes: anything that tripped you up (a setting to enable, a timeout, a Windows quirk).
   ````

3. In the table, set the status to "Tested" with the version and OS, and link your
   section.
4. Open a pull request, and mention the issue it closes. A screenshot of a tool call is
   welcome.

The checks in CI don't cover client setups, so a short, exact section from someone who has
run it is what helps the next person most.

### VS Code (GitHub Copilot)

Tested with VS Code 1.142.0 and GitHub Copilot Chat 0.41 on Windows 11 26H2, pubchem-finder 2.0.1.

Config file: `.vscode/mcp.json` in the workspace (created by the **MCP: Add Server** command, or by hand).

```json
{
  "servers": {
    "pubchem-finder": {
      "type": "stdio",
      "command": "uvx",
      "args": ["pubchem-finder", "mcp"]
    }
  }
}
```

Pictures from `draw_molecules`: not shown. The model replied "[Generated image shown above]", but no image appeared in the chat.

Notes:
- Starting the server from `mcp.json` right after creating it failed with "Cannot start an MCP server before the session runtime has been created. Send a message first." Sending any message in an agent-mode chat first fixed it.
