# Using pubchem-finder with MCP clients

`pubchem-finder mcp` is a standard [Model Context Protocol](https://modelcontextprotocol.io)
server over stdio: the client starts it as a local process and talks to it through stdin and
stdout. Any client that can launch a local (stdio) MCP server can use it, not just Claude.
For clients that connect to a URL instead, `pubchem-finder mcp --transport http` serves the
same tools over streamable HTTP ([Over HTTP](#over-http)). CI checks both transports with
the official MCP Python client on Linux, macOS and Windows.

This page collects tested setups. Claude Code and Claude Desktop are in the
[README](https://github.com/MysterionRise/pubchem-finder#use-it-from-claude); for the
others, [help is welcome](#add-your-client).

## What every client needs

| Setting | Value |
|---|---|
| Transport | stdio (a local command); or streamable HTTP, see [Over HTTP](#over-http) |
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
| VS Code (GitHub Copilot agent mode) | Not tried yet | Help wanted: [#25](https://github.com/MysterionRise/pubchem-finder/issues/25) |
| Windsurf | Not tried yet | Help wanted: [#26](https://github.com/MysterionRise/pubchem-finder/issues/26) |
| Zed | Not tried yet | Help wanted: [#27](https://github.com/MysterionRise/pubchem-finder/issues/27) |
| JetBrains IDEs (AI Assistant) | Not tried yet | Help wanted: [#28](https://github.com/MysterionRise/pubchem-finder/issues/28) |
| Google Antigravity / Gemini CLI | Not tried yet | Help wanted: [#29](https://github.com/MysterionRise/pubchem-finder/issues/29) |
| OpenAI Codex CLI | Not tried yet | Help wanted: [#30](https://github.com/MysterionRise/pubchem-finder/issues/30) |
| Cline | Not tried yet | Help wanted: [#31](https://github.com/MysterionRise/pubchem-finder/issues/31) |
| LM Studio (local models) | Not tried yet | Help wanted: [#32](https://github.com/MysterionRise/pubchem-finder/issues/32) |

Using a client that isn't listed? Open an issue or a pull request in the same format.

Clients that only connect to remote servers (web-based chat apps and their connectors)
need a URL: see [Over HTTP](#over-http). All of this is tracked in
[#23](https://github.com/MysterionRise/pubchem-finder/issues/23).

## Over HTTP

`pubchem-finder mcp --transport http` serves the same 17 tools over
[streamable HTTP](https://modelcontextprotocol.io/specification/2025-06-18/basic/transports#streamable-http)
instead of stdio:

```bash
pubchem-finder mcp --transport http                  # http://127.0.0.1:8000/mcp
pubchem-finder mcp --transport http --port 8765 --home /data/pubchem-finder
```

It prints the URL to stderr and answers `GET /health` with its version. The MCP endpoint
is `/mcp`. Use it when:

- the client connects to a URL rather than starting a process (web-based chat apps,
  agents running elsewhere);
- one machine holds the index for a team: everyone points their client at that URL and
  shares one copy of the index instead of building their own.

Any client that supports streamable HTTP can connect, for example Claude Code with
`claude mcp add --transport http pubchem-finder http://127.0.0.1:8000/mcp`.

**Check that it works:**

```bash
curl -s http://127.0.0.1:8000/health
# {"status":"ok","name":"pubchem-finder","version":"..."}
uv run python scripts/smoke.py --url http://127.0.0.1:8000/mcp   # every tool, from the repo
```

**Beyond your own machine.** The server binds to `127.0.0.1` by default and has no
authentication or TLS of its own. To share it, keep that default and put an HTTPS reverse
proxy with authentication in front of it (Caddy, nginx, a Cloudflare Tunnel...), then bind
to the interface the proxy reaches with `--host`. With Caddy, `/etc/caddy/Caddyfile`:

```
pubchem.example.com {
    basic_auth {
        alice <hash from `caddy hash-password`>
    }
    reverse_proxy 127.0.0.1:8000
}
```

Binding to anything but `127.0.0.1` prints a warning. For loopback binds the server
also rejects requests whose `Host` header is not local (DNS-rebinding protection); for
other binds it accepts any `Host`, so the proxy is the trust boundary. A Docker image and
an `--allowed-host` option are on the roadmap if shared servers see use.

**What a shared server shares:**

- **PubChem's rate limit.** Every user's `pubchem_*` request leaves from the server's IP,
  so they share PubChem's limit of 5 requests per second, and one response cache. The
  local-index tools (`similarity_search`, `substructure_search`) have no such limit, so a
  team server is most useful with a large index.
- **One index in memory.** The first search loads the index; later sessions reuse it.
  A full-PubChem index takes tens of seconds to load (and about 17 GB of RAM), so expect
  the first query after a restart to be slow.

Out of scope for now: built-in OAuth, and a Docker image
([#33](https://github.com/MysterionRise/pubchem-finder/issues/33)).

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
4. Open a pull request and fill in its "Client setup PRs" part, and mention the issue it
   closes. A screenshot of a tool call is welcome.

The checks in CI don't cover client setups, so a short, exact section from someone who has
run it is what helps the next person most.

**What reviewers check:** that the config matches the client's current docs, that the
versions exist, and that the answers in the PR are the tool's: the caffeine lookup is
CID 2519 at 194.19 g/mol, and the aspirin search lists CID 2244 first. Models sometimes
answer from memory without calling a tool; if yours did, watch the client's tool-call
panel, ask again, and say so in Notes. That is useful information for the next person.
