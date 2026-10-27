---
name: review-client-setup
description: Review a "tested setup" pull request for docs/clients.md (the #23 family of issues, one per MCP client). Verifies the screenshot, the versions, the config shape and the tool outputs the contributor reports against what pubchem-finder really returns, then drafts the review for the maintainer to post. Use when asked to review or validate a client-setup PR, or given a PR number that touches docs/clients.md.
---

# Review a client-setup PR

Usage: `/review-client-setup <PR number>`. The PR number is in `$ARGUMENTS`.

These PRs add a `### <Client>` section to docs/clients.md and flip a row in its table to
"Tested". No code changes, so CI proves nothing about them; the review has to check that
the setup was really run and that what it reports is what the server returns. Models
sometimes answer from memory without calling a tool, and contributors (often first-timers)
don't always notice. Be exact and be kind: every one of these PRs is welcome.

## 1. Read the PR

```bash
gh pr view $PR --json title,body,author,files,headRefName,commits,statusCheckRollup
gh pr diff $PR
gh pr list --author <login> --state all --repo MysterionRise/pubchem-finder   # first PR here?
```

- The diff should touch only `docs/clients.md` (a README table row is fine too). Anything
  else needs a reason.
- The new section must follow the template in docs/clients.md ("Add your client"):
  `Tested with <client> <version> on <OS>, pubchem-finder <version>.`, `Config file:`,
  a json/toml block, `Pictures from draw_molecules:`, `Notes:`. The table row must link
  the section with an in-page anchor (`[Setup](#client-name)`).
- Download every screenshot in the PR body (`https://github.com/user-attachments/assets/...`)
  into the scratchpad with `curl -sL -o` and look at it with Read. It should show the
  client's own tool-call UI, not just a chat answer.

## 2. Check that the versions exist and the config shape is current

| Client | Where to check |
|---|---|
| VS Code / Copilot | `gh api repos/microsoft/vscode/releases/latest`, `gh api repos/microsoft/vscode-copilot-chat/releases/latest` |
| Codex CLI | `gh api repos/openai/codex/releases/latest` |
| Gemini CLI | `gh api repos/google-gemini/gemini-cli/releases/latest` |
| Cline | `gh api repos/cline/cline/releases/latest` |
| Zed | `gh api repos/zed-industries/zed/releases/latest` |
| Cursor, Windsurf, JetBrains, LM Studio, OS versions | web search for the version and date |

A version newer than the latest release, or an OS build that doesn't exist, means the
numbers were made up; ask for the real ones. Then open the client's MCP docs (the issue
for that client links them) and confirm the config shape and file location the PR uses
are the current ones.

## 3. Reproduce the four checks with the project itself

The expected answers, from the code rather than memory:

```bash
cd /path/to/pubchem-finder
uv run python -c "from pubchem_finder import service; c = service.lookup('caffeine').compounds[0]; print(c.cid, c.molecular_weight)"
# 2519 194.19
uv run pubchem-finder search "CC(=O)Oc1ccccc1C(=O)O" --threshold 0.7 --top 5   # needs a local index: CID 2244 first
```

- `docs/media/demo-data.json` holds the real top-40 similarity hits for aspirin on
  PubChem's first chunk with their scores. A PR that quotes other CIDs or other scores for
  "close analogues of aspirin" did not get them from the tool.
- For any similarity number quoted for a specific pair, compute it:
  `service.compare(smiles_a, smiles_b).ecfp4_tanimoto` (look the CID up with
  `service.lookup(str(cid), by="cid")` to get its SMILES).
- The tool count is 17 (`uv run python scripts/smoke.py` prints it).

Anything that doesn't match is almost always "the model answered without calling the
tool". That is not a reason to reject the PR; it is a quirk worth a line in the section's
Notes, and the contributor should re-run the question watching the tool-call panel.

## 4. Draft the review

Write the review for the maintainer, in this shape, and print the exact command to post
it. Never post it yourself.

- One sentence of thanks that names what checked out (the screenshot, the config).
- A numbered list of requested changes, each with the exact replacement text.
- Which numbers in the PR body don't match the tool, with the real values and the
  reminder that "CID 2244 first" is the check.
- "Then I'll merge and close #<issue>."

```bash
gh pr review $PR --request-changes --body "..."    # or --approve
```

After the merge: close the client's issue, tick it in #23, add the CHANGELOG line under
Unreleased.

## Welcome reply for "I'd like to work on this"

When someone volunteers on one of the client issues, reply along these lines and assign
them:

> Great, it's yours. docs/clients.md has the settings, the four questions that check the
> setup and the format for your section. Two things that trip people up: run
> `uvx pubchem-finder --version` once in a terminal first (the first start downloads
> packages and can outlast a client's startup timeout), and watch the client's tool-call
> panel when you ask the questions, because models sometimes answer from memory without
> calling the tool; if yours does, say so in Notes. A screenshot of one tool call is
> welcome.
