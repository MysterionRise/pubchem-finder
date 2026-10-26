## What and why

<!-- What this changes and the problem it solves. Link the issue if there is one. -->

## How it was tested

<!-- `make check` passes; anything else you ran (live tests, smoke script, a real index). -->

## Client setup PRs only

<!-- For a docs/clients.md section (issues #24-#32); delete this part otherwise.
     Paste what the client showed, not a summary: the reviewer reproduces these. -->

- Client, version and OS; pubchem-finder version:
- Tools the client lists for the server (expected 17):
- "What is the molecular weight of caffeine?": the tool call and the answer as shown
  (expected `pubchem_lookup`, 194.19 g/mol, CID 2519):
- "Draw aspirin and ibuprofen side by side.": was `draw_molecules` called, and did the
  picture show?
- With a local index, "Find close analogues of aspirin in the local index.": the first CID
  shown (expected 2244):
- Did the client's tool-call panel show each call, or did the model answer on its own?
  (Models sometimes skip the tool. Saying so is useful.)
- A screenshot of one tool call:
