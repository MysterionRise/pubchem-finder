"""``python -m pubchem_finder``: the same CLI as the ``pubchem-finder`` command."""

from .cli import run

# Guarded: multiprocessing's spawn start method re-imports the main module in
# each worker process (index builds use it), which must not start the CLI again.
if __name__ == "__main__":
    run()
