"""The CLI and MCP server can start without RDKit's drawing libraries."""

import subprocess
import sys


def test_entry_points_without_rdkit_drawing():
    script = """
import asyncio
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.modules["rdkit.Chem.Draw"] = None
from mcp.client.client import Client
from typer.testing import CliRunner
from pubchem_finder.cli import app
from pubchem_finder.mcp_server import mcp

runner = CliRunner()
version = runner.invoke(app, ["--version"])
assert version.exit_code == 0, version.output
with TemporaryDirectory() as directory:
    output = Path(directory) / "molecule.png"
    draw = runner.invoke(app, ["draw", "CCO", "-o", str(output)])
    assert draw.exit_code == 1, draw.output
    assert "RDKit drawing is unavailable" in draw.output
    if sys.platform == "linux":
        assert "apt-get install libxrender1 libxext6 libexpat1" in draw.output
    assert not output.exists()

async def check_tool():
    async with Client(mcp) as client:
        response = await client.call_tool("draw_molecules", {"smiles": ["CCO"]})
    assert response.is_error, response
    (error,) = response.content
    assert error.type == "text", error
    assert "RDKit drawing is unavailable:" in error.text, error.text
    if sys.platform == "linux":
        assert "apt-get install libxrender1 libxext6 libexpat1" in error.text, error.text

asyncio.run(check_tool())
"""
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, timeout=90
    )
    assert result.returncode == 0, result.stderr
