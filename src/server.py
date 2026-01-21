"""
Model Context Protocol (MCP) Server for Amber MD Simulations.

This MCP server provides tools for running molecular dynamics simulations
using the Amber software suite.

Tools:
1. amber_run_protein_md: Run a complete MD simulation workflow for a single protein
   - System preparation (solvation, ionization)
   - Energy minimization
   - Heating (NVT)
   - Equilibration (NPT)
   - Production MD (NPT)

2. amber_prepare_system: Prepare a protein system with tleap only
   - Solvation in water box
   - Ion addition for neutralization
   - Outputs topology and coordinate files

3. amber_generate_input_files: Generate Amber input files without running simulations
   - Creates all input files for inspection/customization
   - Useful for dry-run preparation

Reference: https://ambermd.org/
"""

import sys
from pathlib import Path

from fastmcp import FastMCP

# Add current directory to path for imports
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Import tool modules
from tools.simulation import simulation_mcp

# Server definition and mounting
mcp = FastMCP(name="amber_mcp")
mcp.mount(simulation_mcp)

if __name__ == "__main__":
    mcp.run()
