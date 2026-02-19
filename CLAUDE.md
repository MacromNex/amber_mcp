# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Amber MCP is a Model Context Protocol (MCP) server wrapping AmberTools25/Amber24 for molecular dynamics simulations. It exposes three MCP tools (`amber_run_protein_md`, `amber_prepare_system`, `amber_generate_input_files`) that allow AI assistants to drive full protein MD workflows over stdio.

## Architecture

```
MCP Client (stdio)
    → src/server.py          FastMCP("amber_mcp") entrypoint, mounts simulation sub-server
    → src/tools/simulation.py Core logic: SimulationConfig, InputFileGenerator, SimulationRunner
    → env/bin/               Amber binaries (pmemd.cuda, pmemd, sander, tleap)
    → results/               Timestamped output directories
```

**`src/tools/simulation.py`** (~760 lines) contains everything: config validation, input file generation, subprocess execution of Amber binaries, and all three MCP tool functions. `SimulationRunner` sources `env/amber.sh` into subprocess environment and auto-detects the best MD engine (`pmemd.cuda` > `pmemd` > `sander`).

**`scripts/single_protein_simulation.py`** is a standalone CLI mirror of the MCP tools (same logic, argparse interface, no fastmcp dependency).

## Running the MCP Server

```bash
./env/bin/python src/server.py
```

Or install into Claude Code:
```bash
fastmcp install claude-code src/server.py --name amber_mcp --python ./env/bin/python
```

## Building Amber from Source

```bash
./quick_setup.sh                          # full build with CUDA
./quick_setup.sh --no-cuda                # CPU-only
./quick_setup.sh --cuda-path /usr/local/cuda --jobs 8
```

`quick_setup.sh` creates a conda env at `./env`, installs all build dependencies, applies source patches (CMake fixes, TFE GBSA Fortran declarations, numpy<2.0 for parmed), builds with CMake+make, and verifies key binaries.

## Docker

```bash
docker build -t amber_mcp .              # CPU-only multi-stage build
docker run -i amber_mcp                  # runs src/server.py on stdio
```

CI pushes to `ghcr.io` on every push to `main` and version tags via `.github/workflows/docker.yml`.

## Testing

No automated test suite exists. Manual testing is done via the example workflow:

```bash
source ./env/amber.sh
export OPAL_PREFIX=$(dirname $(which pmemd))/..
cd example && ./quick_test.sh             # fast integration check
cd example && ./run_all.sh                # full pipeline on 1l2y.pdb
```

## Key Implementation Details

- **Python deps**: `fastmcp`, `loguru`, `numpy<2.0` (pinned for parmed compatibility). No pyproject.toml; deps installed by `quick_setup.sh`.
- **Python version**: 3.11 (pinned in setup script and Dockerfile)
- **Production timestep**: 2 fs. Step count = `sim_time_ns * 500_000`.
- **Simulation stages**: tleap → minimization with restraints → minimization without restraints → NVT heating (50 ps) → NPT equilibration (500 ps) → NPT production.
- **Force fields**: ff14SB, ff19SB (default). **Water models**: tip3p, opc (default), tip4pew.
- **GPU selection**: `SimulationConfig` parses `gpu_device` (e.g. `"cuda:1"` or `"1"`) and sets `CUDA_VISIBLE_DEVICES`.
- **Environment activation**: `SimulationRunner` reads `env/amber.sh` and injects all exported variables into subprocess `env` dict — it does not use shell sourcing.
