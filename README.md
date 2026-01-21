# Amber MCP

> Model Context Protocol server for Amber molecular dynamics simulation and analysis workflows

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Local Usage (Scripts)](#local-usage-scripts)
- [MCP Server Installation](#mcp-server-installation)
- [Using with Claude Code](#using-with-claude-code)
- [Using with Gemini CLI](#using-with-gemini-cli)
- [Available Tools](#available-tools)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

## Overview

This MCP server provides computational chemistry and molecular dynamics simulation capabilities based on AmberTools25 and Amber24. It offers both standalone scripts and MCP tools for analyzing MD simulations.

### Features
- **Single protein MD simualations**: Run MD simulations for a AlphaFold output PDF file

## Installation

### Create Environment
Please follow the environment setup procedure from `reports/step3_environment.md`. An example workflow is shown below.

```bash
# Navigate to the MCP directory
cd tool-mcps/amber_mcp
./quick_setup.sh
```
## Local Usage (Scripts)

You can use the scripts directly without MCP for local processing.

#### MD Analysis

```bash
# Activate environment
mamba activate ./env

# Run script
python scripts/single_protein_simulation.py example/input/1l2y.pdb -t 10 -o results/md_1l2y
```

## MCP Server Installation

### Option 1: Using fastmcp (Recommended)

```bash
# Activate environment and install MCP server for Claude Code
fastmcp install claude-code src/server.py --name amber_mcp --python  ./env/bin/python
```

## Using with Claude Code

After installing the MCP server, you can use it directly in Claude Code.

### Quick Start

```bash
# Start Claude Code
claude
```

### Example Prompts

#### Tool Discovery
```
What tools are available from amber_mcp?
```

#### Basic MD Analysis
```
Please run the MD analysis for @example/input/1l2y.pdb using cuda:0.
```

## License

Based on AmberTools25/Amber24 molecular dynamics package

## Credits

Based on [AmberTools25/Amber24](https://ambermd.org/) molecular dynamics simulation software