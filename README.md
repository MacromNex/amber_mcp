# Amber MCP

**Amber molecular dynamics simulation and analysis via Docker**

An MCP (Model Context Protocol) server for molecular dynamics simulations using AmberTools25/Amber24 with 1 core tool:
- **single_protein_md_simulation** -- Run full MD simulations for protein PDB structures with GPU acceleration

## Quick Start with Docker

### Approach 1: Pull Pre-built Image from GitHub

The fastest way to get started. A pre-built Docker image is automatically published to GitHub Container Registry on every release.

```bash
# Pull the latest image
docker pull ghcr.io/macromnex/amber_mcp:latest

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add amber_mcp -- docker run -i --rm --user `id -u`:`id -g` --gpus all --ipc=host -v `pwd`:`pwd` ghcr.io/macromnex/amber_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker with NVIDIA GPU support (nvidia-docker2 or Docker 19.03+)
- Claude Code installed

That's it! The Amber MCP server is now available in Claude Code.

---

### Approach 2: Build Docker Image Locally

Build the image yourself and install it into Claude Code. Useful for customization or offline environments.

```bash
# Clone the repository
git clone https://github.com/MacromNex/amber_mcp.git
cd amber_mcp

# Build the Docker image
docker build -t amber_mcp:latest .

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add amber_mcp -- docker run -i --rm --user `id -u`:`id -g` --gpus all --ipc=host -v `pwd`:`pwd` amber_mcp:latest
```

**About the Docker Flags:**
- `-i` -- Interactive mode for Claude Code
- `--rm` -- Automatically remove container after exit
- `` --user `id -u`:`id -g` `` -- Runs the container as your current user
- `--gpus all` -- Grants access to all available GPUs for MD simulations
- `--ipc=host` -- Uses host IPC namespace for shared memory access
- `-v` -- Mounts your project directory

---

## Verify Installation

```bash
claude mcp list
# You should see 'amber_mcp' in the output
```

In Claude Code, you can now use:
- `single_protein_md_simulation` -- Run MD simulations for protein structures

---

## Next Steps

- **Detailed documentation**: See [detail.md](detail.md) for comprehensive guides including local installation, script usage, and environment setup

---

## Usage Examples

### Run an MD Simulation

```
Run a 10 ns MD simulation for @example/input/1l2y.pdb using cuda:0, save results to results/md_1l2y
```

### Analyze a PDB Structure

```
What tools are available from amber_mcp? Then run an MD simulation for my protein structure at /path/to/protein.pdb
```

### AlphaFold Output Processing

```
Run MD simulation for my AlphaFold predicted structure @predictions/model_1.pdb with 5 ns simulation time
```

---

## Troubleshooting

### Docker Issues

**Problem:** Container cannot access GPU
```bash
# Verify NVIDIA Docker runtime is installed
docker run --rm --gpus all nvidia/cuda:12.0-base nvidia-smi

# If this fails, install nvidia-docker2:
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html
```

**Problem:** Permission denied on output files
```bash
# Ensure you're using the --user flag
docker run -i --rm --user `id -u`:`id -g` --gpus all --ipc=host -v `pwd`:`pwd` ghcr.io/macromnex/amber_mcp:latest
```

**Problem:** MCP server not responding
```bash
# Check if the container runs correctly
docker run --rm ghcr.io/macromnex/amber_mcp:latest --help

# Re-register with Claude Code
claude mcp remove amber_mcp
claude mcp add amber_mcp -- docker run -i --rm --user `id -u`:`id -g` --gpus all --ipc=host -v `pwd`:`pwd` ghcr.io/macromnex/amber_mcp:latest
```

**Problem:** File paths not found inside container
```bash
# Make sure to run Claude Code from the directory containing your files
# The -v `pwd`:`pwd` flag mounts only the current directory
cd /path/to/your/project
claude
```

---

## License

Based on AmberTools25/Amber24 molecular dynamics package.

## Credits

Based on [AmberTools25/Amber24](https://ambermd.org/) molecular dynamics simulation software.
