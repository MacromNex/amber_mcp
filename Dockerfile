# syntax=docker/dockerfile:1
# =============================================================================
# Dockerfile for Amber MCP Server
#
# Builds AmberTools25 + PMEMD24 from source and sets up the FastMCP server
# for molecular dynamics simulations.
#
# Prerequisites:
#   The repo/ambertools25_src directory must be present in the build context
#   with the AmberTools25 source code (including PMEMD24 merged).
#
# Build (CPU-only):
#   docker build -t amber_mcp .
#
# Build (with CUDA - requires NVIDIA base image swap):
#   docker build --build-arg ENABLE_CUDA=true -t amber_mcp:cuda .
#
# Run:
#   docker run -i amber_mcp
# =============================================================================

FROM condaforge/miniforge3:latest AS builder

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install OS-level build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    flex \
    bison \
    patch \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---- Dependency layer (cached across rebuilds) ----
# Install conda build dependencies with a persistent package cache
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba install -y -n base \
        python=3.11 \
        cmake \
        gcc \
        gxx \
        gfortran \
        openmpi \
        openmpi-mpicc \
        openmpi-mpicxx \
        openmpi-mpifort \
        netcdf-fortran \
        netcdf4 \
        boost \
        fftw \
        arpack \
        flex \
        bison \
        patch \
        make

# Install Python MCP dependencies with a persistent pip cache
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install fastmcp loguru 'numpy<2.0'

# ---- Build layer ----
# Copy the build script and source code
COPY quick_setup.sh ./
COPY repo/ambertools25_src ./repo/ambertools25_src/
RUN chmod -R a+r /app/repo/

# Build AmberTools25 (CPU-only, skip env setup since we did it above)
RUN bash quick_setup.sh \
        --source-dir /app/repo/ambertools25_src \
        --install-prefix /app/env \
        --no-cuda \
        --no-setup-env \
        --clean

# ---------------------------------------------------------------------------
# Runtime stage - keep only what's needed
# ---------------------------------------------------------------------------
FROM condaforge/miniforge3:latest AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the compiled Amber environment from builder
COPY --from=builder /app/env /app/env

# Copy conda packages (runtime libraries: openmpi, netcdf, fftw, etc.)
COPY --from=builder /opt/conda /opt/conda

# Install Python MCP dependencies in runtime (uses builder's conda/python)
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install fastmcp loguru 'numpy<2.0'

# Copy the MCP server source
COPY src/ ./src/
RUN chmod -R a+r /app/src/

# Create working directories
RUN mkdir -p results jobs tmp

# Set environment variables for Amber
ENV AMBERHOME=/app/env
ENV PMEMDHOME=/app/env
ENV PATH="/app/env/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/env/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH=/app
ENV OPAL_PREFIX=/opt/conda

# Source amber.sh environment on shell entry
RUN echo 'test -f /app/env/amber.sh && source /app/env/amber.sh' >> /etc/bash.bashrc

CMD ["python", "src/server.py"]
