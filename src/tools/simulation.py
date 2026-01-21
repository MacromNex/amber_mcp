"""
Amber MD Simulation Tools for MCP Server.

This module provides MCP tools for running molecular dynamics simulations
using the Amber software suite.

Tools:
1. amber_run_protein_md: Run a complete MD simulation workflow for a single protein
2. amber_prepare_system: Prepare a protein system with tleap (solvation, ionization)
3. amber_generate_input_files: Generate Amber input files without running simulations
"""

import os
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Annotated, Optional

from fastmcp import FastMCP
from loguru import logger

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results"
AMBER_ENV_DIR = PROJECT_ROOT / "env"

# Ensure output directory exists
DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# MCP server instance
simulation_mcp = FastMCP(name="simulation")

# Force field mappings
FORCEFIELD_MAP = {
    "ff14SB": "leaprc.protein.ff14SB",
    "ff14sb": "leaprc.protein.ff14SB",
    "ff19SB": "leaprc.protein.ff19SB",
    "ff19sb": "leaprc.protein.ff19SB",
}

# Water model mappings
WATER_MODEL_MAP = {
    "tip3p": ("leaprc.water.tip3p", "TIP3PBOX"),
    "TIP3P": ("leaprc.water.tip3p", "TIP3PBOX"),
    "opc": ("leaprc.water.opc", "OPCBOX"),
    "OPC": ("leaprc.water.opc", "OPCBOX"),
    "tip4pew": ("leaprc.water.tip4pew", "TIP4PEWBOX"),
    "TIP4PEW": ("leaprc.water.tip4pew", "TIP4PEWBOX"),
}


@dataclass
class SimulationConfig:
    """Configuration for MD simulation parameters."""

    pdb_file: Path
    job_name: str = ""
    sim_time_ns: float = 10.0
    temperature: float = 300.0
    box_buffer: float = 12.0
    salt_conc: float = 0.15
    forcefield: str = "ff19SB"
    water_model: str = "opc"
    use_gpu: bool = True
    gpu_device: str = "cuda:0"  # GPU device selection (e.g., "cuda:0", "cuda:1")
    dry_run: bool = False
    output_dir: Path = field(default_factory=Path)

    # Derived fields
    ff_source: str = field(init=False, default="")
    water_source: str = field(init=False, default="")
    water_box: str = field(init=False, default="")
    cuda_device_id: str = field(init=False, default="0")  # Extracted CUDA device ID

    def __post_init__(self) -> None:
        """Validate inputs and set derived fields."""
        if isinstance(self.pdb_file, str):
            self.pdb_file = Path(self.pdb_file)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_file}")

        self.pdb_file = self.pdb_file.resolve()

        if not self.job_name:
            self.job_name = self.pdb_file.stem

        if not self.output_dir or str(self.output_dir) == ".":
            self.output_dir = DEFAULT_OUTPUT_DIR / f"md_{self.job_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        if self.forcefield not in FORCEFIELD_MAP:
            raise ValueError(f"Unknown force field: {self.forcefield}")
        self.ff_source = FORCEFIELD_MAP[self.forcefield]

        if self.water_model not in WATER_MODEL_MAP:
            raise ValueError(f"Unknown water model: {self.water_model}")
        self.water_source, self.water_box = WATER_MODEL_MAP[self.water_model]

        # Parse GPU device string (e.g., "cuda:0" -> "0", "cuda:1" -> "1")
        if self.gpu_device.lower().startswith("cuda:"):
            self.cuda_device_id = self.gpu_device.split(":")[-1]
        else:
            # Allow direct device ID specification (e.g., "0", "1")
            self.cuda_device_id = self.gpu_device


class InputFileGenerator:
    """Generate Amber input files for MD simulations."""

    def __init__(self, config: SimulationConfig) -> None:
        self.config = config

    def generate_tleap_input(self) -> str:
        """Generate tleap.in for system preparation."""
        return f"""# Load force field
source {self.config.ff_source}
source {self.config.water_source}

# Load protein
mol = loadpdb {self.config.pdb_file}

# Check for problems
check mol

# Solvate with water box
solvatebox mol {self.config.water_box} {self.config.box_buffer}

# Add ions to neutralize the system
addIons2 mol Na+ 0
addIons2 mol Cl- 0

# Save topology and coordinates
saveamberparm mol system.prmtop system.inpcrd

# Save PDB for visualization
savepdb mol system.pdb

quit
"""

    def generate_minimization_input(self, with_restraints: bool = True) -> str:
        """Generate minimization input file."""
        if with_restraints:
            return """Minimization
 &cntrl
   imin=1,
   maxcyc=5000,
   ncyc=2500,
   ntb=1,
   ntr=1,
   restraint_wt=10.0,
   restraintmask='!@H=',
   cut=10.0,
   ntpr=100,
 /
"""
        else:
            return """Minimization (no restraints)
 &cntrl
   imin=1,
   maxcyc=5000,
   ncyc=2500,
   ntb=1,
   ntr=0,
   cut=10.0,
   ntpr=100,
 /
"""

    def generate_heating_input(self) -> str:
        """Generate heating input file (NVT)."""
        nsteps_heat = 25000
        return f"""Heating from 0 to {self.config.temperature} K
 &cntrl
   imin=0,
   irest=0,
   ntx=1,
   ntb=1,
   cut=10.0,
   ntr=1,
   restraint_wt=5.0,
   restraintmask='@CA',
   nstlim={nsteps_heat},
   dt=0.002,
   ntc=2,
   ntf=2,
   tempi=0.0,
   temp0={self.config.temperature},
   ntt=3,
   gamma_ln=2.0,
   ig=-1,
   ntpr=500,
   ntwx=500,
   ntwr=5000,
   iwrap=1,
   nmropt=1,
 /
 &wt type='TEMP0', istep1=0, istep2={nsteps_heat}, value1=0.0, value2={self.config.temperature}, /
 &wt type='END' /
"""

    def generate_equilibration_input(self) -> str:
        """Generate equilibration input file (NPT)."""
        nsteps_equil = 250000
        return f"""Equilibration (NPT)
 &cntrl
   imin=0,
   irest=1,
   ntx=5,
   ntb=2,
   pres0=1.0,
   ntp=1,
   taup=2.0,
   cut=10.0,
   ntr=1,
   restraint_wt=2.0,
   restraintmask='@CA',
   nstlim={nsteps_equil},
   dt=0.002,
   ntc=2,
   ntf=2,
   temp0={self.config.temperature},
   ntt=3,
   gamma_ln=2.0,
   ig=-1,
   ntpr=500,
   ntwx=500,
   ntwr=10000,
   iwrap=1,
 /
"""

    def generate_production_input(self) -> str:
        """Generate production MD input file (NPT)."""
        nsteps_prod = int(self.config.sim_time_ns * 500000)
        return f"""Production MD (NPT)
 &cntrl
   imin=0,
   irest=1,
   ntx=5,
   ntb=2,
   pres0=1.0,
   ntp=1,
   taup=2.0,
   cut=10.0,
   ntr=0,
   nstlim={nsteps_prod},
   dt=0.002,
   ntc=2,
   ntf=2,
   temp0={self.config.temperature},
   ntt=3,
   gamma_ln=2.0,
   ig=-1,
   ntpr=5000,
   ntwx=5000,
   ntwr=50000,
   iwrap=1,
   ioutfm=1,
 /
"""

    def write_all_input_files(self, output_dir: Path) -> list[dict]:
        """Write all input files to the output directory."""
        files = {
            "tleap.in": self.generate_tleap_input(),
            "min.in": self.generate_minimization_input(with_restraints=True),
            "min2.in": self.generate_minimization_input(with_restraints=False),
            "heat.in": self.generate_heating_input(),
            "equil.in": self.generate_equilibration_input(),
            "prod.in": self.generate_production_input(),
        }

        artifacts = []
        for filename, content in files.items():
            filepath = output_dir / filename
            filepath.write_text(content)
            artifacts.append({
                "description": f"Input file: {filename}",
                "path": str(filepath.resolve())
            })
            logger.debug(f"Created {filepath}")

        return artifacts


class SimulationRunner:
    """Run Amber MD simulations."""

    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.md_engine: Optional[str] = None
        self.env: dict = os.environ.copy()

    def setup_environment(self) -> None:
        """Set up Amber environment variables."""
        amber_sh = AMBER_ENV_DIR / "amber.sh"
        if not amber_sh.exists():
            raise RuntimeError(
                f"Amber environment not found at {amber_sh}. "
                "Run quick_setup.sh first."
            )

        cmd = f"source {amber_sh} && env"
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            executable="/bin/bash",
        )

        if result.returncode != 0:
            raise RuntimeError(f"Failed to source amber.sh: {result.stderr}")

        for line in result.stdout.splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                self.env[key] = value

        # Set CUDA_VISIBLE_DEVICES for GPU selection
        if self.config.use_gpu:
            self.env["CUDA_VISIBLE_DEVICES"] = self.config.cuda_device_id
            logger.info(f"GPU device set: CUDA_VISIBLE_DEVICES={self.config.cuda_device_id}")

        logger.info("Amber environment configured")

    def detect_md_engine(self) -> str:
        """Detect the best available MD engine."""
        if self.config.use_gpu:
            result = subprocess.run(
                ["which", "pmemd.cuda"],
                capture_output=True,
                env=self.env,
            )
            if result.returncode == 0:
                self.md_engine = "pmemd.cuda"
                logger.info("Using GPU-accelerated pmemd.cuda")
                return self.md_engine

        result = subprocess.run(
            ["which", "pmemd"],
            capture_output=True,
            env=self.env,
        )
        if result.returncode == 0:
            self.md_engine = "pmemd"
            logger.info("Using CPU pmemd")
            return self.md_engine

        result = subprocess.run(
            ["which", "sander"],
            capture_output=True,
            env=self.env,
        )
        if result.returncode == 0:
            self.md_engine = "sander"
            logger.warning("Using sander (slower than pmemd)")
            return self.md_engine

        raise RuntimeError("No Amber MD engine found")

    def run_command(
        self,
        cmd: list[str],
        step_name: str,
        check_output: Optional[Path] = None,
    ) -> None:
        """Run a command and check for success."""
        logger.info(f"Running {step_name}...")

        result = subprocess.run(
            cmd,
            env=self.env,
            cwd=self.config.output_dir,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            logger.error(f"{step_name} failed!")
            logger.error(f"Stderr: {result.stderr}")
            raise RuntimeError(f"{step_name} failed")

        if check_output and not check_output.exists():
            raise RuntimeError(f"{step_name} output not found: {check_output}")

        logger.success(f"{step_name} completed")

    def run_tleap(self) -> int:
        """Run tleap to prepare the system. Returns atom count."""
        cmd = ["tleap", "-f", "tleap.in"]
        self.run_command(cmd, "tleap")

        prmtop = self.config.output_dir / "system.prmtop"
        inpcrd = self.config.output_dir / "system.inpcrd"

        if not prmtop.exists() or not inpcrd.exists():
            raise RuntimeError("tleap failed to create output files")

        system_pdb = self.config.output_dir / "system.pdb"
        natoms = 0
        if system_pdb.exists():
            natoms = sum(
                1 for line in system_pdb.read_text().splitlines()
                if line.startswith("ATOM")
            )
        return natoms

    def run_minimization(self, stage: int = 1) -> None:
        """Run energy minimization."""
        if stage == 1:
            input_file, output_file = "min.in", "min"
            coords, restart = "system.inpcrd", "min.rst7"
            extra_args = ["-ref", "system.inpcrd"]
            step_name = "Minimization (with restraints)"
        else:
            input_file, output_file = "min2.in", "min2"
            coords, restart = "min.rst7", "min2.rst7"
            extra_args = []
            step_name = "Minimization (no restraints)"

        cmd = [
            self.md_engine, "-O",
            "-i", input_file,
            "-o", f"{output_file}.out",
            "-p", "system.prmtop",
            "-c", coords,
            "-r", restart,
        ] + extra_args

        self.run_command(cmd, step_name, self.config.output_dir / restart)

    def run_heating(self) -> None:
        """Run heating simulation (NVT)."""
        cmd = [
            self.md_engine, "-O",
            "-i", "heat.in",
            "-o", "heat.out",
            "-p", "system.prmtop",
            "-c", "min2.rst7",
            "-r", "heat.rst7",
            "-x", "heat.nc",
            "-ref", "min2.rst7",
        ]
        self.run_command(cmd, f"Heating (0 -> {self.config.temperature} K)",
                         self.config.output_dir / "heat.rst7")

    def run_equilibration(self) -> None:
        """Run equilibration simulation (NPT)."""
        cmd = [
            self.md_engine, "-O",
            "-i", "equil.in",
            "-o", "equil.out",
            "-p", "system.prmtop",
            "-c", "heat.rst7",
            "-r", "equil.rst7",
            "-x", "equil.nc",
            "-ref", "heat.rst7",
        ]
        self.run_command(cmd, "Equilibration (NPT, 500 ps)",
                         self.config.output_dir / "equil.rst7")

    def run_production(self) -> None:
        """Run production MD simulation (NPT)."""
        cmd = [
            self.md_engine, "-O",
            "-i", "prod.in",
            "-o", "prod.out",
            "-p", "system.prmtop",
            "-c", "equil.rst7",
            "-r", "prod.rst7",
            "-x", "prod.nc",
        ]
        self.run_command(cmd, f"Production MD (NPT, {self.config.sim_time_ns} ns)",
                         self.config.output_dir / "prod.rst7")


@simulation_mcp.tool
def amber_run_protein_md(
    pdb_file: Annotated[str, "Path to input PDB file containing the protein structure"],
    job_name: Annotated[str | None, "Job name for output files (default: derived from PDB filename)"] = None,
    sim_time_ns: Annotated[float, "Production simulation time in nanoseconds"] = 10.0,
    temperature: Annotated[float, "Simulation temperature in Kelvin"] = 300.0,
    box_buffer: Annotated[float, "Water box buffer size in Angstroms"] = 12.0,
    salt_conc: Annotated[float, "Salt concentration in Molar (for ion addition)"] = 0.15,
    forcefield: Annotated[str, "Force field to use: 'ff14SB' or 'ff19SB'"] = "ff19SB",
    water_model: Annotated[str, "Water model: 'tip3p', 'opc', or 'tip4pew'"] = "opc",
    use_gpu: Annotated[bool, "Use GPU acceleration with pmemd.cuda if available"] = True,
    gpu_device: Annotated[str, "GPU device to use (e.g., 'cuda:0', 'cuda:1', or just '0', '1')"] = "cuda:0",
    output_dir: Annotated[str | None, "Output directory path (default: results/md_<jobname>_<timestamp>)"] = None,
    dry_run: Annotated[bool, "Generate input files only without running simulations"] = False,
) -> dict:
    """
    Run a complete MD simulation workflow for a single protein.

    This tool performs the full Amber MD simulation pipeline:
    1. System preparation with tleap (solvation, ion addition)
    2. Energy minimization (with and without restraints)
    3. Heating from 0K to target temperature (NVT)
    4. Equilibration at constant pressure (NPT, 500 ps)
    5. Production MD at constant pressure (NPT)

    Input is a PDB file, output is trajectory files and restart files for analysis.
    """
    logger.info("=" * 60)
    logger.info("Amber MD Simulation - Single Protein Workflow")
    logger.info("=" * 60)

    # Validate PDB file
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Create configuration
    config = SimulationConfig(
        pdb_file=pdb_path,
        job_name=job_name or "",
        sim_time_ns=sim_time_ns,
        temperature=temperature,
        box_buffer=box_buffer,
        salt_conc=salt_conc,
        forcefield=forcefield,
        water_model=water_model,
        use_gpu=use_gpu,
        gpu_device=gpu_device,
        dry_run=dry_run,
        output_dir=Path(output_dir) if output_dir else Path("."),
    )

    # Log configuration
    logger.info(f"Input PDB: {config.pdb_file}")
    logger.info(f"Job name: {config.job_name}")
    logger.info(f"Output dir: {config.output_dir}")
    logger.info(f"Force field: {config.forcefield}")
    logger.info(f"Water model: {config.water_model}")
    logger.info(f"Temperature: {config.temperature} K")
    logger.info(f"Box buffer: {config.box_buffer} A")
    logger.info(f"Simulation time: {config.sim_time_ns} ns")
    logger.info(f"GPU device: {config.gpu_device} (CUDA device: {config.cuda_device_id})")
    logger.info(f"Dry run: {config.dry_run}")

    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate input files
    generator = InputFileGenerator(config)
    artifacts = generator.write_all_input_files(config.output_dir)
    logger.success("Input files generated")

    if config.dry_run:
        logger.info("[DRY-RUN] Simulation files generated, skipping execution")
        return {
            "message": f"Dry run completed. Input files generated in {config.output_dir}",
            "artifacts": artifacts,
            "config": {
                "job_name": config.job_name,
                "output_dir": str(config.output_dir),
                "forcefield": config.forcefield,
                "water_model": config.water_model,
                "temperature": config.temperature,
                "sim_time_ns": config.sim_time_ns,
                "gpu_device": config.gpu_device,
            }
        }

    # Run simulation
    runner = SimulationRunner(config)
    runner.setup_environment()
    runner.detect_md_engine()

    # Execute workflow
    natoms = runner.run_tleap()
    logger.success(f"System prepared: {natoms} atoms")
    artifacts.append({
        "description": "Topology file",
        "path": str((config.output_dir / "system.prmtop").resolve())
    })
    artifacts.append({
        "description": "Initial coordinates",
        "path": str((config.output_dir / "system.inpcrd").resolve())
    })

    runner.run_minimization(stage=1)
    runner.run_minimization(stage=2)
    runner.run_heating()
    runner.run_equilibration()
    runner.run_production()

    # Add production outputs to artifacts
    artifacts.append({
        "description": "Production trajectory (NetCDF)",
        "path": str((config.output_dir / "prod.nc").resolve())
    })
    artifacts.append({
        "description": "Final restart file",
        "path": str((config.output_dir / "prod.rst7").resolve())
    })
    artifacts.append({
        "description": "Production output log",
        "path": str((config.output_dir / "prod.out").resolve())
    })

    logger.success("MD simulation completed successfully!")

    return {
        "message": f"MD simulation completed: {config.sim_time_ns} ns at {config.temperature} K",
        "artifacts": artifacts,
        "config": {
            "job_name": config.job_name,
            "output_dir": str(config.output_dir),
            "forcefield": config.forcefield,
            "water_model": config.water_model,
            "temperature": config.temperature,
            "sim_time_ns": config.sim_time_ns,
            "gpu_device": config.gpu_device,
            "natoms": natoms,
            "md_engine": runner.md_engine,
        }
    }


@simulation_mcp.tool
def amber_prepare_system(
    pdb_file: Annotated[str, "Path to input PDB file containing the protein structure"],
    job_name: Annotated[str | None, "Job name for output files"] = None,
    box_buffer: Annotated[float, "Water box buffer size in Angstroms"] = 12.0,
    forcefield: Annotated[str, "Force field: 'ff14SB' or 'ff19SB'"] = "ff19SB",
    water_model: Annotated[str, "Water model: 'tip3p', 'opc', or 'tip4pew'"] = "opc",
    output_dir: Annotated[str | None, "Output directory path"] = None,
) -> dict:
    """
    Prepare a protein system with tleap for MD simulation.

    This tool performs system preparation only:
    1. Load protein structure from PDB
    2. Apply force field parameters
    3. Solvate in a water box
    4. Add neutralizing ions

    Output is topology (.prmtop) and coordinate (.inpcrd) files ready for simulation.
    """
    logger.info("=" * 60)
    logger.info("Amber System Preparation with tleap")
    logger.info("=" * 60)

    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    config = SimulationConfig(
        pdb_file=pdb_path,
        job_name=job_name or "",
        box_buffer=box_buffer,
        forcefield=forcefield,
        water_model=water_model,
        output_dir=Path(output_dir) if output_dir else Path("."),
    )

    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate tleap input only
    generator = InputFileGenerator(config)
    tleap_content = generator.generate_tleap_input()
    tleap_file = config.output_dir / "tleap.in"
    tleap_file.write_text(tleap_content)

    # Run tleap
    runner = SimulationRunner(config)
    runner.setup_environment()
    natoms = runner.run_tleap()

    artifacts = [
        {"description": "tleap input", "path": str(tleap_file.resolve())},
        {"description": "Topology file", "path": str((config.output_dir / "system.prmtop").resolve())},
        {"description": "Coordinate file", "path": str((config.output_dir / "system.inpcrd").resolve())},
        {"description": "Solvated PDB", "path": str((config.output_dir / "system.pdb").resolve())},
    ]

    logger.success(f"System prepared: {natoms} atoms")

    return {
        "message": f"System prepared successfully: {natoms} atoms in {config.water_box} box",
        "artifacts": artifacts,
        "config": {
            "job_name": config.job_name,
            "output_dir": str(config.output_dir),
            "forcefield": config.forcefield,
            "water_model": config.water_model,
            "box_buffer": config.box_buffer,
            "natoms": natoms,
        }
    }


@simulation_mcp.tool
def amber_generate_input_files(
    pdb_file: Annotated[str, "Path to input PDB file (used to determine absolute path in tleap.in)"],
    job_name: Annotated[str | None, "Job name for output files"] = None,
    sim_time_ns: Annotated[float, "Production simulation time in nanoseconds"] = 10.0,
    temperature: Annotated[float, "Simulation temperature in Kelvin"] = 300.0,
    box_buffer: Annotated[float, "Water box buffer size in Angstroms"] = 12.0,
    forcefield: Annotated[str, "Force field: 'ff14SB' or 'ff19SB'"] = "ff19SB",
    water_model: Annotated[str, "Water model: 'tip3p', 'opc', or 'tip4pew'"] = "opc",
    output_dir: Annotated[str | None, "Output directory path"] = None,
) -> dict:
    """
    Generate Amber input files without running simulations.

    Creates all input files needed for a complete MD workflow:
    - tleap.in: System preparation script
    - min.in, min2.in: Minimization inputs
    - heat.in: Heating input (NVT)
    - equil.in: Equilibration input (NPT)
    - prod.in: Production MD input (NPT)

    Use this to inspect or customize input files before running simulations.
    """
    logger.info("=" * 60)
    logger.info("Generating Amber Input Files")
    logger.info("=" * 60)

    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    config = SimulationConfig(
        pdb_file=pdb_path,
        job_name=job_name or "",
        sim_time_ns=sim_time_ns,
        temperature=temperature,
        box_buffer=box_buffer,
        forcefield=forcefield,
        water_model=water_model,
        output_dir=Path(output_dir) if output_dir else Path("."),
    )

    config.output_dir.mkdir(parents=True, exist_ok=True)

    generator = InputFileGenerator(config)
    artifacts = generator.write_all_input_files(config.output_dir)

    logger.success(f"Generated {len(artifacts)} input files in {config.output_dir}")

    return {
        "message": f"Input files generated in {config.output_dir}",
        "artifacts": artifacts,
        "config": {
            "job_name": config.job_name,
            "output_dir": str(config.output_dir),
            "forcefield": config.forcefield,
            "water_model": config.water_model,
            "temperature": config.temperature,
            "sim_time_ns": config.sim_time_ns,
            "nsteps_production": int(config.sim_time_ns * 500000),
        }
    }
