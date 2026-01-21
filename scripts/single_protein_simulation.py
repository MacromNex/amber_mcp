#!/usr/bin/env python3
"""
single_protein_simulation.py - Run a single protein MD simulation with Amber

This script automates the process of running a molecular dynamics simulation
for a single protein molecule, including:
  1. System preparation (tleap)
  2. Energy minimization
  3. Heating
  4. Equilibration (NPT)
  5. Production MD

Usage:
    python single_protein_simulation.py <protein.pdb> [options]

Examples:
    python single_protein_simulation.py protein.pdb
    python single_protein_simulation.py protein.pdb -n my_sim -t 100 -T 310
    python single_protein_simulation.py protein.pdb --forcefield ff14SB --water tip3p
"""

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from loguru import logger

# Configure loguru
logger.remove()
logger.add(
    sys.stderr,
    format="<level>[{level}]</level> {message}",
    level="INFO",
    colorize=True,
)


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
    dry_run: bool = False
    output_dir: Path = field(default_factory=Path)

    # Derived fields (set in __post_init__)
    ff_source: str = field(init=False, default="")
    water_source: str = field(init=False, default="")
    water_box: str = field(init=False, default="")

    def __post_init__(self) -> None:
        """Validate inputs and set derived fields."""
        # Convert to Path if string
        if isinstance(self.pdb_file, str):
            self.pdb_file = Path(self.pdb_file)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

        # Validate PDB file exists
        if not self.pdb_file.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_file}")

        # Get absolute path
        self.pdb_file = self.pdb_file.resolve()

        # Set job name from PDB filename if not specified
        if not self.job_name:
            self.job_name = self.pdb_file.stem

        # Set output directory
        if not self.output_dir or str(self.output_dir) == ".":
            self.output_dir = Path(f"./md_{self.job_name}")

        # Validate and map force field
        if self.forcefield not in FORCEFIELD_MAP:
            raise ValueError(
                f"Unknown force field: {self.forcefield} "
                f"(supported: {', '.join(set(FORCEFIELD_MAP.keys()))})"
            )
        self.ff_source = FORCEFIELD_MAP[self.forcefield]

        # Validate and map water model
        if self.water_model not in WATER_MODEL_MAP:
            raise ValueError(
                f"Unknown water model: {self.water_model} "
                f"(supported: {', '.join(set(WATER_MODEL_MAP.keys()))})"
            )
        self.water_source, self.water_box = WATER_MODEL_MAP[self.water_model]

        # Validate numeric parameters
        if self.sim_time_ns <= 0:
            raise ValueError("Simulation time must be positive")
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")
        if self.box_buffer <= 0:
            raise ValueError("Box buffer must be positive")
        if self.salt_conc < 0:
            raise ValueError("Salt concentration cannot be negative")


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
# addIons2 handles both positive and negative systems automatically
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
   imin=1,           ! Minimization
   maxcyc=5000,      ! Max cycles
   ncyc=2500,        ! Steepest descent cycles, then conjugate gradient
   ntb=1,            ! Constant volume PBC
   ntr=1,            ! Restrain heavy atoms
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
        nsteps_heat = 25000  # 50 ps heating
        return f"""Heating from 0 to {self.config.temperature} K
 &cntrl
   imin=0,           ! MD
   irest=0,          ! New simulation
   ntx=1,            ! Read coordinates only
   ntb=1,            ! Constant volume PBC
   cut=10.0,
   ntr=1,            ! Restrain protein
   restraint_wt=5.0,
   restraintmask='@CA',
   nstlim={nsteps_heat},
   dt=0.002,         ! 2 fs timestep
   ntc=2,            ! SHAKE on hydrogens
   ntf=2,            ! No force calc on H bonds
   tempi=0.0,
   temp0={self.config.temperature},
   ntt=3,            ! Langevin thermostat
   gamma_ln=2.0,
   ig=-1,            ! Random seed
   ntpr=500,
   ntwx=500,
   ntwr=5000,
   iwrap=1,
   nmropt=1,         ! NMR restraints for temperature ramp
 /
 &wt type='TEMP0', istep1=0, istep2={nsteps_heat}, value1=0.0, value2={self.config.temperature}, /
 &wt type='END' /
"""

    def generate_equilibration_input(self) -> str:
        """Generate equilibration input file (NPT)."""
        nsteps_equil = 250000  # 500 ps equilibration
        return f"""Equilibration (NPT)
 &cntrl
   imin=0,
   irest=1,          ! Restart
   ntx=5,            ! Read coordinates and velocities
   ntb=2,            ! Constant pressure PBC
   pres0=1.0,        ! 1 atm
   ntp=1,            ! Isotropic pressure scaling
   taup=2.0,         ! Pressure relaxation time
   cut=10.0,
   ntr=1,            ! Restrain CA atoms
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
        # Calculate steps (2 fs timestep): ns * 1e6 fs / 2 fs = ns * 500000
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
   ntr=0,            ! No restraints
   nstlim={nsteps_prod},
   dt=0.002,
   ntc=2,
   ntf=2,
   temp0={self.config.temperature},
   ntt=3,
   gamma_ln=2.0,
   ig=-1,
   ntpr=5000,        ! Energy output every 10 ps
   ntwx=5000,        ! Trajectory every 10 ps
   ntwr=50000,       ! Restart every 100 ps
   iwrap=1,
   ioutfm=1,         ! NetCDF trajectory format
 /
"""

    def write_all_input_files(self, output_dir: Path) -> None:
        """Write all input files to the output directory."""
        files = {
            "tleap.in": self.generate_tleap_input(),
            "min.in": self.generate_minimization_input(with_restraints=True),
            "min2.in": self.generate_minimization_input(with_restraints=False),
            "heat.in": self.generate_heating_input(),
            "equil.in": self.generate_equilibration_input(),
            "prod.in": self.generate_production_input(),
        }

        for filename, content in files.items():
            filepath = output_dir / filename
            filepath.write_text(content)
            logger.debug(f"Created {filepath}")

        logger.success("Input files created")


class SimulationRunner:
    """Run Amber MD simulations."""

    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.md_engine: Optional[str] = None
        self.amber_env: Optional[Path] = None
        self.env: dict = os.environ.copy()

    def setup_environment(self) -> None:
        """Set up Amber environment variables."""
        # Find amber.sh relative to script location
        script_dir = Path(__file__).parent.parent
        amber_env_dir = script_dir / "env"

        if not (amber_env_dir / "amber.sh").exists():
            raise RuntimeError(
                f"Amber environment not found at {amber_env_dir}/amber.sh. "
                "Run quick_setup.sh first."
            )

        self.amber_env = amber_env_dir

        # Source amber.sh and capture environment
        cmd = f"source {amber_env_dir}/amber.sh && env"
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            executable="/bin/bash",
        )

        if result.returncode != 0:
            raise RuntimeError(f"Failed to source amber.sh: {result.stderr}")

        # Parse environment variables
        for line in result.stdout.splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                self.env[key] = value

        logger.info("Amber environment configured")

    def detect_md_engine(self) -> str:
        """Detect the best available MD engine."""
        if self.config.use_gpu:
            # Check for pmemd.cuda
            result = subprocess.run(
                ["which", "pmemd.cuda"],
                capture_output=True,
                env=self.env,
            )
            if result.returncode == 0:
                self.md_engine = "pmemd.cuda"
                logger.info("Using GPU-accelerated pmemd.cuda")
                return self.md_engine

        # Check for pmemd
        result = subprocess.run(
            ["which", "pmemd"],
            capture_output=True,
            env=self.env,
        )
        if result.returncode == 0:
            self.md_engine = "pmemd"
            logger.info("Using CPU pmemd")
            return self.md_engine

        # Fall back to sander
        result = subprocess.run(
            ["which", "sander"],
            capture_output=True,
            env=self.env,
        )
        if result.returncode == 0:
            self.md_engine = "sander"
            logger.warning("Using sander (slower than pmemd)")
            return self.md_engine

        raise RuntimeError("No Amber MD engine found (pmemd.cuda, pmemd, or sander)")

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
            logger.error(f"Command: {' '.join(cmd)}")
            logger.error(f"Stderr: {result.stderr}")
            raise RuntimeError(f"{step_name} failed. Check output files for details.")

        if check_output and not check_output.exists():
            raise RuntimeError(
                f"{step_name} completed but output file not found: {check_output}"
            )

        logger.success(f"{step_name} completed")

    def run_tleap(self) -> None:
        """Run tleap to prepare the system."""
        cmd = ["tleap", "-f", "tleap.in"]
        self.run_command(cmd, "tleap")

        # Verify outputs
        prmtop = self.config.output_dir / "system.prmtop"
        inpcrd = self.config.output_dir / "system.inpcrd"

        if not prmtop.exists() or not inpcrd.exists():
            log_file = self.config.output_dir / "tleap.log"
            if log_file.exists():
                logger.error(f"tleap log:\n{log_file.read_text()}")
            raise RuntimeError("tleap failed. Check tleap.log for details.")

        # Count atoms in system
        system_pdb = self.config.output_dir / "system.pdb"
        if system_pdb.exists():
            natoms = sum(
                1 for line in system_pdb.read_text().splitlines()
                if line.startswith("ATOM")
            )
            logger.success(f"System prepared: {natoms} atoms")

    def run_minimization(self, stage: int = 1) -> None:
        """Run energy minimization."""
        if stage == 1:
            input_file = "min.in"
            output_file = "min"
            coords = "system.inpcrd"
            restart = "min.rst7"
            extra_args = ["-ref", "system.inpcrd"]
            step_name = "Minimization (with restraints)"
        else:
            input_file = "min2.in"
            output_file = "min2"
            coords = "min.rst7"
            restart = "min2.rst7"
            extra_args = []
            step_name = "Minimization (no restraints)"

        cmd = [
            self.md_engine,
            "-O",
            "-i", input_file,
            "-o", f"{output_file}.out",
            "-p", "system.prmtop",
            "-c", coords,
            "-r", restart,
        ] + extra_args

        self.run_command(
            cmd,
            step_name,
            check_output=self.config.output_dir / restart,
        )

    def run_heating(self) -> None:
        """Run heating simulation (NVT)."""
        cmd = [
            self.md_engine,
            "-O",
            "-i", "heat.in",
            "-o", "heat.out",
            "-p", "system.prmtop",
            "-c", "min2.rst7",
            "-r", "heat.rst7",
            "-x", "heat.nc",
            "-ref", "min2.rst7",
        ]

        self.run_command(
            cmd,
            f"Heating (0 -> {self.config.temperature} K)",
            check_output=self.config.output_dir / "heat.rst7",
        )

    def run_equilibration(self) -> None:
        """Run equilibration simulation (NPT)."""
        cmd = [
            self.md_engine,
            "-O",
            "-i", "equil.in",
            "-o", "equil.out",
            "-p", "system.prmtop",
            "-c", "heat.rst7",
            "-r", "equil.rst7",
            "-x", "equil.nc",
            "-ref", "heat.rst7",
        ]

        self.run_command(
            cmd,
            "Equilibration (NPT, 500 ps)",
            check_output=self.config.output_dir / "equil.rst7",
        )

    def run_production(self) -> None:
        """Run production MD simulation (NPT)."""
        cmd = [
            self.md_engine,
            "-O",
            "-i", "prod.in",
            "-o", "prod.out",
            "-p", "system.prmtop",
            "-c", "equil.rst7",
            "-r", "prod.rst7",
            "-x", "prod.nc",
        ]

        self.run_command(
            cmd,
            f"Production MD (NPT, {self.config.sim_time_ns} ns)",
            check_output=self.config.output_dir / "prod.rst7",
        )

    def run_all(self) -> None:
        """Execute the full MD workflow."""
        # Setup
        self.setup_environment()
        self.detect_md_engine()

        # Create output directory
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Generate input files
        generator = InputFileGenerator(self.config)
        generator.write_all_input_files(self.config.output_dir)

        if self.config.dry_run:
            logger.info("[DRY-RUN] Would run the following simulations:")
            print("  1. Minimization with restraints")
            print("  2. Minimization without restraints")
            print(f"  3. Heating (0 -> {self.config.temperature} K)")
            print("  4. Equilibration (NPT, 500 ps)")
            print(f"  5. Production (NPT, {self.config.sim_time_ns} ns)")
            print()
            logger.info(f"Files created in: {self.config.output_dir}")
            return

        # Run simulations
        self.run_tleap()
        self.run_minimization(stage=1)
        self.run_minimization(stage=2)
        self.run_heating()
        self.run_equilibration()
        self.run_production()

        # Print summary
        self.print_summary()

    def print_summary(self) -> None:
        """Print simulation summary and analysis hints."""
        print()
        print("=" * 46)
        print("  MD Simulation Complete!")
        print("=" * 46)
        print()
        print(f"Output files in: {self.config.output_dir}")
        print()
        print("Key files:")
        print("  - system.prmtop    : Topology file")
        print("  - system.inpcrd    : Initial coordinates")
        print("  - prod.rst7        : Final restart file")
        print("  - prod.nc          : Production trajectory")
        print("  - prod.out         : Production output")
        print()
        print("Analysis commands:")
        print("  # Load trajectory in cpptraj")
        print("  cpptraj -p system.prmtop -y prod.nc")
        print()
        print("  # Calculate RMSD")
        print("  cpptraj << EOF")
        print("  parm system.prmtop")
        print("  trajin prod.nc")
        print("  rms first @CA out rmsd.dat")
        print("  run")
        print("  EOF")
        print()
        print("  # Extract frames as PDB")
        print("  cpptraj -p system.prmtop -y prod.nc -x frames.pdb")
        print()


def print_config(config: SimulationConfig) -> None:
    """Print simulation configuration."""
    print()
    print("=" * 46)
    print("  Amber MD Simulation Setup")
    print("=" * 46)
    print()
    print(f"Input PDB:       {config.pdb_file}")
    print(f"Job name:        {config.job_name}")
    print(f"Output dir:      {config.output_dir}")
    print(f"Force field:     {config.forcefield}")
    print(f"Water model:     {config.water_model}")
    print(f"Temperature:     {config.temperature} K")
    print(f"Box buffer:      {config.box_buffer} \u00c5")
    print(f"Salt conc:       {config.salt_conc} M")
    print(f"Simulation:      {config.sim_time_ns} ns")
    print()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run a single protein MD simulation with Amber",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s protein.pdb
  %(prog)s protein.pdb -n my_sim -t 100 -T 310
  %(prog)s protein.pdb --forcefield ff14SB --water tip3p
  %(prog)s protein.pdb --dry-run
""",
    )

    parser.add_argument(
        "pdb_file",
        type=Path,
        help="Input PDB file",
    )
    parser.add_argument(
        "-n", "--name",
        dest="job_name",
        default="",
        help="Job name (default: from PDB filename)",
    )
    parser.add_argument(
        "-t", "--time",
        dest="sim_time_ns",
        type=float,
        default=10.0,
        help="Production simulation time in ns (default: 10)",
    )
    parser.add_argument(
        "-T", "--temp",
        dest="temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin (default: 300)",
    )
    parser.add_argument(
        "-b", "--box",
        dest="box_buffer",
        type=float,
        default=12.0,
        help="Box buffer size in Angstrom (default: 12)",
    )
    parser.add_argument(
        "-s", "--salt",
        dest="salt_conc",
        type=float,
        default=0.15,
        help="Salt concentration in M (default: 0.15)",
    )
    parser.add_argument(
        "-f", "--forcefield",
        default="ff19SB",
        choices=["ff14SB", "ff19SB"],
        help="Force field (default: ff19SB)",
    )
    parser.add_argument(
        "-w", "--water",
        dest="water_model",
        default="opc",
        choices=["tip3p", "opc", "tip4pew"],
        help="Water model (default: opc)",
    )
    parser.add_argument(
        "-c", "--cpu",
        dest="use_cpu",
        action="store_true",
        help="Force CPU execution (no GPU)",
    )
    parser.add_argument(
        "-d", "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Generate files but don't run simulations",
    )
    parser.add_argument(
        "-o", "--outdir",
        dest="output_dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: ./md_<name>)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose/debug output",
    )

    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = parse_args()

    # Configure logging level
    if args.verbose:
        logger.remove()
        logger.add(
            sys.stderr,
            format="<level>[{level}]</level> {message}",
            level="DEBUG",
            colorize=True,
        )

    try:
        # Create configuration
        config = SimulationConfig(
            pdb_file=args.pdb_file,
            job_name=args.job_name,
            sim_time_ns=args.sim_time_ns,
            temperature=args.temperature,
            box_buffer=args.box_buffer,
            salt_conc=args.salt_conc,
            forcefield=args.forcefield,
            water_model=args.water_model,
            use_gpu=not args.use_cpu,
            dry_run=args.dry_run,
            output_dir=args.output_dir if str(args.output_dir) != "." else Path("."),
        )

        # Print configuration
        print_config(config)

        # Run simulation
        runner = SimulationRunner(config)
        runner.run_all()

        return 0

    except FileNotFoundError as e:
        logger.error(str(e))
        return 1
    except ValueError as e:
        logger.error(str(e))
        return 1
    except RuntimeError as e:
        logger.error(str(e))
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user")
        return 130


if __name__ == "__main__":
    sys.exit(main())
