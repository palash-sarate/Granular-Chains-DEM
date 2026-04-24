import os
import sys
import argparse
from typing import List, Optional

from simulation.orchestrator import SimulationOrchestrator

def main():
    parser = argparse.ArgumentParser(description="LAMMPS Chain Simulation CLI")
    parser.add_argument("--lammps-exe", default="lmp", help="Path to LAMMPS executable (default: lmp)")
    
    subparsers = parser.add_subparsers(dest="command", help="Simulation commands", required=True)

    # 1. run_hopper_fill
    p_fill = subparsers.add_parser("run_hopper_fill", help="Pre-fill a hopper with relaxed molecular chains")
    p_fill.add_argument("--source_dir", default="chain_data/relaxed/N4", help="Directory containing relaxed chains")
    p_fill.add_argument("--fill_template", default="in.hopper_fill", help="LAMMPS input template name")
    p_fill.add_argument("--n_fill", type=int, default=10, help="Number of chains to insert")
    p_fill.add_argument("--relax_steps", type=int, default=100000, help="Number of relaxation steps")
    p_fill.add_argument("--run_name", help="Custom run name")
    p_fill.add_argument("--seed", type=int, help="Random seed")
    p_fill.add_argument("--dt", type=float, default=1e-6, help="Timestep size")
    p_fill.add_argument("--mol_dir", default="chain_data/molecules_temp", help="Temporary molecules directory")
    p_fill.add_argument("--setup_inc", default="simulation_geometries/2D_hopper.inc", help="Geometry setup file")
    p_fill.add_argument("--dump_inc", default="simulation_templates/default_dump.inc", help="Dump settings file")
    p_fill.add_argument("--viscosity", type=float, default=0.001, help="Simulated viscosity")
    p_fill.add_argument("--N", type=int, default=4, help="Chain length (beads)")
    p_fill.add_argument("--outdir", help="Output directory override")
    p_fill.add_argument("--num_procs", type=int, help="Number of MPI processes")
    p_fill.add_argument("--num_threads", type=int, default=1, help="Number of OpenMP threads")
    p_fill.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_fill.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    # 2. resume_hopper_fill
    p_resume_h = subparsers.add_parser("resume_hopper_fill", help="Resume hopper filling from a restart file")
    p_resume_h.add_argument("--restart_path", required=True, help="Path to .bin restart file")
    p_resume_h.add_argument("--source_dir", default="chain_data/relaxed/N4")
    p_resume_h.add_argument("--fill_template", default="in.hopper_fill_resume")
    p_resume_h.add_argument("--n_fill", type=int, default=10)
    p_resume_h.add_argument("--relax_steps", type=int, default=100000)
    p_resume_h.add_argument("--run_name")
    p_resume_h.add_argument("--seed", type=int)
    p_resume_h.add_argument("--dt", type=float, default=1e-6)
    p_resume_h.add_argument("--mol_dir", default="chain_data/molecules_temp")
    p_resume_h.add_argument("--setup_inc", default="simulation_geometries/2D_hopper.inc")
    p_resume_h.add_argument("--dump_inc", default="simulation_templates/default_dump.inc")
    p_resume_h.add_argument("--viscosity", type=float, default=0.001)
    p_resume_h.add_argument("--N", type=int, default=4)
    p_resume_h.add_argument("--outdir")
    p_resume_h.add_argument("--num_procs", type=int)
    p_resume_h.add_argument("--num_threads", type=int, default=1)
    p_resume_h.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_resume_h.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    # 3. run_flop_simulation
    p_flop = subparsers.add_parser("run_flop_simulation", help="Run a single chain-flop mobility simulation")
    p_flop.add_argument("--N", type=int, default=4)
    p_flop.add_argument("--run_steps", type=int, default=50000)
    p_flop.add_argument("--viscosity", type=float, default=0.001)
    p_flop.add_argument("--dt", type=float, default=1e-6)
    p_flop.add_argument("--num_procs", type=int)
    p_flop.add_argument("--num_threads", type=int, default=1)
    p_flop.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_flop.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    # 4. resume_flop_simulation
    p_resume_f = subparsers.add_parser("resume_flop_simulation", help="Resume a chain-flop simulation")
    p_resume_f.add_argument("--N", type=int, default=4)
    p_resume_f.add_argument("--run_steps", type=int, default=50000)
    p_resume_f.add_argument("--viscosity", type=float, default=0.001)
    p_resume_f.add_argument("--dt", type=float, default=1e-6)
    p_resume_f.add_argument("--resume_token", required=True, help="Timestep token of the restart file")
    p_resume_f.add_argument("--num_procs", type=int)
    p_resume_f.add_argument("--num_threads", type=int, default=1)
    p_resume_f.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_resume_f.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    # 5. generate_relaxed_library
    p_lib = subparsers.add_parser("generate_relaxed_library", help="Generate relaxed chain states library")
    p_lib.add_argument("--n_beads", type=int, default=4)
    p_lib.add_argument("--n_states", type=int, default=10)
    p_lib.add_argument("--forced", action="store_true", help="Force re-generation if check exists")
    p_lib.add_argument("--dump_inc", default="simulation_templates/default_dump.inc", help="Dump settings file")
    p_lib.add_argument("--inParallel", type=int, default=1, dest="n_parallel", help="Number of simulations to run in parallel")
    p_lib.add_argument("--num_procs", type=int)
    p_lib.add_argument("--num_threads", type=int, default=1)
    p_lib.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_lib.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    # 6. generate_linear_chains
    p_chains = subparsers.add_parser("generate_linear_chains", help="Generate initial linear chain data files")
    p_chains.add_argument("--Ns", default="4,6,8", help="Comma-separated chain lengths")
    p_chains.add_argument("--orientation", default="horz", choices=["horz", "vert"])
    p_chains.add_argument("--output_dir_name", default="linear_x")

    # 7. run_flop_batch
    p_batch = subparsers.add_parser("run_flop_batch", help="Run a batch of flop simulations")
    p_batch.add_argument("--Ns", default="4,6", help="Comma-separated chain lengths")
    p_batch.add_argument("--run_steps", default="50000", help="Comma-separated run steps")
    p_batch.add_argument("--viscosities", default="0.001", help="Comma-separated viscosities")
    p_batch.add_argument("--dt", type=float, default=1e-6)
    p_batch.add_argument("--num_procs", type=int)
    p_batch.add_argument("--num_threads", type=int, default=1)
    p_batch.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable KOKKOS acceleration")
    p_batch.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable INTEL acceleration")

    args = parser.parse_args()
    
    # Initialize orchestrator
    orchestrator = SimulationOrchestrator(lammps_executable=args.lammps_exe)
    
    # Extract command and arguments
    cmd_name = args.command
    cmd_args = vars(args)
    
    # Remove metadata args not passed to methods
    cmd_args.pop("command")
    cmd_args.pop("lammps_exe")
    
    # Handle list-based arguments for run_flop_batch and generate_linear_chains
    if cmd_name == "run_flop_batch":
        cmd_args["Ns"] = [int(n.strip()) for n in str(cmd_args["Ns"]).split(",")]
        cmd_args["run_steps"] = [int(s.strip()) for s in str(cmd_args["run_steps"]).split(",")]
        cmd_args["viscosities"] = [float(v.strip()) for v in str(cmd_args["viscosities"]).split(",")]
    
    if cmd_name == "generate_linear_chains":
        cmd_name = "generate_chains"
    
    # Execute the method
    print(f"Executing command: {cmd_name}")
    try:
        method = getattr(orchestrator, cmd_name)
        method(**cmd_args)
    except Exception as e:
        print(f"Error executing simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()