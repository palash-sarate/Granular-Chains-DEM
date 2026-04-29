import os
import sys
import argparse
from typing import List, Optional
import json

from simulation.orchestrator import SimulationOrchestrator

def main():
    parser = argparse.ArgumentParser(description="LAMMPS Chain Simulation CLI")
    parser.add_argument("--lammps-exe", default="lmp", help="Path to LAMMPS executable (default: lmp)")
    
    subparsers = parser.add_subparsers(dest="command", help="Simulation commands", required=True)

    # 1. run_hopper_fill
    p_fill = subparsers.add_parser("run_hopper_fill", help="Pre-fill a hopper with relaxed molecular chains")
    p_fill.add_argument("--source_dir", default="chain_data/relaxed/N4")
    p_fill.add_argument("--fill_template", default="in.hopper_fill")
    p_fill.add_argument("--n_fill", type=int, default=10)
    p_fill.add_argument("--relax_steps", type=int, default=100000)
    p_fill.add_argument("--run_name")
    p_fill.add_argument("--seed", type=int)
    p_fill.add_argument("--dt", type=float, default=1e-6)
    p_fill.add_argument("--mol_dir", default="chain_data/molecules_temp")
    p_fill.add_argument("--setup_inc", default="simulation_geometries/2D_hopper.inc")
    p_fill.add_argument("--dump_inc", default="simulation_templates/default_dump.inc")
    p_fill.add_argument("--viscosity", type=float, default=0.001)
    p_fill.add_argument("--N", type=int, default=4)
    p_fill.add_argument("--outdir")
    p_fill.add_argument("--num_procs", type=int)
    p_fill.add_argument("--num_threads", type=int, default=1)
    p_fill.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_fill.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 2. resume_hopper_fill
    p_resume_h = subparsers.add_parser("resume_hopper_fill", help="Resume hopper filling from a restart file")
    p_resume_h.add_argument("--restart_path", required=True)
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
    p_resume_h.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_resume_h.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 3. run_flop_simulation
    p_flop = subparsers.add_parser("run_flop_simulation", help="Run a single chain-flop mobility simulation")
    p_flop.add_argument("--N", type=int, default=4)
    p_flop.add_argument("--run_steps", type=int, default=50000)
    p_flop.add_argument("--viscosity", type=float, default=0.001)
    p_flop.add_argument("--dt", type=float, default=1e-6)
    p_flop.add_argument("--num_procs", type=int)
    p_flop.add_argument("--num_threads", type=int, default=1)
    p_flop.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_flop.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 4. resume_flop_simulation
    p_resume_f = subparsers.add_parser("resume_flop_simulation", help="Resume a chain-flop simulation")
    p_resume_f.add_argument("--N", type=int, default=4)
    p_resume_f.add_argument("--run_steps", type=int, default=50000)
    p_resume_f.add_argument("--viscosity", type=float, default=0.001)
    p_resume_f.add_argument("--dt", type=float, default=1e-6)
    p_resume_f.add_argument("--resume_token", required=True)
    p_resume_f.add_argument("--num_procs", type=int)
    p_resume_f.add_argument("--num_threads", type=int, default=1)
    p_resume_f.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_resume_f.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 5. generate_relaxed_library
    p_lib = subparsers.add_parser("generate_relaxed_library", help="Generate relaxed chain states library")
    p_lib.add_argument("--n_beads", default="4")
    p_lib.add_argument("--n_states", type=int, default=10)
    p_lib.add_argument("--forced", action="store_true")
    p_lib.add_argument("--template", default="in.relax_3d_gen")
    p_lib.add_argument("--output_dir", default="chain_data/relaxed")
    p_lib.add_argument("--run_name_prefix")
    p_lib.add_argument("--simulation_name", default="Relax_Library_Gen")
    p_lib.add_argument("--dump_inc", default="simulation_templates/default_dump.inc")
    p_lib.add_argument("--inParallel", type=int, default=1, dest="n_parallel")
    p_lib.add_argument("--num_procs", type=int)
    p_lib.add_argument("--num_threads", type=int, default=1)
    p_lib.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_lib.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 6. generate_linear_chains
    p_chains = subparsers.add_parser("generate_linear_chains", help="Generate initial linear chain data files")
    p_chains.add_argument("--Ns", default="4,6,8")
    p_chains.add_argument("--orientation", default="horz", choices=["horz", "vert"])
    p_chains.add_argument("--output_dir_name", default="linear_x")

    # 7. run_flop_batch
    p_batch = subparsers.add_parser("run_flop_batch", help="Run a batch of flop simulations")
    p_batch.add_argument("--Ns", default="4,6")
    p_batch.add_argument("--run_steps", default="50000")
    p_batch.add_argument("--viscosities", default="0.001")
    p_batch.add_argument("--dt", type=float, default=1e-6)
    p_batch.add_argument("--num_procs", type=int)
    p_batch.add_argument("--num_threads", type=int, default=1)
    p_batch.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_batch.add_argument("--no-intel", action="store_false", dest="use_intel")

    # 8. run_grid_batch_relaxation
    p_grid = subparsers.add_parser("run_grid_batch_relaxation", help="Run a super-simulation for batch relaxation")
    p_grid.add_argument("--n_beads", default="4,8,12,24,48,100")
    p_grid.add_argument("--n_states", type=int, default=50)
    p_grid.add_argument("--output_dir", default="chain_data/relaxed_grid")
    p_grid.add_argument("--spacing", type=float, default=0.5)
    p_grid.add_argument("--dt", type=float, default=1e-6)
    p_grid.add_argument("--num_procs", type=int, default=1)
    p_grid.add_argument("--num_threads", type=int, default=1)
    p_grid.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_grid.add_argument("--no-intel", action="store_false", dest="use_intel")
    p_grid.add_argument("--lepton_file", default="simulation_templates/lepton.inc")
    p_grid.add_argument("--dump_file", default="simulation_templates/quiet_dump.inc")
    p_grid.add_argument("--simulation", default="Grid_Relaxation_Batch")
    p_grid.add_argument("--seed", type=int, default=None)
    p_grid.add_argument("--motion_steps", type=int, default=500000)
    p_grid.add_argument("--explore_steps", type=int, default=200000)
    p_grid.add_argument("--viscous_relax_steps", type=int, default=300000)
    p_grid.add_argument("--temperature", type=float, default=1e15)

    # 9. run_grid_hopper_filling
    p_grid_h = subparsers.add_parser("run_grid_hopper_filling", help="Run a super-simulation for batch hopper filling")
    p_grid_h.add_argument("--n_hoppers", type=int, default=4)
    p_grid_h.add_argument("--n_fill", default="10")
    p_grid_h.add_argument("--N", default="4")
    p_grid_h.add_argument("--spacing", type=float, default=2.0)
    p_grid_h.add_argument("--relax_steps", type=int, default=500000)
    p_grid_h.add_argument("--dt", type=float, default=1e-6)
    p_grid_h.add_argument("--output_dir", default="chain_data/grid_filled")
    p_grid_h.add_argument("--source_dir", help="Directory containing relaxed chains")
    p_grid_h.add_argument("--hopper_template_data", default="simulation_geometries/2D_hopper.inc")
    p_grid_h.add_argument("--lepton_file", default="simulation_templates/lepton.inc")
    p_grid_h.add_argument("--dump_file", default="simulation_templates/quiet_dump.inc")
    p_grid_h.add_argument("--viscosity", type=float, default=0.001)
    p_grid_h.add_argument("--num_procs", type=int, default=1)
    p_grid_h.add_argument("--num_threads", type=int, default=1)
    p_grid_h.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_grid_h.add_argument("--mode", type=str, default="2D_stacked", choices=["2D_stacked", "2D_worst_case", "3D_grid"],
                            help="Pouring mode: 2D_stacked (smart bbox), 2D_worst_case (linear spacing), or 3D_grid")
    p_grid_h.add_argument("--simulation", default="Grid_Hopper_Filling")
    p_grid_h.add_argument("--template", default="in.grid_hopper_fill")
    p_grid_h.add_argument("--no-vtk", action="store_false", dest="generate_vtk", help="Skip VTK mesh generation")
    p_grid_h.set_defaults(generate_vtk=True)
    p_grid_h.add_argument("--geometry_vars", type=str, help="JSON string for per-hopper geometry variables")

    # 10. resume_grid_hopper_filling
    p_res_grid = subparsers.add_parser("resume_grid_hopper_filling", help="Resume a super-simulation for batch hopper filling")
    p_res_grid.add_argument("--restart_path", required=True)
    p_res_grid.add_argument("--relax_steps", type=int, default=500000)
    p_res_grid.add_argument("--dt", type=float, default=1e-6)
    p_res_grid.add_argument("--output_dir", default="chain_data/grid_filled")
    p_res_grid.add_argument("--lepton_file", default="simulation_templates/lepton.inc")
    p_res_grid.add_argument("--dump_file", default="simulation_templates/quiet_dump.inc")
    p_res_grid.add_argument("--viscosity", type=float, default=0.001)
    p_res_grid.add_argument("--num_procs", type=int, default=1)
    p_res_grid.add_argument("--num_threads", type=int, default=1)
    p_res_grid.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_res_grid.add_argument("--template", default="in.grid_hopper_fill_resume")
    p_res_grid.add_argument("--seed", type=int, default=None)

    # 11. run_grid_hopper_flow
    p_flow = subparsers.add_parser("run_grid_hopper_flow", help="Transition filled hoppers to oscillatory flow")
    p_flow.add_argument("--source_dir", required=True, help="Directory of a Grid_Hopper_Filling run")
    p_flow.add_argument("--run_steps", type=int, default=1000000)
    p_flow.add_argument("--freq", type=str, default="10.0", help="Freq (single float or JSON list)")
    p_flow.add_argument("--amp", type=str, default="0.01", help="Amp (single float or JSON list)")
    p_flow.add_argument("--osc_dir", choices=['x','y','z'], default='z')
    p_flow.add_argument("--dt", type=float, default=1e-6)
    p_flow.add_argument("--output_dir", default="chain_data/grid_flow")
    p_flow.add_argument("--num_procs", type=int, default=1)
    p_flow.add_argument("--num_threads", type=int, default=1)
    p_flow.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_flow.add_argument("--lepton_file", default="simulation_templates/lepton.inc")
    p_flow.add_argument("--dump_file", default="simulation_templates/quiet_dump.inc")
    p_flow.add_argument("--viscosity", type=float, default=0.001)
    p_flow.add_argument("--simulation", default="Grid_Hopper_Flow")
    p_flow.add_argument("--template", default="in.grid_hopper_flow")
    p_flow.add_argument("--seed", type=int, default=None)

    # 12. resume_grid_hopper_flow
    p_res_flow = subparsers.add_parser("resume_grid_hopper_flow", help="Resume a grid flow simulation")
    p_res_flow.add_argument("--restart_path", required=True)
    p_res_flow.add_argument("--run_steps", type=int, default=1000000)
    p_res_flow.add_argument("--dt", type=float, default=1e-6)
    p_res_flow.add_argument("--num_procs", type=int, default=1)
    p_res_flow.add_argument("--num_threads", type=int, default=1)
    p_res_flow.add_argument("--no-kokkos", action="store_false", dest="use_kokkos")
    p_res_flow.add_argument("--lepton_file", default="simulation_templates/lepton.inc")
    p_res_flow.add_argument("--dump_file", default="simulation_templates/quiet_dump.inc")
    p_res_flow.add_argument("--viscosity", type=float, default=0.001)
    p_res_flow.add_argument("--template", default="in.grid_hopper_flow_resume")
    p_res_flow.add_argument("--seed", type=int, default=None)

    args = parser.parse_args()
    
    orchestrator = SimulationOrchestrator(lammps_executable=args.lammps_exe)
    
    cmd_name = args.command
    cmd_args = vars(args).copy()
    
    # Remove metadata
    cmd_args.pop("command", None)
    cmd_args.pop("lammps_exe", None)
    
    # Special parsing
    if cmd_name == "run_flop_batch":
        cmd_args["Ns"] = [int(n.strip()) for n in str(cmd_args["Ns"]).split(",")]
        cmd_args["run_steps"] = [int(s.strip()) for s in str(cmd_args["run_steps"]).split(",")]
        cmd_args["viscosities"] = [float(v.strip()) for v in str(cmd_args["viscosities"]).split(",")]
    
    if cmd_name == "run_grid_hopper_filling":
        # Parse N and n_fill into lists
        if isinstance(cmd_args["N"], str):
            cmd_args["N"] = [int(n.strip()) for n in cmd_args["N"].split(",")]
        else:
            cmd_args["N"] = [cmd_args["N"]]
            
        if isinstance(cmd_args["n_fill"], str):
            cmd_args["n_fill"] = [int(n.strip()) for n in cmd_args["n_fill"].split(",")]
        else:
            cmd_args["n_fill"] = [cmd_args["n_fill"]]

        # Handle geometry_vars JSON
        if cmd_args.get("geometry_vars"):
            try:
                cmd_args["geometry_vars"] = json.loads(cmd_args["geometry_vars"])
            except Exception as e:
                print(f"Error parsing --geometry_vars JSON: {e}")
                sys.exit(1)
        else:
            cmd_args["geometry_vars"] = {}

    if cmd_name == "run_grid_hopper_flow":
        # Parse freq and amp
        for key in ["freq", "amp"]:
            val = cmd_args[key]
            if isinstance(val, str):
                if "[" in val:
                    try:
                        cmd_args[key] = json.loads(val)
                    except:
                        cmd_args[key] = [float(x.strip()) for x in val.strip("[]").split(",")]
                elif "," in val:
                    cmd_args[key] = [float(x.strip()) for x in val.split(",")]
                else:
                    cmd_args[key] = float(val)

    if cmd_name == "generate_linear_chains":
        cmd_name = "generate_chains"
    
    # Execute
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