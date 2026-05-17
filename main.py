import os
import sys
import argparse
import json

from utilities.env_loader import load_env
from simulation.generalized_engine import GeneralizedEngine

def main():
    env = load_env()
    parser = argparse.ArgumentParser(description="Generalized Pulse Simulation Engine")
    parser.add_argument("--lammps-exe", default=env.get("LAMMPS_EXECUTABLE", "lmp"), help="Path to LAMMPS executable")
    
    subparsers = parser.add_subparsers(dest="command", help="Simulation commands", required=True)

    # run_generalized
    p_gen = subparsers.add_parser("run_generalized", help="Run a generalized configuration-driven LAMMPS simulation stage")
    p_gen.add_argument("--schema", required=True, help="Path to simulation YAML or JSON schema file")
    p_gen.add_argument("--stage", required=True, help="Stage ID to run from the schema")
    p_gen.add_argument("--run-name", dest="run_name", help="Custom run directory name override")
    p_gen.add_argument("--parent-dir", dest="parent_dir", help="Parent directory (required for branch or resume modes)")
    p_gen.add_argument("--params", help="JSON string of parameter overrides")
    p_gen.add_argument("--num_procs", type=int, default=1, help="Number of processes (MPI)")
    p_gen.add_argument("--num_threads", type=int, default=1, help="Number of threads per process (OpenMP)")
    p_gen.add_argument("--no-kokkos", action="store_false", dest="use_kokkos", help="Disable Kokkos acceleration")
    p_gen.add_argument("--no-intel", action="store_false", dest="use_intel", help="Disable Intel acceleration")
    p_gen.add_argument("--inplace", action="store_true", help="Resume simulation directly in parent folder")

    args = parser.parse_args()
    
    cmd_name = args.command
    cmd_args = vars(args).copy()
    
    # Remove metadata
    cmd_args.pop("command", None)
    cmd_args.pop("lammps_exe", None)
    
    if cmd_name == "run_generalized":
        params_dict = {}
        if cmd_args.get("params"):
            try:
                params_dict = json.loads(cmd_args["params"])
            except Exception as e:
                print(f"Error parsing --params JSON: {e}")
                sys.exit(1)
        
        engine = GeneralizedEngine()
        
        # Override LAMMPS executable path from CLI if provided
        if args.lammps_exe:
            engine.lammps_exe = args.lammps_exe
            
        print(f"Executing generalized stage '{cmd_args['stage']}' using schema '{cmd_args['schema']}'")
        try:
            engine.run_stage(
                schema_path=cmd_args["schema"],
                stage_id=cmd_args["stage"],
                run_name=cmd_args["run_name"],
                parent_dir=cmd_args["parent_dir"],
                param_overrides=params_dict,
                num_procs=cmd_args["num_procs"],
                num_threads=cmd_args["num_threads"],
                use_kokkos=cmd_args["use_kokkos"],
                use_intel=cmd_args["use_intel"],
                inplace=cmd_args["inplace"]
            )
            print("Generalized stage execution completed successfully.")
            sys.exit(0)
        except Exception as e:
            print(f"Error executing generalized stage: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()