import os
import sys
import argparse
import json
from pathlib import Path

from utilities.env_loader import load_env
from simulation.generalized_engine import GeneralizedEngine, SchemaValidationError

def run_interactive_wizard(lammps_exe: str):
    """
    A beautiful, colorized, step-by-step interactive CLI wizard to run simulations.
    Designed in vanilla Python to avoid external dependencies while providing a premium experience.
    """
    # ANSI escape codes for beautiful styling
    c_blue = "\033[94m"
    c_green = "\033[92m"
    c_yellow = "\033[93m"
    c_red = "\033[91m"
    c_bold = "\033[1m"
    c_end = "\033[0m"

    print(f"\n{c_blue}{c_bold}===================================================={c_end}")
    print(f"{c_blue}{c_bold}       ⚡ WELCOME TO THE PULSE SIMULATION WIZARD     {c_end}")
    print(f"{c_blue}{c_bold}===================================================={c_end}\n")

    # 1. Scan for available schemas
    schema_dir = Path("simulation_schemas")
    if not schema_dir.exists():
        print(f"{c_red}Error: 'simulation_schemas' folder not found at project root.{c_end}")
        return

    schemas = sorted([f for f in schema_dir.glob("*") if f.suffix.lower() in [".yaml", ".yml", ".json"]])
    if not schemas:
        print(f"{c_yellow}No schemas found in 'simulation_schemas' folder.{c_end}")
        return

    print(f"{c_bold}Available Simulation Schemas:{c_end}")
    for idx, schema_path in enumerate(schemas, 1):
        print(f"  {c_blue}[{idx}]{c_end} {schema_path.name}")

    while True:
        try:
            choice = input(f"\nSelect a schema number (1-{len(schemas)}): ").strip()
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(schemas):
                selected_schema = schemas[choice_idx]
                break
            print(f"{c_red}Invalid choice. Please select a valid number.{c_end}")
        except ValueError:
            print(f"{c_red}Please enter an integer number.{c_end}")

    print(f"\nLoading {c_green}{selected_schema.name}{c_end}...")
    engine = GeneralizedEngine()
    engine.lammps_exe = lammps_exe
    
    try:
        schema = engine.load_schema(str(selected_schema))
    except Exception as e:
        print(f"{c_red}Failed to load schema: {e}{c_end}")
        return

    sim_type = schema.get("simulation_type", "Generalized_Simulation")
    print(f"Simulation Type: {c_bold}{sim_type}{c_end}")

    # 2. Normalize and list stages
    stages_raw = schema.get("stages", [])
    stages_list = []
    if isinstance(stages_raw, dict):
        for k, v in stages_raw.items():
            if isinstance(v, dict):
                v_copy = v.copy()
                v_copy["id"] = k
                stages_list.append(v_copy)
    elif isinstance(stages_raw, list):
        stages_list = stages_raw

    if not stages_list:
        print(f"{c_red}No stages defined in this schema.{c_end}")
        return

    print(f"\n{c_bold}Available Stages:{c_end}")
    for idx, stage in enumerate(stages_list, 1):
        s_id = stage.get("id", f"stage_{idx}")
        s_mode = stage.get("mode", "start")
        s_desc = stage.get("name", "No description")
        print(f"  {c_blue}[{idx}]{c_end} {c_bold}{s_id}{c_end} ({s_mode}) - {s_desc}")

    while True:
        try:
            choice = input(f"\nSelect a stage number (1-{len(stages_list)}): ").strip()
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(stages_list):
                selected_stage_dict = stages_list[choice_idx]
                selected_stage = selected_stage_dict.get("id")
                break
            print(f"{c_red}Invalid choice. Please select a valid number.{c_end}")
        except ValueError:
            print(f"{c_red}Please enter an integer number.{c_end}")

    # 3. Handle parent directory for resume/branch modes
    parent_dir = None
    mode = selected_stage_dict.get("mode", "start")
    if mode in ["resume", "branch"] or selected_stage_dict.get("parent_stage"):
        print(f"\n{c_yellow}This stage is in '{mode}' mode and requires a parent run directory.{c_end}")
        # Look in dumping_yard/simulation_type/ for candidates
        sim_dump_dir = Path(engine.dumping_yard) / sim_type
        candidates = []
        if sim_dump_dir.exists():
            candidates = sorted(
                [d for d in sim_dump_dir.iterdir() if d.is_dir() and (d / "metadata.json").exists()],
                key=lambda x: x.stat().st_mtime,
                reverse=True
            )
        
        if candidates:
            print(f"Recent compatible parent runs:")
            for idx, cand in enumerate(candidates[:5], 1):
                print(f"  [{idx}] {cand.name}")
            print(f"  [o] Other (enter path manually)")
            
            cand_choice = input(f"\nSelect a parent run or press Enter for the most recent [{candidates[0].name}]: ").strip()
            if not cand_choice:
                parent_dir = str(candidates[0])
            elif cand_choice.lower() == "o":
                parent_dir = input("Enter full path to parent run folder: ").strip()
            else:
                try:
                    c_idx = int(cand_choice) - 1
                    if 0 <= c_idx < len(candidates):
                        parent_dir = str(candidates[c_idx])
                    else:
                        parent_dir = input("Enter full path to parent run folder: ").strip()
                except ValueError:
                    parent_dir = input("Enter full path to parent run folder: ").strip()
        else:
            parent_dir = input("Enter full path to parent run folder: ").strip()

    # 4. Prompt for parameter overrides
    param_overrides = {}
    
    # Collect available overrides from stage parameters
    stage_params = selected_stage_dict.get("parameters", selected_stage_dict.get("params", []))
    available_overrides = []
    if isinstance(stage_params, dict):
        for k, v in stage_params.items():
            available_overrides.append({"name": k, "default": v})
    elif isinstance(stage_params, list):
        for p in stage_params:
            if isinstance(p, dict) and "name" in p:
                available_overrides.append({"name": p["name"], "default": p.get("default")})

    if available_overrides:
        print(f"\n{c_bold}Configure Stage Overrides (Press Enter to keep default):{c_end}")
        for override in available_overrides:
            val = input(f"  {override['name']} [default: {override['default']}]: ").strip()
            if val:
                # Try parsing as float or int
                try:
                    if "." in val:
                        param_overrides[override["name"]] = float(val)
                    else:
                        param_overrides[override["name"]] = int(val)
                except ValueError:
                    param_overrides[override["name"]] = val

    # 5. Configure hardware counts
    print(f"\n{c_bold}Configure Hardware Scaling:{c_end}")
    
    num_procs = 1
    np_val = input("  Number of MPI Processes (default 1): ").strip()
    if np_val:
        try: num_procs = int(np_val)
        except: pass
        
    num_threads = 1
    nt_val = input("  OpenMP Threads per Process (default 1): ").strip()
    if nt_val:
        try: num_threads = int(nt_val)
        except: pass

    use_kokkos = True
    kk_val = input("  Enable KOKKOS Acceleration (Y/n, default Y): ").strip()
    if kk_val.lower() == "n":
        use_kokkos = False

    use_intel = True
    in_val = input("  Enable INTEL Acceleration (Y/n, default Y): ").strip()
    if in_val.lower() == "n":
        use_intel = False

    # Custom run name
    run_name = None
    rn_val = input("\nCustom Run Name Override (Press Enter for default naming): ").strip()
    if rn_val:
        run_name = rn_val

    # 6. Confirm and run!
    print(f"\n{c_green}{c_bold}Configuration Confirmed!{c_end}")
    print(f"  Schema: {selected_schema.name}")
    print(f"  Stage:  {selected_stage}")
    print(f"  Mode:   {mode}")
    if parent_dir:
        print(f"  Parent: {parent_dir}")
    if param_overrides:
        print(f"  Params: {param_overrides}")
    print(f"  Scaling: MPI={num_procs}, OpenMP={num_threads}, Kokkos={use_kokkos}, Intel={use_intel}")

    launch = input(f"\n{c_bold}Launch simulation now? (Y/n, default Y): {c_end}").strip()
    if launch.lower() == "n":
        print(f"{c_yellow}Launch cancelled.{c_end}")
        return

    print(f"\n{c_blue}Invoking Generalized Engine...{c_end}\n")
    try:
        engine.run_stage(
            schema_path=str(selected_schema),
            stage_id=selected_stage,
            run_name=run_name,
            parent_dir=parent_dir,
            param_overrides=param_overrides,
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel,
            inplace=False
        )
        print(f"\n{c_green}{c_bold}✓ Simulation Stage Completed Successfully!{c_end}\n")
    except SchemaValidationError as sve:
        print(f"{c_red}{sve}{c_end}")
    except Exception as e:
        print(f"\n{c_red}❌ Error during execution: {e}{c_end}\n")


def main():
    env = load_env()
    parser = argparse.ArgumentParser(description="Generalized Pulse Simulation Engine")
    parser.add_argument("--lammps-exe", default=env.get("LAMMPS_EXECUTABLE", "lmp"), help="Path to LAMMPS executable")
    
    # Subcommands
    subparsers = parser.add_subparsers(dest="command", help="Simulation commands")

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

    # interactive subcommand
    subparsers.add_parser("interactive", help="Start the interactive step-by-step terminal wizard")

    args = parser.parse_args()
    
    # Launch interactive wizard if no command is specified, or 'interactive' is selected
    if not args.command or args.command == "interactive":
        run_interactive_wizard(args.lammps_exe)
        sys.exit(0)

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
        except SchemaValidationError as sve:
            print(sve)
            sys.exit(1)
        except Exception as e:
            print(f"Error executing generalized stage: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()