import argparse
import sys
import os
import pandas as pd

# Add project root to sys.path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pipeline.pipeline_engine import PipelineEngine

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=str, required=True, help="Absolute path to the run directory")
    parser.add_argument("--stage", type=str, required=True, choices=["load_data", "bond_angle_calc", "time_duration", "mass_flow_rate"])
    args = parser.parse_args()
    
    run_dir = os.path.abspath(args.run)
    
    # Check if stage is applicable to this run type
    applicable = PipelineEngine.get_applicable_stages(run_dir)
    if args.stage not in applicable:
        print(f"Error: Stage '{args.stage}' is not applicable for this run type.")
        sys.exit(1)
        
    pipeline_dir = os.path.join(run_dir, "results_pipeline")
    os.makedirs(pipeline_dir, exist_ok=True)
    
    df_status = PipelineEngine.load_status(run_dir)
    
    df_status.loc[args.stage, "status"] = "running"
    df_status.loc[args.stage, "error"] = ""
    PipelineEngine.save_status(df_status, run_dir)
    
    try:
        if args.stage == "load_data":
            from Pulse.pipeline.load_data import run_load_data
            out_path = os.path.join(pipeline_dir, "load_data.parquet")
            run_load_data(run_dir, out_path)
            
        elif args.stage == "bond_angle_calc":
            from Pulse.pipeline.bond_angle_calc import run_bond_angle_calc
            atoms_path = os.path.join(pipeline_dir, "load_data.parquet")
            bonds_path = os.path.join(pipeline_dir, "bonds.parquet")
            angles_path = os.path.join(pipeline_dir, "angles.parquet")
            run_bond_angle_calc(run_dir, atoms_path, bonds_path, angles_path)
            
        elif args.stage == "time_duration":
            from Pulse.pipeline.time_duration import run_time_duration
            atoms_path = os.path.join(pipeline_dir, "load_data.parquet")
            out_path = os.path.join(pipeline_dir, "time_duration.parquet")
            run_time_duration(run_dir, atoms_path, out_path)
            
        elif args.stage == "mass_flow_rate":
            from Pulse.pipeline.mass_flow_rate import run_mass_flow_rate
            atoms_path = os.path.join(pipeline_dir, "load_data.parquet")
            out_path = os.path.join(pipeline_dir, "mass_flow_rate.parquet")
            run_mass_flow_rate(run_dir, atoms_path, out_path)
            
        df_status.loc[args.stage, "status"] = "completed"
    except Exception as e:
        import traceback
        traceback.print_exc()
        df_status.loc[args.stage, "status"] = "failed"
        df_status.loc[args.stage, "error"] = str(e)
        
    PipelineEngine.save_status(df_status, run_dir)

if __name__ == "__main__":
    main()
