import argparse
import sys
import os
import pandas as pd

# Add project root to sys.path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

def load_status(status_path):
    stages = ["load_atoms", "time_duration", "mass_flow_rate"]
    if os.path.exists(status_path):
        try:
            df = pd.read_parquet(status_path)
            df.set_index('stage', inplace=True)
            for s in stages:
                if s not in df.index:
                    df.loc[s] = ["not_started", "", ""]
            return df
        except:
            pass
            
    df = pd.DataFrame([
        {"stage": "load_atoms", "status": "not_started", "job_id": "", "error": ""},
        {"stage": "time_duration", "status": "not_started", "job_id": "", "error": ""},
        {"stage": "mass_flow_rate", "status": "not_started", "job_id": "", "error": ""}
    ]).set_index('stage')
    return df

def save_status(df, status_path):
    df.reset_index().to_parquet(status_path, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=str, required=True, help="Absolute path to the run directory")
    parser.add_argument("--stage", type=str, required=True, choices=["load_atoms", "time_duration", "mass_flow_rate"])
    args = parser.parse_args()
    
    run_dir = os.path.abspath(args.run)
    pipeline_dir = os.path.join(run_dir, "results_pipeline")
    os.makedirs(pipeline_dir, exist_ok=True)
    
    status_path = os.path.join(pipeline_dir, "pipeline_status.parquet")
    df_status = load_status(status_path)
    
    df_status.loc[args.stage, "status"] = "running"
    df_status.loc[args.stage, "error"] = ""
    save_status(df_status, status_path)
    
    try:
        if args.stage == "load_atoms":
            from Pulse.pipeline.load_atoms import run_load_atoms
            out_path = os.path.join(pipeline_dir, "load_atoms.parquet")
            run_load_atoms(run_dir, out_path)
            
        elif args.stage == "time_duration":
            from Pulse.pipeline.time_duration import run_time_duration
            atoms_path = os.path.join(pipeline_dir, "load_atoms.parquet")
            out_path = os.path.join(pipeline_dir, "time_duration.parquet")
            run_time_duration(run_dir, atoms_path, out_path)
            
        elif args.stage == "mass_flow_rate":
            from Pulse.pipeline.mass_flow_rate import run_mass_flow_rate
            atoms_path = os.path.join(pipeline_dir, "load_atoms.parquet")
            out_path = os.path.join(pipeline_dir, "mass_flow_rate.parquet")
            run_mass_flow_rate(run_dir, atoms_path, out_path)
            
        df_status.loc[args.stage, "status"] = "completed"
    except Exception as e:
        import traceback
        traceback.print_exc()
        df_status.loc[args.stage, "status"] = "failed"
        df_status.loc[args.stage, "error"] = str(e)
        
    save_status(df_status, status_path)

if __name__ == "__main__":
    main()
