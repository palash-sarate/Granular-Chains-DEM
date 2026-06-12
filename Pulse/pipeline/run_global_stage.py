import argparse
import sys
import os
import pandas as pd

# Add project root to sys.path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

def load_status(status_path):
    stages = ["pool_time_durations"]
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
        {"stage": "pool_time_durations", "status": "not_started", "job_id": "", "error": ""}
    ]).set_index('stage')
    return df

def save_status(df, status_path):
    df.reset_index().to_parquet(status_path, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=str, required=True, choices=["pool_time_durations"])
    parser.add_argument("--runs", type=str, nargs="+", required=True, help="List of run directories to pool")
    args = parser.parse_args()
    
    global_dir = os.path.join(ROOT_DIR, "dumping_yard", "Results_Pipeline")
    os.makedirs(global_dir, exist_ok=True)
    
    status_path = os.path.join(global_dir, "pipeline_status.parquet")
    df_status = load_status(status_path)
    
    df_status.loc[args.stage, "status"] = "running"
    df_status.loc[args.stage, "error"] = ""
    save_status(df_status, status_path)
    
    try:
        if args.stage == "pool_time_durations":
            from Pulse.pipeline.pool_time_durations import run_pool_time_durations
            out_path = os.path.join(global_dir, "time_durations.parquet")
            run_pool_time_durations(args.runs, out_path)
            
        df_status.loc[args.stage, "status"] = "completed"
    except Exception as e:
        import traceback
        traceback.print_exc()
        df_status.loc[args.stage, "status"] = "failed"
        df_status.loc[args.stage, "error"] = str(e)
        
    save_status(df_status, status_path)

if __name__ == "__main__":
    main()
