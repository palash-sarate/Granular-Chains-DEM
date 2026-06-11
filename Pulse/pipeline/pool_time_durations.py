import os
import pandas as pd

def run_pool_time_durations(run_dirs, output_path):
    print(f"Pooling time duration results from {len(run_dirs)} runs...")
    rows = []
    for r_dir in run_dirs:
        stage_path = os.path.join(r_dir, "results_pipeline", "time_duration.parquet")
        if os.path.exists(stage_path):
            try:
                df = pd.read_parquet(stage_path)
                if not df.empty:
                     rows.append(df)
            except Exception as e:
                print(f"Error reading {stage_path}: {e}")
                
    if not rows:
        df_out = pd.DataFrame(columns=["run_name", "N", "orifice_width", "frequency", "empty_timestep", "time_taken"])
    else:
        df_out = pd.concat(rows, ignore_index=True)
        
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_out.to_parquet(output_path, index=False)
    print(f"Saved global pooled durations to {output_path}")
