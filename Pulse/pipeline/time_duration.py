import os
import json
import pandas as pd

def run_time_duration(run_path, atoms_path, output_path):
    print(f"Calculating time duration taken to empty the hopper for: {run_path}")
    df = pd.read_parquet(atoms_path)
        
    meta_names = ["grid_metadata.json", "metadata.json", "sim_metadata.json"]
    meta = {}
    for name in meta_names:
        p = os.path.join(run_path, name)
        if os.path.exists(p):
            try:
                with open(p, 'r') as f_meta:
                    meta = json.load(f_meta)
                    break
            except:
                pass
                
    from Pulse.pulse_core import PBSManager
    lineage = PBSManager.load_lineage()
    run_info = lineage.get(run_path, {})
    
    N = run_info.get("N") or meta.get("N") or 4
    params = run_info.get("params", {}) or meta.get("params", {})
    dt = float(params.get("dt", 1e-6))
    
    orifice_width = params.get("orifice_width") or params.get("orifice")
    if not orifice_width and "geometry_vars" in params:
        g_vars = params["geometry_vars"]
        if isinstance(g_vars, list) and g_vars:
            orifice_width = g_vars[0].get("orifice_width")
    if not orifice_width:
        orifice_width = 0.05
        
    frequency = params.get("freq") or params.get("frequency") or 0.0
    
    df_reset = df.reset_index()
    in_hopper = df_reset[df_reset['z'] > 0]
    if in_hopper.empty:
        empty_timestep = 0
        time_taken = 0.0
    else:
        all_timesteps = sorted(df_reset['timestep'].unique())
        ts_in_hopper = set(in_hopper['timestep'].unique())
        
        empty_timestep = all_timesteps[-1]
        for ts in reversed(all_timesteps):
            if ts in ts_in_hopper:
                empty_timestep = ts
                break
        
        time_taken = float(empty_timestep * dt)
        
    results = {
        "run_name": os.path.basename(run_path),
        "N": int(N),
        "orifice_width": float(orifice_width),
        "frequency": float(frequency),
        "empty_timestep": int(empty_timestep),
        "time_taken": time_taken
    }
    
    df_out = pd.DataFrame([results])
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_out.to_parquet(output_path, index=False)
    print(f"Saved time duration results to {output_path}")
