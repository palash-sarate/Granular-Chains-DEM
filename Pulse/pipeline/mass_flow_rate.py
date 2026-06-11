import os
import pandas as pd

def run_mass_flow_rate(run_path, atoms_path, output_path):
    print(f"Calculating mass flow rate details for: {run_path}")
    df = pd.read_parquet(atoms_path)
    
    df_reset = df.reset_index()
    
    records = []
    for ts, group in df_reset.groupby('timestep'):
        in_hopper = group[group['z'] > 0]
        num_beads = len(in_hopper)
        total_mass = float(in_hopper['mass'].sum()) if 'mass' in in_hopper.columns else 0.0
        records.append({
            "timestep": int(ts),
            "mass_in_hopper": total_mass,
            "number_of_beads_in_hopper": int(num_beads)
        })
        
    df_out = pd.DataFrame(records)
    df_out.sort_values(by='timestep', inplace=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_out.to_parquet(output_path, index=False)
    print(f"Saved mass flow rate data to {output_path}")
