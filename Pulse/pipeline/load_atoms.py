import os
import sys
import pandas as pd
from analysis.data_manager import SimulationData

def run_load_atoms(run_path, output_path):
    print(f"Loading and caching raw simulation data for: {run_path}")
    sim_data = SimulationData(run_path)
    sim_data.load_metadata()
    
    num_batches = sim_data.get_batch_count()
    if num_batches == 0:
        raise ValueError(f"No simulation dump files found in {run_path}")
        
    all_atoms_list = []
    
    for b_idx in range(num_batches):
        print(f"Loading batch {b_idx + 1} of {num_batches}...")
        batch = sim_data.load_batch(b_idx)
        df_atoms = batch.get('atoms')
        if df_atoms is not None and not df_atoms.empty:
            df_reset = df_atoms.reset_index()
            cols_to_keep = [c for c in ['timestep', 'id', 'mol', 'z', 'mass'] if c in df_reset.columns]
            df_sub = df_reset[cols_to_keep].copy()
            all_atoms_list.append(df_sub)
            
    if not all_atoms_list:
        raise ValueError("No atom frames were parsed successfully.")
        
    df_all = pd.concat(all_atoms_list)
    for col in ['id', 'mol', 'z', 'mass', 'timestep']:
        if col in df_all.columns:
            df_all[col] = pd.to_numeric(df_all[col], errors='coerce')
            
    if 'mass' not in df_all.columns:
        df_all['mass'] = 7.1e-6
        
    df_all.dropna(subset=['id', 'z', 'timestep'], inplace=True)
    df_all.set_index(['timestep', 'id'], inplace=True)
    df_all.sort_index(inplace=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_all.to_parquet(output_path)
    print(f"Saved simplified atoms coordinates to: {output_path}")
