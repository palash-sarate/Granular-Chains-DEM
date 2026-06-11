import os
import shutil
import time
import pandas as pd
import numpy as np
from analysis.data_manager import SimulationData

def run_bench():
    sim_dir = 'dumping_yard/Hopper_Fill/Grid_Fill_1H_S374952'
    cache_dir = os.path.join(sim_dir, 'sim_cache')
    
    # Clean the cache folder to measure raw speed (first load)
    if os.path.exists(cache_dir):
        print(f"Cleaning existing cache in {cache_dir}...")
        shutil.rmtree(cache_dir)
        
    print(f"\n--- Processing & Benchmarking {sim_dir} ---")
    
    # 1. Initialize simulation data source
    sim = SimulationData(sim_dir)
    
    # 2. Benchmark raw parsing + calculation + cache writing (First load)
    t0 = time.time()
    data_dict = sim.load_data(force_reload=True)
    t1 = time.time()
    time_raw = t1 - t0
    
    df_atoms = data_dict['atoms']
    df_bonds = data_dict['bonds']
    df_angles = data_dict['angles']
    
    print(f"First Load (Raw parse, calculation, and cache write): {time_raw*1000:.2f} ms")
    print(f"  Atoms shape: {df_atoms.shape if df_atoms is not None else 'None'}")
    print(f"  Bonds shape: {df_bonds.shape if df_bonds is not None else 'None'}")
    print(f"  Angles shape: {df_angles.shape if df_angles is not None else 'None'}")
    
    # Verify values are valid
    if df_bonds is not None and not df_bonds.empty:
        print(f"  Sample Bond lengths: min={df_bonds['dist'].min():.6f}, max={df_bonds['dist'].max():.6f}, mean={df_bonds['dist'].mean():.6f}")
    if df_angles is not None and not df_angles.empty:
        # angles are calculated in radians
        angles_deg = df_angles['theta'] * 180.0 / np.pi
        print(f"  Sample Bond angles (deg): min={angles_deg.min():.2f}, max={angles_deg.max():.2f}, mean={angles_deg.mean():.2f}")
        
    # 3. Benchmark cached loading (Second load)
    t0 = time.time()
    data_dict_cached = sim.load_batch(0, force_reload=False)
    t1 = time.time()
    time_cached = t1 - t0
    print(f"Second Load (Cached retrieval from pickle): {time_cached*1000:.2f} ms")
    print(f"Speedup from cache: {time_raw / time_cached:.2f}x")
    
    # Check that cache folder has been populated
    cache_files = os.listdir(cache_dir)
    print(f"Cache files written: {cache_files}")

if __name__ == '__main__':
    run_bench()
