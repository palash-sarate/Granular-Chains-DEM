import os
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Force headless backend to prevent Errno 5
import matplotlib.pyplot as plt
from pathlib import Path

import streamlit as st

def get_latest_dump(run_dir):
    """Finds the most recent dump file in the 'chain' directory."""
    chain_dir = os.path.join(run_dir, "chain")
    if not os.path.exists(chain_dir):
        return None
    
    dump_files = glob.glob(os.path.join(chain_dir, "*.dump"))
    if not dump_files:
        return None
    
    # Sort by step number in filename if possible
    try:
        dump_files.sort(key=lambda x: int(os.path.basename(x).replace("chain_", "").replace(".dump", "")))
    except:
        dump_files.sort(key=os.path.getmtime)
        
    return dump_files[-1]

def parse_dump_last_frame(dump_path):
    """Parses x, y coordinates from a LAMMPS dump file."""
    skip = 0
    cols = []
    try:
        with open(dump_path, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith("ITEM: ATOMS"):
                    cols = line.strip().split()[2:]
                    skip = i + 1
                    break
        
        # We only need x, y for a 2D preview. 
        # Check if they exist.
        needed = [c for c in cols if c in ['x', 'y', 'z', 'mol', 'diameter', 'type']]
        df = pd.read_csv(dump_path, skiprows=skip, sep=r'\s+', names=cols, usecols=needed, engine='c')
        return df
    except Exception as e:
        return None

def generate_snapshot(run_dir, output_dir="temp/snapshots", force=False):
    """Generates a PNG snapshot of the simulation state."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    run_name = os.path.basename(run_dir)
    output_path = os.path.join(output_dir, f"{run_name}.png")
    
    # Check if we should re-generate (if dump is newer than snapshot)
    dump_path = get_latest_dump(run_dir)
    if not dump_path:
        return None
    
    if not force and os.path.exists(output_path):
        if os.path.getmtime(output_path) > os.path.getmtime(dump_path):
            return output_path


    # Generate new snapshot
    df = parse_dump_last_frame(dump_path)
    if df is None or df.empty:
        return None
    
    # Set fixed Y limits for snapshot
    spacing = 0.2


    plt.figure(figsize=(6, 4))

    # Color by molecule if possible, else by type
    color_col = 'mol' if 'mol' in df.columns else ('type' if 'type' in df.columns else None)
    
    # Plot Y-Z plane (X pointing out)
    y_data = df['y']
    z_data = df['z'] if 'z' in df.columns else df['x'] # Fallback if 2D
    
    if color_col:
        plt.scatter(y_data, z_data, 
                    c=df[color_col], cmap='tab20', s=0.05, alpha=0.7)
    else:
        plt.scatter(y_data, z_data, 
                    color='blue', s=0.05, alpha=0.7)



    
    plt.axis('auto')
    plt.xlim(-spacing, spacing) # Set Y limits as requested (horizontal axis)



    plt.title(f"Preview: {run_name} (Y-Z plane)")
    plt.xlabel("Y")
    plt.ylabel("Z")

    plt.tight_layout()
    
    plt.savefig(output_path, dpi=100)
    plt.close()
    
    return output_path
