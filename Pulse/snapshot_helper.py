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
    """Efficiently parses only the LAST frame of a LAMMPS dump file."""
    try:
        file_size = os.path.getsize(dump_path)
        # Read the last 1MB of the file (should be enough for a frame of several thousand atoms)
        # If it's not enough, we can increase this or implement a multi-step seek.
        read_size = min(file_size, 1024 * 1024 * 5) # 5MB buffer
        
        with open(dump_path, 'rb') as f:
            f.seek(file_size - read_size)
            chunk = f.read(read_size).decode('utf-8', errors='ignore')
            
            # Find the last occurrence of "ITEM: ATOMS"
            atom_start = chunk.rfind("ITEM: ATOMS")
            if atom_start == -1:
                # If not found in the last 5MB, we might need to read more, but for now fallback
                return None
            
            # Extract column names
            lines = chunk[atom_start:].splitlines()
            cols = lines[0].strip().split()[2:]
            
            # The remaining lines are the atom data
            data_lines = lines[1:]
            
            # Parse into DataFrame
            from io import StringIO
            df = pd.read_csv(StringIO("\n".join(data_lines)), sep=r'\s+', names=cols, engine='c')
            return df
    except Exception as e:
        print(f"Error parsing last frame: {e}")
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
