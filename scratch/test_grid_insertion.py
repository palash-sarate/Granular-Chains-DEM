import os
import sys
import numpy as np
from pathlib import Path

# Add project root to sys.path
ROOT_DIR = Path(__file__).parent.parent.absolute()
sys.path.append(str(ROOT_DIR))

from simulation.runner import SimulationRunner
from simulation.grid_hopper_manager import GridHopperManager
from simulation.hopper_manager import HopperManager

def main():
    # 1. Setup paths
    source_dir = ROOT_DIR / "chain_data" / "relaxed_2D_x" / "N4"
    temp_dir = ROOT_DIR / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    mol_dir = temp_dir / "molecules"
    mol_dir.mkdir(exist_ok=True)
    
    # 2. Setup Managers
    # We don't need a real runner just to generate the file, but we need it for the dummy run
    runner = SimulationRunner(lammps_executable="lmp") 
    h_manager = HopperManager(runner)
    grid_manager = GridHopperManager(runner)
    
    print(f"--- Preparing molecule templates from {source_dir} ---")
    # This generates molecule files m1.mol, m2.mol... and an include file
    mol_inc_path = h_manager.prepare_molecules(str(source_dir), str(mol_dir))
    n_templates = len(list(source_dir.glob("*.data")))
    print(f"Generated templates for {n_templates} molecules.")

    # 3. Setup Metadata for multiple hoppers
    n_hoppers = 4
    spacing = 2.0
    
    # Generate metadata for a 2x2 grid (or similar)
    metadata = {}
    nx_h = 2
    for i in range(n_hoppers):
        ix = i % nx_h
        iy = i // nx_h
        metadata[i] = {'offset': np.array([ix * spacing, iy * spacing, 0.0])}
    
    print(f"--- Generating insertions.inc for {n_hoppers} hoppers, 6000 chains each (N=4) ---")
    n_fill = 6000
    seed = 12345
    N = 4
    
    insertion_path_str, z_max = grid_manager._generate_grid_insertion_file(
        temp_dir, n_hoppers, n_fill, n_templates, seed, N, spacing, metadata, mode="2D_stacked"
    )
    print(f"Generated {insertion_path_str}, z_max = {z_max}")

    # 4. Create Dummy LAMMPS Script to generate a data file for visualization
    dummy_in = temp_dir / "dummy_visualize.in"
    
    # Calculate box to fit all hoppers
    x_max_h = (n_hoppers - 1) % nx_h * spacing + 1.0
    y_max_h = (n_hoppers - 1) // nx_h * spacing + 1.0
    
    with open(dummy_in, 'w') as f:
        f.write(f"""
units lj
atom_style hybrid molecular sphere
boundary p p p

# Box expanded to fit all hoppers
region world block -1.0 {x_max_h} -1.0 {y_max_h} 0.0 {z_max + 0.5}
create_box 1 world bond/types 1 angle/types 1 extra/bond/per/atom 5 extra/angle/per/atom 5 extra/special/per/atom 10

# Load molecule templates
include {mol_inc_path}

# Perform insertions
include {insertion_path_str}

mass 1 1.0
set type 1 diameter 0.002

# Write out for visualization
write_data {temp_dir / 'inserted_state.data'}
""")

    print(f"--- Running dummy LAMMPS simulation to create data file ---")
    cmd = f"lmp -in {dummy_in}"
    import subprocess
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"SUCCESS: Data file created at {temp_dir / 'inserted_state.data'}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR running LAMMPS:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")

    # 5. Cleanup
    print(f"--- Cleaning up temporary molecule files ---")
    import shutil
    shutil.rmtree(mol_dir)

if __name__ == "__main__":
    main()
