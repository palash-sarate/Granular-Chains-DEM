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
    temp_dir = ROOT_DIR / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    mol_dir = temp_dir / "molecules"
    mol_dir.mkdir(exist_ok=True)
    
    # 2. Setup Managers
    # We don't need a real runner just to generate the file, but we need it for the dummy run
    runner = SimulationRunner(lammps_executable="lmp") 
    h_manager = HopperManager(runner)
    grid_manager = GridHopperManager(runner)

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
    
    N_list = [4, 12, 24, 48]
    n_fill_list = [14400 // N for N in N_list]
    
    from simulation.molecule_converter import convert_data_to_molecule
    mol_ranges = {}
    mol_bboxes = {}
    combined_inc_lines = []
    
    current_mol_id = 1
    
    for N in N_list:
        source_dir = ROOT_DIR / "chain_data" / "relaxed_2D_x" / f"N{N}"
        print(f"--- Preparing molecule templates from {source_dir} (N={N}) ---")
        
        data_files = list(source_dir.glob("*.data"))
        if not data_files:
            print(f"Warning: No data files found for N={N}!")
            continue
            
        mol_ranges[N] = {'start': current_mol_id, 'count': len(data_files)}
        
        for data_file in data_files:
            mol_id = current_mol_id
            mol_filename = f"mol_{mol_id}.mol"
            output_mol = mol_dir / mol_filename
            bbox = convert_data_to_molecule(str(data_file), str(output_mol))
            mol_bboxes[mol_id] = bbox
            mol_rel_path = str(output_mol).replace("\\", "/")
            combined_inc_lines.append(f"molecule m{mol_id} {mol_rel_path}")
            current_mol_id += 1
            
    mol_inc_path = mol_dir / "molecules.inc"
    with open(mol_inc_path, 'w') as f:
        f.write("\n".join(combined_inc_lines))
        
    print(f"--- Generating insertions.inc for {n_hoppers} hoppers ---")
    seed = 12345
    
    insertion_path_str, z_max = grid_manager._generate_grid_insertion_file(
        temp_dir, n_hoppers, n_fill_list, mol_ranges, seed, N_list, spacing, metadata, mol_bboxes=mol_bboxes, mode="2D_stacked"
    )
    print(f"Generated {insertion_path_str}, z_max = {z_max}")

    # 4. Create Dummy LAMMPS Script to generate a data file for visualization
    dummy_in = temp_dir / "dummy_visualize_mixed.in"
    
    # Calculate box to fit all hoppers
    x_max_h = (n_hoppers - 1) % nx_h * spacing + 1.0
    y_max_h = (n_hoppers - 1) // nx_h * spacing + 1.0

    with open(dummy_in, 'w') as f:
        f.write(f"""
units si
atom_style hybrid sphere molecular
boundary p p p

comm_modify vel yes

# Box expanded to fit all hoppers
region world block -1.0 {x_max_h} -1.0 {y_max_h} 0.0 {z_max + 0.5}
create_box 1 world bond/types 1 angle/types 1 extra/bond/per/atom 5 extra/angle/per/atom 5 extra/special/per/atom 10

# Include actual Lepton potentials
include {ROOT_DIR}/simulation_templates/lepton.inc

# Load molecule templates
include {mol_inc_path}

# Perform insertions
include {insertion_path_str}

mass 1 2.96e-05
set type 1 diameter 0.002

neighbor 0.005 bin
neigh_modify delay 0 every 1 check yes

thermo 1000

fix 1 all nve/sphere
run 5000

# Write out for visualization
write_data {temp_dir / 'inserted_state_mixed.data'}
""")

    print(f"--- Running dummy LAMMPS simulation for mixed hoppers ---")
    cmd = f"lmp -in {dummy_in}"
    import subprocess
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"SUCCESS: Data file created at {temp_dir / 'inserted_state_mixed.data'}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR running LAMMPS:\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")

    # 5. Cleanup
    print(f"--- Cleaning up temporary molecule files ---")
    import shutil
    shutil.rmtree(mol_dir)

if __name__ == "__main__":
    main()
