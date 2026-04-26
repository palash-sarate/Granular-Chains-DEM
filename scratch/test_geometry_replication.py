import sys
import os
from pathlib import Path

# Add project root to path
sys.path.append(os.getcwd())

from simulation.grid_hopper_manager import GridHopperManager

class MockRunner:
    pass

def test_replication():
    # Setup paths
    temp_dir = Path("temp")
    temp_dir.mkdir(exist_ok=True)
    
    input_inc = Path("simulation_geometries/2D_hopper.inc")
    n_hoppers = 4
    spacing = 2.0
    
    manager = GridHopperManager(MockRunner())
    
    print(f"Testing replication of {input_inc} for {n_hoppers} hoppers...")
    out_path_str, metadata = manager._generate_replicated_geometry(input_inc, n_hoppers, spacing, temp_dir)
    
    print(f"Generated file: {out_path_str}")
    print(f"Metadata: {metadata}")
    
    # Check a few lines of the output
    with open(out_path_str, 'r') as f:
        lines = f.readlines()
        print("\n--- First 20 lines of output ---")
        for line in lines[:20]:
            print(line.strip())

if __name__ == "__main__":
    test_replication()
