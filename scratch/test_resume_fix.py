import os
import sys
import json
from pathlib import Path

# Add root to path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from simulation.grid_hopper_manager import GridHopperManager

def test_resume_geometry_regeneration():
    parent_path = "dumping_yard/Hopper_Fill_Resume/Grid_Fill_1H_S469898_S906801_S672044_S388652_S172787_S623109"
    seed = 999999
    
    print(f"--- Testing Resume from {parent_path} ---")
    
    # Mock Runner
    class MockRunner:
        def run(self, config, clean_dir=True):
            print(f"MOCK RUN: Config generated at {config.output_dir}")
            print(f"Geometry Inc used: {config.extra_vars['geometry_inc']}")
            
            # Check if the file exists in the NEW directory
            local_inc = Path(config.extra_vars['geometry_inc'])
            if local_inc.exists():
                print(f"SUCCESS: Local geometry file found at {local_inc}")
            else:
                print(f"FAILURE: Local geometry file NOT found at {local_inc}")
    
    manager = GridHopperManager(runner=MockRunner())
    
    try:
        manager.resume_grid_filling(
            restart_path=parent_path,
            seed=seed,
            relax_steps=100,
            num_procs=1
        )
        
    except Exception as e:
        print(f"Error during test: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_resume_geometry_regeneration()
