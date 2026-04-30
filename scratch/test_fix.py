import os
import shutil
from pathlib import Path
from simulation.orchestrator import SimulationOrchestrator

def test_unique_mol_dirs():
    orchestrator = SimulationOrchestrator(lammps_executable="echo") # Use echo to avoid running real lammps
    
    # Test Grid Hopper Filling
    print("Testing Grid Hopper Filling...")
    try:
        # Mocking some parameters
        orchestrator.run_grid_hopper_filling(
            n_hoppers=1,
            n_fill=1,
            N=4,
            seed=12345,
            source_dir="chain_data/relaxed_2D_x",
            simulation="Test_Race_Condition"
        )
    except Exception as e:
        # It might fail because it tries to run LAMMPS or convert molecules
        # but we care about the directory creation before the run call
        print(f"Captured expected run-time exception or info: {e}")

    # Check dumping_yard
    base_dir = Path("dumping_yard/Test_Race_Condition")
    runs = list(base_dir.glob("Grid_Fill_*_S12345"))
    if runs:
        job_dir = runs[0]
        mol_dir = job_dir / "molecules"
        print(f"Checking {mol_dir}...")
        if mol_dir.exists():
            print("SUCCESS: Unique molecule directory created in job folder.")
            # Check molecules.inc content if it exists
            inc_file = mol_dir / "molecules.inc"
            if inc_file.exists():
                with open(inc_file, 'r') as f:
                    content = f.read()
                    print(f"molecules.inc content:\n{content}")
                    if str(mol_dir) in content:
                        print("SUCCESS: molecules.inc uses unique paths.")
                    else:
                        print("FAILURE: molecules.inc does not use unique paths.")
        else:
            print("FAILURE: Unique molecule directory NOT created.")
    else:
        print("FAILURE: Job directory NOT created.")

    # Cleanup
    if base_dir.exists():
        shutil.rmtree(base_dir)

if __name__ == "__main__":
    test_unique_mol_dirs()
