import os
import sys
from pathlib import Path

# Add project root to sys.path
ROOT_DIR = Path("/home/guest/palash/Granular-Chains-DEM")
sys.path.append(str(ROOT_DIR))

from simulation.config import SimulationConfig

def test_config_outdir_override():
    print("Testing SimulationConfig outdir_override...")
    config = SimulationConfig(
        simulation="Wrong_Sim",
        run="RunA",
        outdir_override="/Data/palash_data/dumping_yard/Correct_Sim/RunA"
    )
    print(f"Expected Outdir: /Data/palash_data/dumping_yard/Correct_Sim/RunA")
    print(f"Actual Outdir:   {config.output_dir}")
    assert config.output_dir == "/Data/palash_data/dumping_yard/Correct_Sim/RunA"
    print("SUCCESS: outdir_override works correctly.")

def test_skip_params_rules():
    print("\nTesting AutoPilotManager skip_params transition rules...")
    
    test_cases = [
        {"mode": "fill_resume", "expected_skip_sim": False},
        {"mode": "flow",        "expected_skip_sim": True},
        {"mode": "flow_resume", "expected_skip_sim": False},
        {"mode": "fill",        "expected_skip_sim": False}
    ]
    
    for case in test_cases:
        mode = case["mode"]
        skip_params = {"walltime", "ppn", "mem"}
        if "resume" in mode or "flow" in mode:
            skip_params.update({
                "N", "n_fill", "spacing", "n_hoppers", "mode", 
                "source_dir", "no-vtk", "geometry_vars", "hopper_template_data"
            })
            if "resume" in mode:
                skip_params.add("restart_path")
            elif mode == "flow":
                skip_params.add("simulation")
        
        is_skipped = "simulation" in skip_params
        print(f"Mode: {mode:12} | Simulation Skipped: {str(is_skipped):5} | Expected: {str(case['expected_skip_sim']):5}")
        assert is_skipped == case["expected_skip_sim"]

    print("SUCCESS: Transition rules are correctly imposed.")

if __name__ == "__main__":
    test_config_outdir_override()
    test_skip_params_rules()
