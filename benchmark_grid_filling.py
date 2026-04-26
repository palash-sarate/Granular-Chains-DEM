import os
import sys
import time
import glob
import pandas as pd
from simulation.orchestrator import SimulationOrchestrator

# Add project root to sys.path to support imports
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

def parse_lammps_log(log_path):
    """
    Extracts step and CPU time from the LAMMPS log file.
    Expects thermo_style custom step cpu ...
    """
    results = []
    if not os.path.exists(log_path):
        print(f"Warning: Log file not found at {log_path}")
        return results
    
    with open(log_path, 'r') as f:
        capture = False
        for line in f:
            # Look for the thermo header
            if "Step CPU PotEng KinEng TotEng" in line or "Step CPU" in line:
                capture = True
                continue
            if capture:
                if "Loop time" in line:
                    capture = False
                    continue
                parts = line.split()
                # Ensure it's a data line (starts with step number)
                if len(parts) >= 2 and parts[0].isdigit():
                    try:
                        results.append({
                            "step": int(parts[0]),
                            "cpu_time": float(parts[1])
                        })
                    except ValueError:
                        continue
    return results

def run_benchmark():
    orchestrator = SimulationOrchestrator(lammps_executable="lmp")
    
    # Configuration constants
    RELAX_STEPS = 10000
    N_FILL_PER_HOPPER = 2000  # Load per hopper to test scaling
    SOURCE_DIR = "chain_data/relaxed_2D_x"
    HOPPER_TEMPLATE = "simulation_geometries/2D_hopper.data"
    BENCHMARK_DUMP = "simulation_templates/benchmark_filling_dump.inc"
    SIM_NAME = "Grid_Filling_Benchmarks"
    
    all_metrics = []

    # --- Test 1: Hopper Scaling (Fixed num_procs = 8) ---
    # Vary n_hoppers to see how filling the GPU affects performance
    print("\n" + "="*60)
    print("SCALING TEST 1: Hopper Scaling (num_procs=8)")
    print("="*60)
    
    for n_h in [1, 5, 10]:
        print(f"\n>>> Running Benchmark: n_hoppers={n_h}, num_procs=8")
        orchestrator.run_grid_hopper_filling(
            n_hoppers=n_h,
            n_fill=N_FILL_PER_HOPPER,
            N=4,
            relax_steps=RELAX_STEPS,
            num_procs=8,
            dump_file=BENCHMARK_DUMP,
            simulation=SIM_NAME,
            source_dir=SOURCE_DIR,
            hopper_template_data=HOPPER_TEMPLATE
        )
        
        # Locate the log file
        log_pattern = f"dumping_yard/{SIM_NAME}/Grid_Fill_{n_h}H_N4_S*/lammps.log"
        logs = glob.glob(log_pattern)
        if logs:
            latest_log = max(logs, key=os.path.getmtime)
            timings = parse_lammps_log(latest_log)
            
            # Calculate intervals
            prev_cpu = 0
            for t in timings:
                interval_time = t["cpu_time"] - prev_cpu
                t.update({
                    "n_hoppers": n_h,
                    "num_procs": 8,
                    "test_type": "hopper_scaling",
                    "interval_time": interval_time
                })
                prev_cpu = t["cpu_time"]
                all_metrics.append(t)
            print(f"Captured {len(timings)} timing points for {n_h} hoppers.")

    # --- Test 2: Process Scaling (Fixed n_hoppers = 5) ---
    # Vary num_procs to find the sweet spot for MPI vs GPU work
    print("\n" + "="*60)
    print("SCALING TEST 2: Process Scaling (n_hoppers=5)")
    print("="*60)
    
    for n_p in [1, 4, 8, 12, 16]:
        print(f"\n>>> Running Benchmark: n_hoppers=5, num_procs={n_p}")
        orchestrator.run_grid_hopper_filling(
            n_hoppers=5,
            n_fill=N_FILL_PER_HOPPER,
            N=4,
            relax_steps=RELAX_STEPS,
            num_procs=n_p,
            dump_file=BENCHMARK_DUMP,
            simulation=SIM_NAME,
            source_dir=SOURCE_DIR,
            hopper_template_data=HOPPER_TEMPLATE
        )
        
        log_pattern = f"dumping_yard/{SIM_NAME}/Grid_Fill_5H_N4_S*/lammps.log"
        logs = glob.glob(log_pattern)
        if logs:
            latest_log = max(logs, key=os.path.getmtime)
            timings = parse_lammps_log(latest_log)
            
            prev_cpu = 0
            for t in timings:
                interval_time = t["cpu_time"] - prev_cpu
                t.update({
                    "n_hoppers": 5,
                    "num_procs": n_p,
                    "test_type": "process_scaling",
                    "interval_time": interval_time
                })
                prev_cpu = t["cpu_time"]
                all_metrics.append(t)
            print(f"Captured {len(timings)} timing points for {n_p} processes.")

    # Save to CSV
    if all_metrics:
        df = pd.DataFrame(all_metrics)
        output_csv = "grid_filling_benchmark_results.csv"
        df.to_csv(output_csv, index=False)
        print(f"\n{'='*60}")
        print(f"BENCHMARKING COMPLETE")
        print(f"Results saved to: {output_csv}")
        print(f"{'='*60}")
    else:
        print("\nError: No timing data captured. Check LAMMPS log files.")

if __name__ == "__main__":
    run_benchmark()
