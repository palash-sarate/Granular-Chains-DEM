import os
import glob
import pandas as pd
import re
import numpy as np

def parse_lammps_log(log_path):
    results = []
    mem_per_proc = 0.0
    if not os.path.exists(log_path):
        return results, mem_per_proc
    
    with open(log_path, 'r') as f:
        capture = False
        for line in f:
            # Use the robust whitespace-insensitive check
            if "Step" in line and "CPU" in line:
                capture = True
                continue
            if capture:
                if "Loop time" in line:
                    capture = False
                    continue
                parts = line.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    try:
                        results.append({
                            "step": int(parts[0]),
                            "cpu_time": float(parts[1])
                        })
                    except ValueError:
                        continue
            
            if "Memory usage per processor =" in line:
                match = re.search(r"Memory usage per processor =\s+([\d\.]+)\s+Mbytes", line)
                if match:
                    mem_per_proc = float(match.group(1))
                    
    return results, mem_per_proc

def analyze():
    base_dir = "dumping_yard/Grid_Filling_Benchmarks"
    # Find all Grid_Fill directories
    dirs = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and "Grid_Fill" in d]
    
    # Sort by creation time (mtime is usually a good proxy on these systems)
    dirs.sort(key=os.path.getmtime)
    
    print(f"Found {len(dirs)} benchmark directories.")
    
    # Sequence of benchmarks from the script:
    # 1. Hopper Scaling (n_p=8): n_h = 1, 5, 10
    # 2. Process Scaling (n_h=5): n_p = 1, 4, 8, 12, 16
    
    test_sequence = [
        {"n_h": 1,  "n_p": 8,  "test": "hopper_scaling"},
        {"n_h": 5,  "n_p": 8,  "test": "hopper_scaling"},
        {"n_h": 10, "n_p": 8,  "test": "hopper_scaling"},
        {"n_h": 5,  "n_p": 1,  "test": "process_scaling"},
        {"n_h": 5,  "n_p": 4,  "test": "process_scaling"},
        {"n_h": 5,  "n_p": 8,  "test": "process_scaling"},
        {"n_h": 5,  "n_p": 12, "test": "process_scaling"},
        {"n_h": 5,  "n_p": 16, "test": "process_scaling"},
    ]
    
    all_metrics = []
    
    for i, d in enumerate(dirs):
        if i >= len(test_sequence):
            print(f"Warning: More directories found than expected in sequence. Skipping {d}")
            break
            
        config = test_sequence[i]
        log_path = os.path.join(d, "lammps.log")
        
        print(f"Analyzing {d} as {config['test']} (n_h={config['n_h']}, n_p={config['n_p']})...")
        
        timings, lammps_mem = parse_lammps_log(log_path)
        
        if not timings:
            print(f"  Warning: No timing data found in {log_path}")
            continue
            
        prev_cpu = 0
        for t in timings:
            interval_time = t["cpu_time"] - prev_cpu
            t.update({
                "n_hoppers": config["n_h"],
                "num_procs": config["n_p"],
                "n_fill": 2000,
                "test_type": config["test"],
                "interval_time": interval_time,
                "is_setup_step": (t["step"] == 0),
                "lammps_mem_per_proc_mb": lammps_mem,
                "dir_name": os.path.basename(d)
            })
            prev_cpu = t["cpu_time"]
            all_metrics.append(t)
            
    if all_metrics:
        df = pd.DataFrame(all_metrics)
        output_csv = "grid_filling_benchmark_results_restored.csv"
        df.to_csv(output_csv, index=False)
        print(f"\nAnalysis complete. Results saved to: {output_csv}")
    else:
        print("\nError: No metrics could be extracted.")

if __name__ == "__main__":
    analyze()
