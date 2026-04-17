import os
import time
import subprocess
import csv
from typing import List, Dict, Any
from simulation.orchestrator import SimulationOrchestrator

def parse_lammps_log(log_path: str) -> Dict[str, Any]:
    """Parses loop time and performance from LAMMPS log file."""
    results = {"loop_time": None, "perf_tau_sec": None, "perf_day": None}
    if not os.path.exists(log_path):
        return results

    with open(log_path, 'r') as f:
        lines = f.readlines()
        for line in reversed(lines):
            # Example: Loop time of 10.123 on 128 procs for 5000 steps with 4800 atoms
            if "Loop time of" in line:
                parts = line.split()
                results["loop_time"] = float(parts[3])
            # Example: Performance: 1234.567 tau/day, 0.012 day/tau
            if "Performance:" in line:
                # This format varies, let's try a safer parse
                parts = line.split()
                for i, p in enumerate(parts):
                    if "tau/day" in p:
                        results["perf_tau_day"] = float(parts[i-1])
                    if "ns/day" in p:
                        results["perf_ns_day"] = float(parts[i-1])
    return results

def run_benchmark():
    orch = SimulationOrchestrator()
    
    # Test Matrix
    # (label, use_kokkos, use_intel, num_threads)
    scenarios = [
        ("CPU-Standard-24T", False, False, 24),
        ("CPU-Standard-48T", False, False, 48),
        ("CPU-Standard-64T", False, False, 64),
        ("CPU-Standard-128T", False, False, 128),
        ("CPU-KOKKOS-64T", True, False, 64),
        ("CPU-INTEL-64T", False, True, 64),
        ("GPU-KOKKOS-1T", True, False, 1),
        ("GPU-KOKKOS-32T", True, False, 32),
        ("GPU-KOKKOS-64T", True, False, 64),
        ("GPU-KOKKOS-128T", True, False, 128),
    ]

    results_file = "benchmark_results.csv"
    fieldnames = ["scenario", "use_kokkos", "use_intel", "threads", "loop_time", "perf_tau_day", "perf_ns_day"]

    with open(results_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for label, kokkos, intel, threads in scenarios:
            print(f"\n>>> Running Benchmark: {label} ...")
            
            # Use run_hopper_fill as the benchmark task
            # 100 chains of N=48 poured into hopper
            # 5000 relax steps to measure performance
            outdir = f"benchmarking/{label}"
            
            # Clean old log if exists
            log_path = f"{outdir}/lammps.log"
            if os.path.exists(log_path):
                os.remove(log_path)

            start_t = time.time()
            try:
                orch.run_hopper_fill(
                    N=48,
                    n_fill=100,
                    relax_steps=5000,
                    num_threads=threads,
                    use_kokkos=kokkos,
                    use_intel=intel,
                    outdir=outdir
                )
            except Exception as e:
                print(f"Error in {label}: {e}")
                continue
            
            end_t = time.time()
            wall_t = end_t - start_t
            
            # Parse results
            metrics = parse_lammps_log(log_path)
            
            row = {
                "scenario": label,
                "use_kokkos": kokkos,
                "use_intel": intel,
                "threads": threads,
                "loop_time": metrics.get("loop_time"),
                "perf_tau_day": metrics.get("perf_tau_day"),
                "perf_ns_day": metrics.get("perf_ns_day")
            }
            writer.writerow(row)
            csvfile.flush()
            
            print(f"Finished {label} in {wall_t:.2f}s. Loop time: {metrics.get('loop_time')}s")

    print(f"\nBenchmarking complete. Results saved to {results_file}")

if __name__ == "__main__":
    run_benchmark()
