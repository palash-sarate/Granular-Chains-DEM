import os
import time
import subprocess
import csv
from multiprocessing import Process
from typing import List, Dict, Any
from simulation.orchestrator import SimulationOrchestrator

def parse_lammps_log(log_path: str) -> Dict[str, Any]:
    """Parses loop time and performance from LAMMPS log file."""
    results = {"loop_time": None, "perf_tau_day": None, "perf_ns_day": None}
    if not os.path.exists(log_path):
        return results

    try:
        with open(log_path, 'r') as f:
            lines = f.readlines()
            for line in reversed(lines):
                if "Loop time of" in line:
                    parts = line.split()
                    results["loop_time"] = float(parts[3])
                if "Performance:" in line:
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if "tau/day" in p:
                            results["perf_tau_day"] = float(parts[i-1])
                        if "ns/day" in p:
                            results["perf_ns_day"] = float(parts[i-1])
    except:
        pass
    return results

def run_single_instance(run_id, N, n_fill, relax_steps, threads, use_kokkos, use_intel, label):
    """Worker function for parallel benchmarking."""
    orch = SimulationOrchestrator()
    outdir = f"benchmarking/parallel/{label}/run_{run_id}"
    
    try:
        orch.run_hopper_fill(
            N=N,
            n_fill=n_fill,
            relax_steps=relax_steps,
            num_threads=threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel,
            run_name=f"bench_{label}_{run_id}",
            outdir=outdir
        )
    except Exception as e:
        print(f"Error in {label} run {run_id}: {e}")

def run_benchmarks():
    orch = SimulationOrchestrator()
    results_file = "benchmark_results.csv"
    
    # --------------------------------------------------------------------------
    # PART 1: Sequential Scaling (1 Job at a time, different thread counts)
    # --------------------------------------------------------------------------
    seq_scenarios = [
        # (label, use_kokkos, use_intel, threads)
        ("CPU-Pure-64T", False, False, 64),
        ("KOKKOS-GPU-1T", True, False, 1),
        ("KOKKOS-GPU-16T", True, False, 16),
        ("KOKKOS-GPU-32T", True, False, 32),
        ("KOKKOS-GPU-64T", True, False, 64),
    ]

    fieldnames = ["type", "parallel_jobs", "threads_per_job", "use_kokkos", "use_intel", "total_wall_time", "throughput_states_per_hr", "loop_time"]
    
    with open(results_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        print("\n=== PART 1: SEQUENTIAL SCALING (N=4) ===")
        for label, kokkos, intel, threads in seq_scenarios:
            print(f">>> Testing Sequential: {label}...")
            outdir = f"benchmarking/sequential/{label}"
            log_path = f"{outdir}/lammps.log"
            if os.path.exists(log_path): os.remove(log_path)
            
            start_t = time.time()
            orch.run_hopper_fill(N=4, n_fill=200, relax_steps=10000, num_threads=threads, use_kokkos=kokkos, use_intel=intel, outdir=outdir)
            end_t = time.time()
            wall_t = end_t - start_t
            
            metrics = parse_lammps_log(log_path)
            writer.writerow({
                "type": "Sequential", "parallel_jobs": 1, "threads_per_job": threads, 
                "use_kokkos": kokkos, "use_intel": intel, "total_wall_time": round(wall_t, 2),
                "throughput_states_per_hr": round(3600 / wall_t, 2),
                "loop_time": metrics.get("loop_time")
            })
            csvfile.flush()

        # --------------------------------------------------------------------------
        # PART 2: Parallel Throughput (M Jobs at once, T threads each)
        # --------------------------------------------------------------------------
        print("\n=== PART 2: PARALLEL THROUGHPUT (N=4, Total Cores=64) ===")
        par_scenarios = [
            # (n_para, n_threads, use_kokkos, use_intel)
            (2, 32, True, False),  # 2 jobs on GPU, 32T each
            (2, 32, False, False), # 2 jobs on Pure CPU
            (4, 16, True, False),  # 4 jobs on GPU
            (4, 16, False, False), # 4 jobs on Pure CPU
            (8, 8, True, False),   # 8 jobs on GPU
            (8, 8, False, False),  # 8 jobs on Pure CPU
        ]

        for n_para, n_threads, kokkos, intel in par_scenarios:
            mode = "KOKKOS" if kokkos else "CPU"
            label = f"{n_para}x{n_threads}_{mode}"
            print(f">>> Testing Parallel: {label}...")
            
            start_t = time.time()
            processes = []
            for i in range(n_para):
                p = Process(target=run_single_instance, args=(i, 4, 200, 10000, n_threads, kokkos, intel, label))
                processes.append(p)
                p.start()
            
            for p in processes:
                p.join()
                
            end_t = time.time()
            wall_t = end_t - start_t
            
            throughput = (n_para / wall_t) * 3600
            writer.writerow({
                "type": "Parallel", "parallel_jobs": n_para, "threads_per_job": n_threads, 
                "use_kokkos": kokkos, "total_wall_time": round(wall_t, 2),
                "throughput_states_per_hr": round(throughput, 2),
                "loop_time": None
            })
            csvfile.flush()
            print(f"Finished {label} in {wall_t:.2f}s. Total Throughput: {throughput:.2f} states/hr")

    print(f"\n✅ All Benchmarks Complete. Results: {results_file}")

if __name__ == "__main__":
    run_benchmarks()
