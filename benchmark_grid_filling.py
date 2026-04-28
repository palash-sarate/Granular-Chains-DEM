import os
import sys
import time
import glob
import pandas as pd
import threading
import psutil
import subprocess
import re
from simulation.orchestrator import SimulationOrchestrator

# Add project root to sys.path to support imports
_script_dir = os.path.dirname(os.path.abspath(__file__))
_root_dir = os.path.dirname(_script_dir)
if _root_dir not in sys.path:
    sys.path.insert(0, _root_dir)

class ResourceMonitor(threading.Thread):
    def __init__(self, interval=0.5):
        super().__init__()
        self.interval = interval
        self.stopped = threading.Event()
        self.cpu_samples = []
        self.mem_samples = [] # MB
        self.gpu_util_samples = [] # List of max utils per poll
        self.gpu_mem_samples = [] # List of max mems per poll
        
    def get_gpu_metrics(self):
        try:
            # query gpu_util and mem_used for ALL gpus
            res = subprocess.check_output(
                ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader,nounits"],
                encoding='utf-8'
            )
            utils = []
            mems = []
            for line in res.strip().split('\n'):
                u, m = line.split(',')
                utils.append(float(u))
                mems.append(float(m))
            # Return MAX across all GPUs to capture the one being used by LAMMPS
            return max(utils), max(mems)
        except Exception:
            return 0.0, 0.0

    def run(self):
        while not self.stopped.is_set():
            self.cpu_samples.append(psutil.cpu_percent(interval=None))
            self.mem_samples.append(psutil.virtual_memory().used / (1024*1024))
            g_util, g_mem = self.get_gpu_metrics()
            self.gpu_util_samples.append(g_util)
            self.gpu_mem_samples.append(g_mem)
            time.sleep(self.interval)

    def stop(self):
        self.stopped.set()
        
    def get_stats(self):
        stats = {}
        if self.cpu_samples:
            stats['avg_cpu_util'] = sum(self.cpu_samples) / len(self.cpu_samples)
            stats['max_cpu_util'] = max(self.cpu_samples)
        if self.mem_samples:
            stats['max_ram_usage_mb'] = max(self.mem_samples)
        if self.gpu_util_samples:
            stats['avg_gpu_util'] = sum(self.gpu_util_samples) / len(self.gpu_util_samples)
            stats['max_gpu_util'] = max(self.gpu_util_samples)
        if self.gpu_mem_samples:
            stats['max_vram_usage_mb'] = max(self.gpu_mem_samples)
        return stats

def parse_lammps_log(log_path):
    results = []
    mem_per_proc = 0.0
    if not os.path.exists(log_path):
        return results, mem_per_proc
    
    with open(log_path, 'r') as f:
        capture = False
        for line in f:
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

def run_benchmark():
    orchestrator = SimulationOrchestrator(lammps_executable="lmp")
    
    RELAX_STEPS = 10000
    HOPPER_TEMPLATE = "simulation_geometries/2D_hopper.inc"
    BENCHMARK_DUMP = "simulation_templates/benchmark_filling_dump.inc"
    SIM_NAME = "Grid_Filling_Benchmarks"
    SIM_TEMPLATE = "in.grid_hopper_fill"
    
    all_metrics = []

    # --- Benchmark Combinations ---
    combos = [
        {"n_h": 6, "n_p": 16, "n_t": 1, "n_fill": 6000, "N": 4},
        {"n_h": 7, "n_p": 16, "n_t": 1, "n_fill": 6000, "N": 4},
        {"n_h": 8, "n_p": 16, "n_t": 1, "n_fill": 6000, "N": 4},
        {"n_h": 9, "n_p": 16, "n_t": 1, "n_fill": 6000, "N": 4},
    ]

    print("\n" + "="*60)
    print(f"RUNNING {len(combos)} BENCHMARK COMBINATIONS")
    print("="*60)

    for combo in combos:
        n_h = combo["n_h"]
        n_p = combo["n_p"]
        n_t = combo.get("n_t", 1)
        n_fill = combo["n_fill"]
        N = combo["N"]
        
        # Determine source directory for this N
        current_source_dir = f"chain_data/relaxed_2D_x/N{N}"
        
        print(f"\n>>> Running Combo: n_hoppers={n_h}, num_procs={n_p}, num_threads={n_t}, N={N}, n_fill={n_fill}")
        monitor = ResourceMonitor()
        monitor.start()
        
        python_setup_time = orchestrator.run_grid_hopper_filling(
            n_hoppers=n_h, n_fill=n_fill, N=N,
            relax_steps=RELAX_STEPS, num_procs=n_p, num_threads=n_t,
            dump_file=BENCHMARK_DUMP, simulation=SIM_NAME,
            source_dir=current_source_dir, hopper_template_data=HOPPER_TEMPLATE,
            template=SIM_TEMPLATE
        )
        
        monitor.stop()
        monitor.join()
        resource_stats = monitor.get_stats()
        
        # Note: We match the directory naming convention used by GridHopperManager
        log_pattern = f"dumping_yard/{SIM_NAME}/Grid_Fill_{n_h}H_N{N}_P{n_p}T{n_t}_S*/lammps.log"
        logs = glob.glob(log_pattern)
        if logs:
            latest_log = max(logs, key=os.path.getmtime)
            timings, lammps_mem = parse_lammps_log(latest_log)
            
            prev_cpu = 0
            for t in timings:
                interval_time = t["cpu_time"] - prev_cpu
                t.update({
                    "n_hoppers": n_h, "num_procs": n_p, "num_threads": n_t,
                    "n_fill": n_fill,
                    "chain_length": N,
                    "test_type": "combo_benchmark",
                    "source_dir": current_source_dir,
                    "hopper_template": HOPPER_TEMPLATE,
                    "sim_name": SIM_NAME,
                    "sim_template": SIM_TEMPLATE,
                    "interval_time": interval_time,
                    "python_setup_time": python_setup_time,
                    "is_setup_step": (t["step"] == 0),
                    "lammps_mem_per_proc_mb": lammps_mem
                })
                t.update(resource_stats)
                prev_cpu = t["cpu_time"]
                all_metrics.append(t)
            print(f"Captured {len(timings)} timing points (Setup: {python_setup_time:.2f}s) for {n_h} hoppers.")

    if all_metrics:
        df = pd.DataFrame(all_metrics)
        output_csv = "grid_filling_benchmark_results.csv"
        file_exists = os.path.isfile(output_csv)
        df.to_csv(output_csv, mode='a', index=False, header=not file_exists)
        print(f"\n{'='*60}\nBENCHMARKING COMPLETE\nResults saved to: {output_csv}\n{'='*60}")

if __name__ == "__main__":
    run_benchmark()
