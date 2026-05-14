import os
import subprocess
import re
import json
import time
from datetime import datetime
from typing import List, Dict, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import psutil
import sys

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class PBSManager:
    METADATA_FILE = "Pulse/pulse_metadata.json"

    @staticmethod
    def load_lineage() -> Dict:
        """Loads the simulation lineage data."""
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lineage.json")
        if os.path.exists(path):
            try:
                with open(path, 'r') as f:
                    return json.load(f)
            except Exception:
                return {}
        return {}

    @staticmethod
    def save_lineage(data: Dict):
        """Saves simulation lineage data."""
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lineage.json")
        with open(path, 'w') as f:
            json.dump(data, f, indent=4)

    @staticmethod
    def load_metadata() -> Dict:
        """Loads persistent job metadata."""
        if os.path.exists(PBSManager.METADATA_FILE):
            try:
                with open(PBSManager.METADATA_FILE, 'r') as f:
                    return json.load(f)
            except Exception:
                return {}
        return {}

    @staticmethod
    def save_metadata(data: Dict):
        """Saves job metadata to file."""
        os.makedirs(os.path.dirname(PBSManager.METADATA_FILE), exist_ok=True)
        with open(PBSManager.METADATA_FILE, 'w') as f:
            json.dump(data, f, indent=4)

    @staticmethod
    def format_mem(mem_str: str) -> str:
        """Converts a memory string (e.g. 1339240kb) to a human-readable format (MB/GB)."""
        if not mem_str or mem_str == "N/A":
            return mem_str
        
        # Extract numeric value
        match = re.search(r'(\d+)', mem_str.lower())
        if not match:
            return mem_str
        
        kb = int(match.group(1))
        if kb > 1024 * 1024:
            return f"{kb / (1024 * 1024):.2f} GB"
        elif kb > 1024:
            return f"{kb / 1024:.2f} MB"
        return f"{kb} KB"

    @staticmethod
    def get_jobs(user: Optional[str] = None) -> List[Dict]:
        """Gets a list of jobs, optionally filtered by user."""
        try:
            # Step 1: Get the list of Job IDs for the user
            # We use qstat -u as it's the most reliable way to get a user's job list across PBS versions
            if user:
                id_cmd = ["qstat", "-u", user]
            else:
                id_cmd = ["qstat"]
                
            id_res = subprocess.run(id_cmd, capture_output=True, text=True, timeout=10)
            if id_res.returncode != 0:
                return []
            
            # Extract IDs from the table (first column)
            job_ids = []
            lines = id_res.stdout.splitlines()
            for line in lines:
                parts = line.split()
                # Job IDs usually look like '1234.master' or '1234'
                if parts and ("." in parts[0] or parts[0].isdigit()):
                    job_ids.append(parts[0])
            
            if not job_ids:
                return []
            
            # Step 2: Fetch detailed info for these specific IDs
            # We limit to batches of 50 to avoid command line length limits
            all_jobs = []
            for i in range(0, len(job_ids), 50):
                batch = job_ids[i:i+50]
                detail_cmd = ["qstat", "-f"] + batch
                detail_res = subprocess.run(detail_cmd, capture_output=True, text=True, timeout=15)
                
                if detail_res.returncode == 0:
                    all_jobs.extend(PBSManager.parse_qstat_f(detail_res.stdout, user))
                
            return all_jobs
        except Exception as e:
            print(f"Error in get_jobs: {e}")
            return []

    @staticmethod
    def parse_qstat_f(output: str, filter_user: Optional[str] = None) -> List[Dict]:
        """Parses the output of qstat -f into a list of dictionaries."""
        jobs = []
        current_job = {}
        
        # Split by "Job Id:" to get individual job blocks
        blocks = re.split(r'Job Id:\s+', output)
        
        for block in blocks:
            if not block.strip():
                continue
                
            lines = block.splitlines()
            job_id = lines[0].strip()
            job_data = {"id": job_id}
            
            # Use a regex to catch "Attribute = Value"
            # Some values might span multiple lines if they are indented
            last_key = None
            for line in lines[1:]:
                match = re.match(r'^\s+([\w.]+)\s+=\s+(.*)$', line)
                if match:
                    key, value = match.groups()
                    job_data[key] = value
                    last_key = key
                elif last_key and (line.startswith("    ") or line.startswith("\t")):
                    # Continuation line (handles both spaces and tabs)
                    job_data[last_key] += line.strip()
            
            # Post-process memory fields
            for key in ["resources_used.mem", "resources_used.vmem"]:
                if key in job_data:
                    job_data[key] = PBSManager.format_mem(job_data[key])
            
            # Filter by user if requested
            if filter_user:
                owner = job_data.get("Job_Owner", "")
                if filter_user not in owner:
                    continue
                    
            jobs.append(job_data)
        return jobs

    @staticmethod
    def get_utilization(job_data: Dict) -> Dict:
        """Extracts utilization metrics from job data."""
        return {
            "cpu_percent": job_data.get("resources_used.cpupercent", "0"),
            "mem_used": job_data.get("resources_used.mem", "0kb"),
            "vmem_used": job_data.get("resources_used.vmem", "0kb"),
            "walltime": job_data.get("resources_used.walltime", "00:00:00")
        }

    @staticmethod
    def submit_job(script_path: str, extra_args: List[str] = None) -> str:
        """Submits a job and returns the Job ID."""
        cmd = ["qsub"]
        if extra_args:
            cmd.extend(extra_args)
        cmd.append(script_path)
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            raise Exception(f"Submission failed: {result.stderr}")

    @staticmethod
    def delete_job(job_id: str):
        """Deletes a job."""
        subprocess.run(["qdel", job_id], check=True)

    @staticmethod
    def trace_job(job_id: str, cache: Dict = None) -> Dict[str, str]:
        """Runs tracejob and parses the resource usage summary. Uses cache if provided."""
        # Check cache first
        if cache and job_id in cache:
            return cache[job_id]

        try:
            # Use -n 7 to look back a week
            result = subprocess.run(["tracejob", "-n", "7", job_id], capture_output=True, text=True, timeout=15)
            if result.returncode != 0:
                return {}
            
            output = result.stdout
            stats = {}
            
            # Extract resources_used line
            # Format: 04/24/2026 17:07:42  S    Exit_status=0 resources_used.cpupercent=1429 resources_used.cput=01:50:10 resources_used.mem=1339240kb ...
            all_matches = re.findall(r'resources_used\.(\w+)=([\w:]+)', output)
            for key, value in all_matches:
                if key in ["mem", "vmem"]:
                    stats[key] = PBSManager.format_mem(value)
                else:
                    stats[key] = value
                
            # Also catch Exit_status
            exit_match = re.search(r'Exit_status=(\d+)', output)
            if exit_match:
                stats["exit_status"] = exit_match.group(1)
            
            # Extract Owner and Job Name (often in the "Job Queued" line)
            owner_match = re.search(r'owner\s*=\s*([\w@.]+)', output)
            if owner_match:
                stats["owner"] = owner_match.group(1).split('@')[0] # Just the username
                
            name_match = re.search(r'job name\s*=\s*([\w_-]+)', output)
            if name_match:
                stats["job_name"] = name_match.group(1)
                
            return stats
        except Exception:
            return {}

    @staticmethod
    def get_job_history(user: str, days: int = 1) -> List[Dict]:
        """
        Attempts to find recently completed jobs for a user.
        Since qstat -x is often disabled, we search for common PBS output files
        to discover Job IDs, then trace them.
        """
        job_ids = set()
        
        # 1. Search for common PBS output files like *.o[0-9]* or *.e[0-9]*
        # We search in home and project directories
        search_paths = [f"/home/{user}", os.getcwd()]
        for path in search_paths:
            if not os.path.exists(path):
                continue
            try:
                # Find files matching the pattern (e.g., job.o1234)
                # We use a simple listdir + regex to avoid slow find
                for filename in os.listdir(path):
                    match = re.search(r'\.[oe](\d+)$', filename)
                    if match:
                        job_ids.add(match.group(1))
            except Exception:
                continue
        
        # 2. Trace found jobs
        history = []
        cache = PBSManager.load_metadata()
        
        # Merge new IDs with cache
        for jid in sorted(list(job_ids), reverse=True)[:10]: # Limit to last 10
            trace_data = PBSManager.trace_job(jid, cache=cache)
            if trace_data:
                # Update cache
                cache[jid] = trace_data
                if trace_data.get("owner") == user or trace_data.get("owner") == "N/A":
                    history.append({"id": jid, **trace_data})
        
        PBSManager.save_metadata(cache)
        return history

    @staticmethod
    def scan_range(start: int, end: int, user_filter: str, progress_callback=None) -> List[Dict]:
        """Scans a range of Job IDs and updates metadata."""
        cache = PBSManager.load_metadata()
        total = end - start + 1
        
        for i, jid_int in enumerate(range(start, end + 1)):
            jid = str(jid_int)
            if jid not in cache:
                trace_data = PBSManager.trace_job(jid)
                if trace_data:
                    cache[jid] = trace_data
                else:
                    # Mark as not found to avoid re-parsing
                    cache[jid] = {"owner": "N/A", "status": "not_found"}
            
            if progress_callback:
                progress_callback(i + 1, total)
                
            # Periodic save every 50 jobs
            if (i + 1) % 50 == 0:
                PBSManager.save_metadata(cache)
                
        PBSManager.save_metadata(cache)
        
        # Filter for the requested user
        results = []
        for jid, data in cache.items():
            if data.get("owner") == user_filter:
                results.append({"id": jid, **data})
        return sorted(results, key=lambda x: int(x["id"]) if x["id"].isdigit() else 0, reverse=True)

    @staticmethod
    def delete_cached_job(job_id: str):
        """Removes a job from the metadata cache."""
        cache = PBSManager.load_metadata()
        if job_id in cache:
            del cache[job_id]
            PBSManager.save_metadata(cache)

    @staticmethod
    def get_active_seeds(user: Optional[str] = None) -> List[str]:
        """Returns a list of seeds currently running in PBS jobs."""
        jobs = PBSManager.get_jobs(user=user)
        active_seeds = []
        for job in jobs:
            job_name = job.get('Job_Name', '')
            # Match convention: [Prefix][ParentSeed]_[ChildSeed]
            # e.g., S374952_481365 or R481365_956135
            match = re.match(r"^[A-Z]+(\d{3,6})_(\d{3,6})$", job_name)
            if match:
                parent_seed, child_seed = match.groups()
                active_seeds.append(child_seed)
        return active_seeds



    @staticmethod
    def get_node_temperatures() -> Dict[str, float]:
        """Fetches CPU and GPU temperatures from the system."""
        temps = {}
        
        # 1. GPU Temperature (via nvidia-smi)
        try:
            gpu_out = subprocess.getoutput("nvidia-smi --query-gpu=temperature.gpu --format=csv,noheader,nounits")
            if gpu_out and gpu_out.isdigit():
                temps["gpu"] = float(gpu_out)
        except Exception:
            pass
            
        # 2. CPU Temperature (via hwmon)
        try:
            # We look for Tdie or similar in hwmon
            # Find which hwmon is the CPU
            for i in range(10):
                path = f"/sys/class/hwmon/hwmon{i}/name"
                if os.path.exists(path):
                    with open(path, 'r') as f:
                        name = f.read().strip()
                    
                    # k10temp is common for AMD, coretemp for Intel
                    if name in ["k10temp", "coretemp"]:
                        # Try to find a better sensor than Tctl (which has offsets)
                        # We look for Tccd1 or Tdie
                        target_sensor = "temp1_input" # Default fallback
                        for j in range(1, 10):
                            label_path = f"/sys/class/hwmon/hwmon{i}/temp{j}_label"
                            if os.path.exists(label_path):
                                with open(label_path, 'r') as f:
                                    label = f.read().strip()
                                if label in ["Tccd1", "Tdie"]:
                                    target_sensor = f"temp{j}_input"
                                    break
                                elif label == "Tctl":
                                    # If Tctl is all we have, we'll use it but it's less ideal
                                    pass

                        temp_path = f"/sys/class/hwmon/hwmon{i}/{target_sensor}"
                        if os.path.exists(temp_path):
                            with open(temp_path, 'r') as f:
                                temps["cpu"] = float(f.read().strip()) / 1000.0
                        break
        except Exception:
            pass
            
        return temps

class SimulationMonitor:
    @staticmethod
    def estimate_eta(directory: str, target_relative_steps: int) -> Dict:
        """
        Estimates simulation ETA by treating 'target' as the number of steps 
        to be performed in the CURRENT run (relative steps).
        """
        if not os.path.exists(directory):
            return {"error": f"Directory not found: {directory}"}

        # 1. Collect all dump files
        search_dirs = [directory]
        chain_sub = os.path.join(directory, "chain")
        if os.path.isdir(chain_sub):
            search_dirs.insert(0, chain_sub)
            
        raw_files = []
        try:
            for s_dir in search_dirs:
                for f in os.listdir(s_dir):
                    match = re.match(r'chain_(\d+)\.dump', f)
                    if match:
                        timestep = int(match.group(1))
                        raw_files.append((timestep, os.path.getmtime(os.path.join(s_dir, f))))
                if raw_files: break
        except Exception as e:
            return {"error": f"Error scanning directory: {str(e)}"}

        if not raw_files:
            return {"error": "No dump files found matching 'chain_*.dump'"}

        raw_files.sort() # Sort by absolute timestep
        
        # 2. Convert to Relative Timesteps (0 is start of this job)
        initial_abs_ts = raw_files[0][0]
        files_rel = [(f[0] - initial_abs_ts, f[1]) for f in raw_files]
        
        current_rel_ts = files_rel[-1][0]
        current_abs_ts = raw_files[-1][0]
        current_time = time.time()
        
        # 3. Check Progress
        if current_rel_ts >= target_relative_steps:
            return {
                "status": "Target Reached",
                "current_timestep": current_abs_ts,
                "current_relative_step": current_rel_ts,
                "target_relative_steps": target_relative_steps,
                "progress_percent": 100.0,
                "time_remaining_hr": 0.0,
                "time_elapsed_hr": (files_rel[-1][1] - files_rel[0][1]) / 3600.0,
                "completion_time": datetime.fromtimestamp(files_rel[-1][1]).strftime("%Y-%m-%d %H:%M:%S"),
                "cost_per_10k_steps": 0.0,
                "total_files": len(files_rel),
                "data_points": len(files_rel),
                "last_updated": datetime.now().strftime("%H:%M:%S"),
                "message": f"Reached target relative duration: {target_relative_steps:,}"
            }

        # 4. Sampling & Quadratic Extrapolation
        total_files = len(files_rel)
        if total_files > 100:
            step = total_files // 50 # Aim for ~50-100 points
            sampled_rel = files_rel[::step]
        else:
            sampled_rel = files_rel
            
        if files_rel[-1] not in sampled_rel:
            sampled_rel.append(files_rel[-1])
            sampled_rel.sort()

        if len(sampled_rel) < 2:
            return {"error": "Insufficient data: Need at least 2 dump files."}

        x_arr = np.array([f[0] for f in sampled_rel])
        y_arr = np.array([f[1] for f in sampled_rel])

        degree = 2 if len(sampled_rel) >= 3 else 1
        coeffs = np.polyfit(x_arr, y_arr, degree)
        poly = np.poly1d(coeffs)
        
        estimated_completion_time = poly(target_relative_steps)
        
        # Linear safety fallback
        last_x, last_y = files_rel[-1]
        prev_x, prev_y = files_rel[-min(len(files_rel), 5)]
        if last_x != prev_x:
            latest_rate = (last_y - prev_y) / (last_x - prev_x)
            linear_eta = last_y + (target_relative_steps - last_x) * latest_rate
            if estimated_completion_time < current_time or estimated_completion_time < linear_eta:
                estimated_completion_time = linear_eta

        # 5. Result Generation
        time_remaining_hr = max(0, estimated_completion_time - current_time) / 3600.0
        time_elapsed_hr = (current_time - files_rel[0][1]) / 3600.0
        
        deriv = np.polyder(poly)
        current_cost_per_step = deriv(current_rel_ts)
        step_interval = files_rel[-1][0] - files_rel[-2][0] if len(files_rel) >= 2 else 1000

        return {
            "status": "Active",
            "current_timestep": current_abs_ts,
            "current_relative_step": current_rel_ts,
            "target_relative_steps": target_relative_steps,
            "step_interval": step_interval,
            "progress_percent": (current_rel_ts / target_relative_steps) * 100,
            "time_remaining_hr": time_remaining_hr,
            "time_elapsed_hr": time_elapsed_hr,
            "completion_time": datetime.fromtimestamp(estimated_completion_time).strftime("%Y-%m-%d %H:%M:%S"),
            "cost_per_10k_steps": current_cost_per_step * 10000 / 60.0,
            "cost_per_dump": (current_cost_per_step * step_interval) / 60.0,
            "cost_history": (deriv(x_arr) * 10000 / 60.0).tolist(),
            "total_files": total_files,
            "data_points": len(sampled_rel),
            "last_updated": datetime.now().strftime("%H:%M:%S")
        }

class SyncManager:
    SYNC_LOG = "Pulse/sync.log"
    SYNC_LOCK = "Pulse/sync.lock"
    SYNC_STOP = "Pulse/.sync_stop"
    # SYNC_SCRIPT = "Pulse/sync_drive.sh"

    @staticmethod
    def start_sync(user: str = "guest"):
        """Starts the sync process in a background thread."""
        if SyncManager.is_running():
            return False, "Sync is already running."

        # Clear any old stop requests
        stop_path = os.path.join(ROOT_DIR, SyncManager.SYNC_STOP)
        if os.path.exists(stop_path):
            os.remove(stop_path)

        import threading
        thread = threading.Thread(target=SyncManager.run_sync_cycle, args=(user,))
        thread.daemon = True
        
        # We still use a lock file to track it globally
        lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
        with open(lock_path, "w") as f:
            f.write(str(os.getpid())) # Note: in-thread sync uses current PID
            
        thread.start()
        return True, "Sync started in background."

    @staticmethod
    def start_sync_pbs(user: str = "guest"):
        """Submits the sync process as a PBS job."""
        if SyncManager.is_running():
            return False, "Sync is already running."

        # Clear any old stop requests
        stop_path = os.path.join(ROOT_DIR, SyncManager.SYNC_STOP)
        if os.path.exists(stop_path):
            os.remove(stop_path)

        job_script = f"""#!/bin/bash
#PBS -N Sync_Cycle
#PBS -q workq
#PBS -l nodes=master:ppn=8
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o {os.path.join(ROOT_DIR, SyncManager.SYNC_LOG)}

cd $PBS_O_WORKDIR

# Ensure all child processes are killed on exit
trap 'kill 0' EXIT

# Activate conda environment if available
if [ -f /home/guest/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/guest/miniconda3/etc/profile.d/conda.sh
    conda activate gchain
fi

python Pulse/run_sync_job.py --user {user}
"""
        
        script_path = os.path.join(ROOT_DIR, "Pulse/temp_sync.pbs")
        with open(script_path, "w") as f:
            f.write(job_script)
        
        try:
            result = subprocess.run(["qsub", script_path], capture_output=True, text=True, check=True)
            job_id = result.stdout.strip()
            
            # Write Job ID to lock file (prepended with PBS: to distinguish)
            lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
            with open(lock_path, "w") as f:
                f.write(f"PBS:{job_id}")
            
            return True, f"Sync submitted as PBS job: {job_id}"
        except Exception as e:
            return False, f"Failed to submit PBS job: {e}"
        finally:
            if os.path.exists(script_path):
                os.remove(script_path)

    @staticmethod
    def run_sync_cycle(user: str = "guest"):
        """The main sync logic: Zip, Move, Update Lineage, and Cleanup."""
        log_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOG)
        lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
        stop_path = os.path.join(ROOT_DIR, SyncManager.SYNC_STOP)
        
        try:
            with open(log_path, "w") as log:
                def log_msg(msg):
                    print(msg)
                    log.write(f"{datetime.now().strftime('%H:%M:%S')} - {msg}\n")
                    log.flush()

                def check_stop():
                    if os.path.exists(stop_path):
                        log_msg("STOP REQUESTED. Terminating sync cycle...")
                        return True
                    return False

                log_msg("Starting Sync Cycle...")
                
                # 0. Rescan Lineage from disk to get latest state
                try:
                    from Pulse.lineage_tracker import scan_dumping_yard
                    scan_dumping_yard()
                    log_msg("Lineage rescan complete.")
                except Exception as e:
                    log_msg(f"Warning: Lineage rescan failed: {e}")

                # 1. Load Lineage
                lineage = PBSManager.load_lineage()
                if not lineage:
                    log_msg("No lineage data found. Skipping.")
                    return

                # 2. Identify Ongoing Runs (to skip)
                active_seeds = PBSManager.get_active_seeds(user=user)
                log_msg(f"Active seeds in PBS: {active_seeds}")

                # 3. Identify Syncable Runs
                # Rules: Not currently running, has local data, and not already synced (or dirty)
                syncable = []
                for path, info in lineage.items():
                    seed = str(info.get('params', {}).get('seed', ''))
                    if seed in active_seeds:
                        continue
                    
                    if not os.path.exists(path):
                        continue
                    
                    # Check if dirty or never synced
                    status = info.get('sync_status', 'Local')
                    if status != 'Synced':
                        syncable.append(path)
                
                log_msg(f"Found {len(syncable)} runs to check for sync.")
                
                # 4. Cloud Inventory Check (Optimization)
                # Get list of existing zips on Drive to avoid redundant compression
                log_msg("Fetching cloud inventory from Google Drive...")
                cloud_zips = set()
                try:
                    res = subprocess.run([
                        "rclone", "lsf", 
                        "gdrive:Granular-Chains-DEM/dumping_yard_zips"
                    ], capture_output=True, text=True, timeout=120)
                    if res.returncode == 0:
                        cloud_zips = set(res.stdout.splitlines())
                        log_msg(f"Found {len(cloud_zips)} files already on Drive.")
                except Exception as e:
                    log_msg(f"Warning: Could not fetch cloud inventory: {e}")

                # 5. Parallel Compression
                ZIP_DIR = os.path.join(ROOT_DIR, "Zips_dir")
                os.makedirs(ZIP_DIR, exist_ok=True)
                
                processed_zips = []
                lineage_changed = False
                
                def compress_task(path):
                    info = lineage[path]
                    rel_path = info.get('rel_path', os.path.basename(path))
                    unique_name = rel_path.replace('/', '_').replace('\\', '_').replace(' ', '_')
                    zip_name = unique_name + ".tar.gz"
                    
                    if zip_name in cloud_zips and info.get('sync_status') == 'Synced':
                        return path, None, True
                        
                    zip_path = os.path.join(ZIP_DIR, zip_name)
                    try:
                        subprocess.run(["tar", "-czf", zip_path, "-C", os.path.dirname(path), os.path.basename(path)], check=True)
                        return path, zip_path, False
                    except Exception as e:
                        if os.path.exists(zip_path): os.remove(zip_path) # Cleanup partial file
                        raise e

                log_msg(f"Starting parallel compression with 8 workers...")
                with ThreadPoolExecutor(max_workers=8) as executor:
                    futures = {executor.submit(compress_task, path): path for path in syncable}
                    
                    completed = 0
                    for future in as_completed(futures):
                        if check_stop(): break
                        
                        try:
                            path, zip_path, already_synced = future.result()
                            info = lineage[path]
                            completed += 1
                            
                            if already_synced:
                                log_msg(f"Skipping {info['name']} (Already in Cloud Inventory)")
                                lineage[path]['sync_status'] = 'Synced'
                                lineage[path]['sync_time'] = lineage[path].get('sync_time') or datetime.now().isoformat()
                                lineage_changed = True
                            else:
                                log_msg(f"PROGRESS_COMPRESS: {completed} / {len(syncable)}")
                                log_msg(f"Finished compressing {info['name']}.")
                                if zip_path: processed_zips.append((path, zip_path))
                        except Exception as e:
                            log_msg(f"ERROR: Failed to compress {futures[future]}: {e}")
                            continue # Don't crash the whole cycle for one bad folder
                
                if lineage_changed:
                    PBSManager.save_lineage(lineage)

                # 6. Bulk Move
                if processed_zips:
                    log_msg("Starting bulk upload to Google Drive...")
                    cmd = [
                        "rclone", "move", ZIP_DIR, 
                        "gdrive:Granular-Chains-DEM/dumping_yard_zips",
                        "--transfers", "16",
                        "--drive-chunk-size", "128M",
                        "--buffer-size", "256M",
                        "--progress"
                    ]
                    
                    # We use Popen so we can pipe output to our log for the dashboard
                    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    for line in proc.stdout:
                        if os.path.exists(stop_path):
                            log_msg("STOP REQUESTED. Killing rclone...")
                            import signal
                            proc.terminate()
                            break
                        log.write(line)
                        log.flush()
                    proc.wait()
                    
                    if proc.returncode == 0:
                        log_msg("Upload successful.")
                        # Update Lineage Metadata
                        for path, z_path in processed_zips:
                            if not os.path.exists(z_path): # Verified moved
                                lineage[path]['sync_status'] = 'Synced'
                                lineage[path]['sync_time'] = datetime.now().isoformat()
                        PBSManager.save_lineage(lineage)
                    else:
                        log_msg(f"Upload failed with code {proc.returncode}")
                
                # 6. Smart Cleanup (Refined Logic)
                if not os.path.exists(stop_path):
                    log_msg("Starting Smart Cleanup...")
                
                # Pre-calculate children and ghost status
                child_map = {} # parent_path -> list of child_paths
                for path, info in lineage.items():
                    p = info.get('parent')
                    if p:
                        if p not in child_map: child_map[p] = []
                        child_map[p].append(path)
                
                # Find Ghost children in PBS (jobs that haven't written metadata yet)
                jobs = PBSManager.get_jobs(user=user)
                # Map all known seeds to their paths
                seed_to_path = {str(info.get('params', {}).get('seed', '')): path for path, info in lineage.items() if info.get('params', {}).get('seed')}
                
                ghost_parents = set()
                for job in jobs:
                    job_name = job.get('Job_Name', '')
                    # If any job in the queue contains a seed we know, protect that path
                    for seed, path in seed_to_path.items():
                        if seed in job_name:
                            ghost_parents.add(path)

                cleaned_count = 0
                for path, info in lineage.items():
                    # Condition 1: Must be Synced
                    if info.get('sync_status') != 'Synced':
                        continue
                    
                    children = child_map.get(path, [])
                    has_ghosts = path in ghost_parents
                    
                    # Condition 2: Must NOT be the last run (Leaf Node Protection)
                    if len(children) == 0 and not has_ghosts:
                        # This is a leaf node, we keep it on disk for convenience
                        continue
                    
                    # Condition 3: ALL children must be Synced
                    # If there are ghosts, it's definitely not safe to delete
                    if has_ghosts:
                        continue
                    
                    all_children_synced = True
                    for c_path in children:
                        if lineage.get(c_path, {}).get('sync_status') != 'Synced':
                            all_children_synced = False
                            break
                    
                    if all_children_synced:
                        if os.path.exists(path):
                            log_msg(f"Cleaning up local files for: {info['name']}")
                            import shutil
                            shutil.rmtree(path, ignore_errors=True)
                            cleaned_count += 1
                
                log_msg(f"Cleanup complete. Removed {cleaned_count} parent directories.")
                log_msg("Sync Cycle Finished Successfully.")

        except Exception as e:
            with open(log_path, "a") as log:
                log.write(f"CRITICAL ERROR: {str(e)}\n")
                import traceback
                log.write(traceback.format_exc())
        finally:
            if os.path.exists(lock_path):
                os.remove(lock_path)
            if os.path.exists(stop_path):
                os.remove(stop_path)

    @staticmethod
    def restore_run(run_path: str):
        """Downloads a run back from Drive and extracts it."""
        lineage = PBSManager.load_lineage()
        if run_path not in lineage:
            return False, "Run not found in lineage."
        
        info = lineage[run_path]
        
        # Use the same unique naming convention as the sync process
        rel_path = info.get('rel_path', os.path.basename(run_path))
        unique_name = rel_path.replace('/', '_').replace('\\', '_').replace(' ', '_')
        zip_name = unique_name + ".tar.gz"
        
        ZIP_DIR = os.path.join(ROOT_DIR, "Zips_dir")
        os.makedirs(ZIP_DIR, exist_ok=True)
        zip_path = os.path.join(ZIP_DIR, zip_name)
        
        try:
            # 1. Download
            subprocess.run([
                "rclone", "copy", 
                f"gdrive:Granular-Chains-DEM/dumping_yard_zips/{zip_name}", 
                ZIP_DIR
            ], check=True)
            
            # 2. Unzip
            # Ensure target directory exists
            target_parent = os.path.dirname(run_path)
            os.makedirs(target_parent, exist_ok=True)
            
            subprocess.run(["tar", "-xzf", zip_path, "-C", target_parent], check=True)
            
            # 3. Cleanup zip
            if os.path.exists(zip_path):
                os.remove(zip_path)
            
            # 4. Update lineage
            # We keep it as Synced because it is still on Drive, but now also local.
            lineage[run_path]['sync_status'] = 'Synced'
            PBSManager.save_lineage(lineage)
            
            return True, f"Successfully restored {info['name']} to {run_path}"
        except Exception as e:
            return False, f"Restore failed: {str(e)}"

    @staticmethod
    def free_restored_space(run_paths: List[str]):
        """
        Clears local folders for runs that are confirmed synced on Drive.
        This is used to clean up runs that were temporarily restored for tasks like visualization.
        """
        lineage = PBSManager.load_lineage()
        freed_count = 0
        for path in run_paths:
            # Safety check: Must be in lineage and marked as Synced
            if path in lineage and lineage[path].get('sync_status') == 'Synced':
                if os.path.exists(path):
                    import shutil
                    shutil.rmtree(path, ignore_errors=True)
                    freed_count += 1
        return freed_count

    @staticmethod
    def is_running():
        """Checks if a sync process is currently running (locally or on PBS)."""
        lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
        if not os.path.exists(lock_path):
            return False
        
        try:
            with open(lock_path, "r") as f:
                lock_content = f.read().strip()
            
            if lock_content.startswith("PBS:"):
                job_id = lock_content.replace("PBS:", "")
                # Check if PBS job is still in queue
                jobs = PBSManager.get_jobs()
                for job in jobs:
                    if job.get('id') == job_id:
                        return True
            else:
                pid = int(lock_content)
                # Check if local process exists
                import psutil
                if psutil.pid_exists(pid):
                    proc = psutil.Process(pid)
                    if proc.is_running() and proc.status() != psutil.STATUS_ZOMBIE:
                        return True
        except Exception:
            pass
        
        # If we got here, the process is not running, so remove lock
        if os.path.exists(lock_path):
            try: os.remove(lock_path)
            except: pass
        return False

    @staticmethod
    def stop_sync():
        """Stops the ongoing sync process (local or PBS)."""
        lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
        if not os.path.exists(lock_path):
            return False, "No sync process found."
        
        try:
            with open(lock_path, "r") as f:
                lock_content = f.read().strip()
            
            if lock_content.startswith("PBS:"):
                job_id = lock_content.replace("PBS:", "")
                subprocess.run(["qdel", job_id], check=True)
                return True, f"Sent qdel for PBS job: {job_id}"
            else:
                # Signal local thread to stop
                stop_path = os.path.join(ROOT_DIR, SyncManager.SYNC_STOP)
                with open(stop_path, "w") as f:
                    f.write("stop")
                return True, "Stop signal sent to local thread."
        except Exception as e:
            return False, f"Error stopping sync: {e}"
        finally:
            if os.path.exists(lock_path):
                try: os.remove(lock_path)
                except: pass

    @staticmethod
    def get_logs(max_lines=100):
        """Reads the latest logs from the sync log file."""
        log_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOG)
        if not os.path.exists(log_path):
            return "No logs found. Start a sync to see progress."
        
        try:
            with open(log_path, "r") as f:
                lines = f.readlines()
                return "".join(lines[-max_lines:])
        except Exception as e:
            return f"Error reading logs: {e}"

    @staticmethod
    def get_running_info():
        """Returns details about the currently running sync process."""
        lock_path = os.path.join(ROOT_DIR, SyncManager.SYNC_LOCK)
        if not os.path.exists(lock_path):
            return None
        try:
            with open(lock_path, "r") as f:
                content = f.read().strip()
            if content.startswith("PBS:"):
                return {"mode": "HPC Job", "id": content.replace("PBS:", "")}
            else:
                return {"mode": "Local Thread", "id": content}
        except:
            return None

class AutoPilotManager:
    LOG_FILE = "Pulse/auto_pilot.log"

    @staticmethod
    def is_running():
        """Checks if the auto_pilot_manager.py script is currently executing."""
        for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
            try:
                cmd = proc.info.get('cmdline')
                if cmd and any("auto_pilot_manager.py" in part for part in cmd):
                    return True
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        return False

    @staticmethod
    def trigger_manual() -> tuple[bool, str]:
        """Triggers a manual run of the auto-pilot manager in the background."""
        if AutoPilotManager.is_running():
            return False, "Auto-Pilot manager is already running."
        
        # Use the current python executable if possible
        python_exe = sys.executable or "/home/guest/miniconda3/envs/gchain/bin/python"
        manager_path = os.path.join(ROOT_DIR, "Pulse", "auto_pilot_manager.py")
        log_path = os.path.join(ROOT_DIR, AutoPilotManager.LOG_FILE)
        
        try:
            # Ensure log file exists and append a separator
            with open(log_path, "a") as f:
                f.write(f"\n--- [MANUAL TRIGGER] {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---\n")
                f.flush()
                # Run as a detached process
                subprocess.Popen(
                    [python_exe, manager_path],
                    stdout=f,
                    stderr=f,
                    cwd=ROOT_DIR,
                    start_new_session=True
                )
            return True, "Auto-Pilot manager triggered in background."
        except Exception as e:
            return False, f"Failed to trigger manager: {str(e)}"
