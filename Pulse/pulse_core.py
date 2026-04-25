import os
import subprocess
import re
import json
import time
from datetime import datetime
from typing import List, Dict, Optional
import numpy as np

class PBSManager:
    METADATA_FILE = "Pulse/pulse_metadata.json"

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
        cmd = ["qstat", "-f"]
        if user:
            # Note: qstat -f doesn't always support -u directly in all versions
            # so we fetch all and filter in python for robustness
            pass
            
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                return []
            return PBSManager.parse_qstat_f(result.stdout, user)
        except Exception:
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
                elif last_key and line.startswith("    "):
                    # Continuation line
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
    def estimate_eta(directory: str, target_timestep: int) -> Dict:
        """
        Estimates the ETA for a running simulation based on dump files.
        Logic:
        1. Scan directory for chain_{timestep}.dump files.
        2. Filter for every 10th file (based on sorted order).
        3. Get modification times.
        4. Extrapolate to target_timestep using a quadratic fit (handles slowdown).
        """
        if not os.path.exists(directory):
            return {"error": f"Directory not found: {directory}"}

        # 1. Collect all chain_*.dump files
        files = []
        try:
            for f in os.listdir(directory):
                match = re.match(r'chain_(\d+)\.dump', f)
                if match:
                    timestep = int(match.group(1))
                    files.append((timestep, os.path.getmtime(os.path.join(directory, f))))
        except Exception as e:
            return {"error": f"Error scanning directory: {str(e)}"}

        if not files:
            return {"error": "No dump files found matching 'chain_*.dump'"}

        files.sort() # Sort by timestep

        # 2. Sample data points
        # Every 10th file (0, 10, 20...)
        sampled_files = files[::10]
        
        # Ensure we always include the latest file for better precision
        if files[-1] not in sampled_files:
            sampled_files.append(files[-1])
            sampled_files.sort()

        if len(sampled_files) < 2:
            if len(files) >= 2:
                sampled_files = files # Use all if we don't have enough for sampling
            else:
                return {"error": "Insufficient data: Need at least 2 dump files."}

        # 3. Prepare for fit
        x = np.array([f[0] for f in sampled_files])
        y = np.array([f[1] for f in sampled_files]) # Timestamps in seconds

        current_time = time.time()
        current_timestep = files[-1][0]
        
        if current_timestep >= target_timestep:
            return {
                "status": "Target Reached",
                "current_timestep": current_timestep,
                "target_timestep": target_timestep,
                "progress_percent": 100.0,
                "time_remaining_hr": 0.0,
                "time_elapsed_hr": (files[-1][1] - files[0][1]) / 3600.0,
                "completion_time": datetime.fromtimestamp(files[-1][1]).strftime("%Y-%m-%d %H:%M:%S"),
                "cost_per_10k_steps": 0.0,
                "cost_history": [],
                "data_points": len(sampled_files),
                "last_updated": datetime.now().strftime("%H:%M:%S"),
                "message": "The simulation has already reached or passed the target timestep."
            }

        # 4. Perform Extrapolation
        # Quadratic fit: y = ax^2 + bx + c
        # We use a 2nd degree polynomial because the computational cost often 
        # increases linearly with time in growing systems, making total time quadratic.
        degree = 2 if len(sampled_files) >= 3 else 1
        
        # Fit
        coeffs = np.polyfit(x, y, degree)
        poly = np.poly1d(coeffs)
        
        # Estimate completion time
        estimated_completion_time = poly(target_timestep)
        
        # Safety Check: If quadratic extrapolation goes "backward" or is unrealistically fast 
        # (can happen if the simulation recently sped up), fallback to latest linear rate.
        last_x, last_y = files[-1]
        prev_x, prev_y = files[-min(len(files), 5)] # Compare with something recent
        
        if last_x != prev_x:
            latest_rate = (last_y - prev_y) / (last_x - prev_x)
            linear_eta = last_y + (target_timestep - last_x) * latest_rate
            
            # If quadratic is faster than linear (rare in slowing sims) or in the past, use linear
            if estimated_completion_time < current_time or estimated_completion_time < linear_eta:
                estimated_completion_time = linear_eta

        # Calculate metrics
        time_remaining_sec = max(0, estimated_completion_time - current_time)
        time_remaining_hr = time_remaining_sec / 3600.0
        
        # Time elapsed since the first discovered dump file
        start_time = files[0][1]
        time_elapsed_hr = (current_time - start_time) / 3600.0
        
        # Real world time
        completion_dt = datetime.fromtimestamp(estimated_completion_time)
        
        # Performance trend (Seconds per 10k steps)
        # We look at the derivative of the polynomial at the current timestep
        # poly'(x) = 2ax + b
        deriv = np.polyder(poly)
        current_cost_per_step = deriv(current_timestep)
        
        return {
            "status": "Active",
            "current_timestep": current_timestep,
            "target_timestep": target_timestep,
            "progress_percent": (current_timestep / target_timestep) * 100,
            "time_remaining_hr": time_remaining_hr,
            "time_elapsed_hr": time_elapsed_hr,
            "completion_time": completion_dt.strftime("%Y-%m-%d %H:%M:%S"),
            "cost_per_10k_steps": current_cost_per_step * 10000 / 60.0, # in minutes
            "cost_history": (deriv(x) * 10000 / 60.0).tolist(), # Convert numpy array to list for JSON/session_state
            "data_points": len(sampled_files),
            "last_updated": datetime.now().strftime("%H:%M:%S")
        }
