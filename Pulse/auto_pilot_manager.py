import os
import json
import subprocess
import datetime
import time

# Paths
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
LINEAGE_FILE = os.path.join(ROOT_DIR, "Pulse", "lineage.json")

QSTAT_PATH = "/opt/pbs/bin/qstat"
QSUB_PATH = "/opt/pbs/bin/qsub"

PBS_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -l nodes=master:ppn=16
#PBS -l mem=16gb
#PBS -m n
#PBS -o {log_dir}/{job_name}.log
#PBS -e {log_dir}/{job_name}_err.log
#PBS -S /bin/bash

# 1. Set working directory
cd $PBS_O_WORKDIR

# 2. Activate environment
source /home/guest/miniconda3/etc/profile.d/conda.sh
conda activate gchain

# 3. Launch Simulation
echo "Starting simulation: {job_name} at: $(date)"

{command}

echo "=============================================================================="
echo "Simulation complete at: $(date)"
echo "=============================================================================="
"""

def get_pbs_jobs():
    """Returns a list of all jobs currently in the PBS queue."""
    try:
        # Standard qstat shows all jobs on this cluster
        result = subprocess.check_output([QSTAT_PATH], text=True)
        lines = result.strip().split("\n")
        jobs = []
        for line in lines:
            if "---" in line or "Job id" in line: continue
            parts = line.split()
            if len(parts) >= 5:
                jobs.append({
                    "id": parts[0],
                    "name": parts[1],
                    "user": parts[2],
                    "status": parts[4] if len(parts) > 4 else "?"
                })
        return jobs
    except Exception as e:
        print(f"Error fetching PBS jobs: {e}")
        return []

def is_student_active(jobs):
    return any(j["user"] != "guest" for j in jobs)

def count_guest_jobs(jobs):
    return len([j for j in jobs if j["user"] == "guest"])

def run_manager():
    # 1. Update Heartbeat
    if not os.path.exists(AUTO_PILOT_FILE):
        return
        
    with open(AUTO_PILOT_FILE, "r") as f:
        data = json.load(f)
    
    data["settings"]["last_heartbeat"] = datetime.datetime.now().isoformat()
    
    if not data["settings"]["enabled"]:
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)
        return

    # 2. Check for Students (Polite Mode)
    jobs = get_pbs_jobs()
    if data["settings"]["polite_mode"] and is_student_active(jobs):
        print("Student detected. Standing down.")
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)
        return

    # 3. Check Concurrency
    active_guest_jobs = count_guest_jobs(jobs)
    max_concurrent = data["settings"].get("max_concurrent", 3)
    
    if active_guest_jobs >= max_concurrent:
        print(f"Concurrency limit reached ({active_guest_jobs}/{max_concurrent})")
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)
        return

    slots_available = max_concurrent - active_guest_jobs

    # 4. Fetch latest lineage
    with open(LINEAGE_FILE, "r") as f:
        lineage = json.load(f)

    # 5. Round Robin Selection
    goals = data.get("goals", {})
    eligible = []
    for path, g in goals.items():
        if path not in lineage: continue
        curr_steps = lineage[path].get("steps", 0)
        if curr_steps < g["target_steps"]:
            # Check if this run is already in queue by name
            run_name = os.path.basename(path)
            if any(run_name in j["name"] for j in jobs):
                continue
            
            eligible.append((path, g, curr_steps))

    # Sort by last_submitted
    eligible.sort(key=lambda x: x[1].get("last_submitted", ""))

    # 6. Submit
    log_dir = os.path.join(ROOT_DIR, "PBS_Output")
    os.makedirs(log_dir, exist_ok=True)
    
    for i in range(min(slots_available, len(eligible))):
        path, g, curr_steps = eligible[i]
        run_name = os.path.basename(path)
        # Generate Command
        mode = g.get("mode", "fill_resume")
        
        # Map modes to main.py commands
        mapping = {
            "fill": "run_grid_hopper_filling",
            "fill_resume": "resume_grid_hopper_filling",
            "flow": "run_grid_hopper_flow",
            "flow_resume": "resume_grid_hopper_flow"
        }
        cmd_name = mapping.get(mode, mode)
        
        cmd_parts = ["python", "-u", "main.py", cmd_name]
        
        # Handle specific mode paths
        if "flow" in mode and "resume" not in mode:
            # Flow study takes source_dir
            cmd_parts.extend(["--source_dir", path])
        else:
            # Resumes take restart_path
            cmd_parts.extend(["--restart_path", path])

        cmd_parts.extend([
            "--relax_steps" if "fill" in mode else "--run_steps", str(g["increment"]),
            "--inplace" if g.get("in_place", True) else ""
        ])
        
        # Add any overrides (dt, viscosity, num_procs, etc.)
        overrides = g.get("params", {})
        if overrides:
            for k, v in overrides.items():
                if v is not None and v != "":
                    cmd_parts.extend([f"--{k}", str(v)])
        
        full_command = " ".join([c for c in cmd_parts if c])

        # Generate PBS script
        job_name = f"AP_{run_name}"
        script_content = PBS_TEMPLATE.format(
            job_name=job_name,
            log_dir=log_dir,
            command=full_command
        )
        
        script_path = os.path.join(ROOT_DIR, "Pulse", f"temp_{job_name}.pbs")
        with open(script_path, "w") as f:
            f.write(script_content)
        
        try:
            res = subprocess.run([QSUB_PATH, script_path], capture_output=True, text=True, check=True)
            job_id = res.stdout.strip()
            print(f"Submitted {job_name} as {job_id}")
            
            data["goals"][path]["last_submitted"] = datetime.datetime.now().isoformat()
            data["goals"][path]["status"] = "Running"
        except Exception as e:
            print(f"Failed to submit {job_name}: {e}")
        finally:
            if os.path.exists(script_path): os.remove(script_path)

    with open(AUTO_PILOT_FILE, "w") as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    run_manager()
