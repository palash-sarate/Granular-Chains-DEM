import os
import json
import subprocess
import datetime
import time
import fcntl

# Paths
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
LINEAGE_FILE = os.path.join(ROOT_DIR, "Pulse", "lineage.json")

QSTAT_PATH = "/opt/pbs/bin/qstat"
QSUB_PATH = "/opt/pbs/bin/qsub"

def save_auto_pilot(data):
    with open(AUTO_PILOT_FILE, "w") as f:
        json.dump(data, f, indent=4)


PBS_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -q workq
#PBS -p {priority}
#PBS -l walltime={walltime}
#PBS -l nodes=master:ppn={ppn}
#PBS -l mem={mem}
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

# 4. Trigger next Auto-Pilot cycle
python3 {manager_script}
"""

def get_pbs_jobs():
    """Returns a list of all jobs currently in the PBS queue."""
    try:
        # Try structured JSON format first (prevents job name truncation!)
        result = subprocess.check_output([QSTAT_PATH, "-f", "-F", "json"], text=True)
        data = json.loads(result)
        jobs = []
        for job_id, details in data.get("Jobs", {}).items():
            # Standardize user name (remove host suffix if present)
            owner = details.get("Job_Owner", "")
            user = owner.split("@")[0] if "@" in owner else owner
            
            jobs.append({
                "id": job_id,
                "name": details.get("Job_Name", ""),
                "user": user,
                "status": details.get("job_state", "?")
            })
        return jobs
    except Exception as json_err:
        try:
            # Fallback: Standard qstat shows all jobs on this cluster (names may be truncated)
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
    # 1. Use a lock file to prevent concurrent manager runs
    lock_path = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.lock")
    lock_f = open(lock_path, "w")
    try:
        fcntl.flock(lock_f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except IOError:
        print("Auto-Pilot is already running. Exiting.")
        return

    # 2. Update Heartbeat
    if not os.path.exists(AUTO_PILOT_FILE):
        return
        
    with open(AUTO_PILOT_FILE, "r") as f:
        data = json.load(f)
    
    data["settings"]["last_heartbeat"] = datetime.datetime.now().isoformat()
    
    if not data["settings"]["enabled"]:
        save_auto_pilot(data)
        return

    # 2. Check for Students (Polite Mode)
    jobs = get_pbs_jobs()
    students_active = is_student_active(jobs)
    global_polite = data["settings"].get("polite_mode", True)
    if students_active and global_polite:
        print("Student detected and global Polite Mode is active. Yielding polite goals.")

    # 3. Check Concurrency
    active_guest_jobs = count_guest_jobs(jobs)
    max_concurrent = data["settings"].get("max_concurrent") or 3
    
    if active_guest_jobs >= max_concurrent:
        print(f"Concurrency limit reached ({active_guest_jobs}/{max_concurrent})")
        save_auto_pilot(data)
        return

    slots_available = max_concurrent - active_guest_jobs

    # 4. Fetch latest lineage
    with open(LINEAGE_FILE, "r") as f:
        lineage = json.load(f)

    # 5. Round Robin Selection
    goals = data.get("goals", {})
    
    # Discovery: Transition NewRoot goals to actual paths once they appear in lineage
    for path, g in list(goals.items()):
        if path.startswith("NewRoot") and g.get("last_submitted"):
            target_seed = g.get("params", {}).get("seed")
            for l_path, l_info in lineage.items():
                l_seed = l_info.get("params", {}).get("seed")
                if l_seed and str(l_seed) == str(target_seed):
                    # Found the newly created run! Hand off the goal to the real path
                    goals[l_path] = g.copy()
                    goals[l_path]["mode"] = "fill_resume" if g["mode"] == "fill" else g["mode"]
                    del goals[path]
                    # Update data store immediately to persist the transition
                    save_auto_pilot(data)
                    break

    eligible = []
    for path, g in goals.items():
        # A goal is polite if the global polite mode is on AND the goal's polite mode is also on (or defaults to on)
        is_goal_polite = global_polite and g.get("polite_mode", True)
        
        if students_active and is_goal_polite:
            # Skip this goal since students are active and this goal is polite
            continue
            
        if path.startswith("NewRoot"):
            if g.get("last_submitted"): continue
            eligible.append((path, g, 0))
            continue
        if path not in lineage: continue
        curr_steps = lineage[path].get("steps") or 0
        target_steps = g.get("target_steps") or 0
        if curr_steps < target_steps:
            # Check if this run is already in queue by name
            run_name = os.path.basename(path)
            if any(run_name in j["name"] for j in jobs):
                continue
            
            eligible.append((path, g, curr_steps))

    # Sort by last_submitted (handle None by converting to empty string)
    eligible.sort(key=lambda x: x[1].get("last_submitted") or "")

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
        if path.startswith("NewRoot"):
            # New fill doesn't have a source path usually, or it's in params
            pass
        elif "flow" in mode and "resume" not in mode:
            # Flow study takes source_dir
            cmd_parts.extend(["--source_dir", path])
        else:
            # Resumes take restart_path
            cmd_parts.extend(["--restart_path", path])

        cmd_parts.extend([
            "--relax_steps" if "fill" in mode else "--run_steps", str(g.get("increment") or 100000),
            "--inplace" if g.get("in_place", True) else ""
        ])
        
        # Ensure a unique seed for NewRoot jobs if not already present
        # This prevents collisions and allows discovery to work
        overrides = g.get("params", {})
        if path.startswith("NewRoot") and "seed" not in overrides:
            try:
                # Use the numeric suffix from the NewRoot path as the seed
                seed_val = int(path.split("_")[-1]) % 1000000
                overrides["seed"] = seed_val
            except (ValueError, IndexError):
                overrides["seed"] = int(time.time() * 1000) % 1000000
        
        # Add any overrides (dt, viscosity, num_procs, etc.)
        # Exclude PBS-specific params and setup-only params for resumes/flows
        skip_params = {"walltime", "ppn", "mem"}
        if "resume" in mode or "flow" in mode:
            # These are handled by the restart file or source_dir
            skip_params.update({
                "N", "n_fill", "spacing", "n_hoppers", "mode", 
                "source_dir", "no-vtk", "geometry_vars", "hopper_template_data"
            })
            if "resume" in mode:
                skip_params.add("restart_path")
            elif mode == "flow":
                # When starting a new flow from a fill, we want to pivot to the 
                # default Grid_Hopper_Flow category rather than inheriting the fill category.
                skip_params.add("simulation")

        overrides = g.get("params", {})
        if overrides:
            for k, v in overrides.items():
                if k not in skip_params and v is not None and v != "":
                    if isinstance(v, (dict, list)):
                        cmd_parts.extend([f"--{k}", f"'{json.dumps(v)}'"])
                    elif isinstance(v, bool):
                        if v: cmd_parts.append(f"--{k}")
                    else:
                        cmd_parts.extend([f"--{k}", str(v)])
        
        full_command = " ".join([c for c in cmd_parts if c])

        # Generate PBS script
        job_name = f"AP_{run_name}"
        goal_polite = g.get("polite_mode", True)
        priority_val = -1024 if goal_polite else 0

        script_content = PBS_TEMPLATE.format(
            job_name=job_name,
            log_dir=log_dir,
            command=full_command,
            manager_script=os.path.abspath(__file__),
            walltime=overrides.get("walltime", "24:00:00"),
            ppn=overrides.get("ppn", 16),
            mem=overrides.get("mem", "16gb"),
            priority=priority_val
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
            
            # Transactional save to prevent duplicate submissions if we crash later
            save_auto_pilot(data)
        except Exception as e:
            print(f"Failed to submit {job_name}: {e}")
        finally:
            if os.path.exists(script_path): os.remove(script_path)

    save_auto_pilot(data)

if __name__ == "__main__":
    run_manager()
