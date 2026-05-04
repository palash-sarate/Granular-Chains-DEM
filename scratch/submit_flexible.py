import os
import json
import subprocess
import argparse
import random
import re
from datetime import datetime

# --- 1. CONFIGURATION & TEMPLATE ---

PBS_TEMPLATE = """#!/bin/bash
#PBS -N {job_name}
#PBS -q workq
#PBS -l walltime={walltime}
#PBS -l nodes=master:ppn={ppn}
#PBS -l mem={mem}
#PBS -m aeb
#PBS -M palashsarate@iisc.ac.in
#PBS -o PBS_Output/{job_name}.log
#PBS -e PBS_Output/{job_name}_err.log
#PBS -S /bin/bash

# 1. Set working directory
cd $PBS_O_WORKDIR

# Clean old logs if they exist
rm -f PBS_Output/{job_name}.log
rm -f PBS_Output/{job_name}_err.log

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

DEFAULTS = {
    "walltime": "24:00:00",
    "ppn": 16,
    "mem": "16gb",
}

# --- 2. DEFINE YOUR STUDY HERE ---

def get_job_list(mode):
    """
    Define your parametric study here. 
    Selects the job generation logic based on the provided mode.
    """
    jobs = []
    
    # Global Parameters
    Ns = [4, 12, 24, 48]
    n_atoms = 14400
    
    # Frequencies/Amplitudes for Flow Study
    freqs = [5.0, 10.0, 15.0]
    amps = [0.01]

    # Map N to the directory containing the results you want to resume or flow from
    # Replace these paths with actual simulation output directories
    # ----------------------------------------------------------------------
    # Currently running resumes
    # ----------------------------------------------------------------------
    # resume_source_dirs = {
    #     4:  "dumping_yard/Hopper_Fill/Grid_Fill_1H_S808086",
    #     12: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S481365",
    #     24: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S469898",
    #     48: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S437513"
    # }
    # ----------------------------------------------------------------------
    # On Hold resumes
    # ----------------------------------------------------------------------
    # resume_source_dirs = {
    #     4:  "dumping_yard/Hopper_Fill/Grid_Fill_1H_S374952",
    #     12: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S935840",
    #     24: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S240587",
    #     48: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S358667"
    #     48: "dumping_yard/Hopper_Fill_Resume/Grid_Fill_1H_S437513_S404527",
    #     24: "dumping_yard/Hopper_Fill_Resume/Grid_Fill_1H_S469898_S180256",
    #     12: "dumping_yard/Hopper_Fill_Resume/Grid_Fill_1H_S481365_S956135",
    #     4:  "dumping_yard/Hopper_Fill_Resume/Grid_Fill_1H_S808086_S234775" 
    # }
    # ----------------------------------------------------------------------
    # Next to resume
    # ----------------------------------------------------------------------
    resume_source_dirs = {
        # 48: "",
        # 24: "",
        # 12: "",
        4:  "dumping_yard/Hopper_Fill/Grid_Fill_1H_S374952"        
    }
    for n in Ns:
        n_fill = n_atoms // n
        
        match mode:
            case "fill":
                # ----------------------------------------------------------------------
                # TYPE 1: FILL (New Hopper Filling Run)
                # ----------------------------------------------------------------------
                child_seed = random.randint(100000, 999999)
                jobs.append({
                    "name": f"S000000_{child_seed}",
                    "type": "fill",
                    "walltime": "48:00:00",
                    "ppn": 16,
                    "mem": "16gb",
                    "params": {
                        "num_procs": 8,
                        "num_threads": 1,
                        "N": n,
                        "seed": child_seed,
                        "n_fill": n_fill,
                        "n_hoppers": 1,
                        "relax_steps": 1000000,
                        "simulation": "Hopper_Fill",
                        "dt": 1e-06,
                        "viscosity": 0.001,
                        "dump_file": "simulation_templates/default_dump.inc",
                        "no-vtk": True,
                        "spacing": 0.5,
                        "mode": "2D_stacked",
                        "source_dir": "chain_data/relaxed_2D_x",
                        "hopper_template_data": "simulation_geometries/2D_hopper_with_orifice_cover.inc"
                    }
                })

            case "fill_resume":
                # ----------------------------------------------------------------------
                # TYPE 2: FILL_RESUME (Continue or add relaxation to a filling run)
                # ----------------------------------------------------------------------
                if n in resume_source_dirs and os.path.exists(resume_source_dirs[n]):
                    match = re.search(r'S(\d+)', resume_source_dirs[n])
                    parent_seed = match.group(1)[-6:] if match else "000000"
                    child_seed = random.randint(100000, 999999)
                    jobs.append({
                        "name": f"R{parent_seed}_{child_seed}",
                        "type": "fill_resume",
                        "walltime": "48:00:00",
                        "params": {
                            "num_procs": 8,
                            "num_threads": 1,
                            "restart_path": f"{resume_source_dirs[n]}",
                            "relax_steps": 1000000,
                            "simulation": f"Hopper_Fill_Resume",
                            "dt": 1e-06,
                            "viscosity": 0.001,
                            "dump_file": "simulation_templates/default_dump.inc",
                            "seed": child_seed,
                        }
                    })

            case "flow":
                if n in resume_source_dirs and os.path.exists(resume_source_dirs[n]):
                    match = re.search(r'S(\d+)', resume_source_dirs[n])
                    parent_seed = match.group(1)[-6:] if match else "000000"
                    for f in freqs:
                        for a in amps:
                            child_seed = random.randint(100000, 999999)
                            jobs.append({
                                "name": f"F{parent_seed}_{child_seed}",
                                "type": "flow",
                                "walltime": "48:00:00",
                                "params": {
                                    "source_dir": resume_source_dirs[n],
                                    "seed": child_seed,
                                    "freq": f,
                                    "amp": a,
                                    "run_steps": 2000000,
                                    "osc_dir": "z",
                                    "simulation": f"Flow_Study_N{n}_F{f}_A{a}"
                                }
                            })
            
            case "flow_resume":
                if n in resume_source_dirs and os.path.exists(resume_source_dirs[n]):
                    match = re.search(r'S(\d+)', resume_source_dirs[n])
                    parent_seed = match.group(1)[-6:] if match else "000000"
                    for f in freqs:
                        for a in amps:
                            child_seed = random.randint(100000, 999999)
                            jobs.append({
                                "name": f"FR{parent_seed}_{child_seed}",
                                "type": "flow_resume",
                                "walltime": "24:00:00",
                                "params": {
                                    "restart_path": resume_source_dirs[n],
                                    "run_steps": 1000000,
                                    "dt": 1e-06,
                                    "seed": child_seed,
                                }
                            })

    return jobs

# --- 3. SUBMISSION ENGINE ---

def generate_command(config):
    cmd_type = config.get("type", "run")
    # Merge defaults with specific params
    full_params = DEFAULTS.copy()
    full_params.update(config.get("params", {}))
    
    mapping = {
        "fill": "run_grid_hopper_filling",
        "fill_resume": "resume_grid_hopper_filling",
        "flow": "run_grid_hopper_flow",
        "flow_resume": "resume_grid_hopper_flow"
    }
    
    command_name = mapping.get(cmd_type, cmd_type)
    base_cmd = f"python -u main.py {command_name}"
    
    cmd_parts = [base_cmd]
    
    # Parameters to ignore for the CLI (they go to PBS header instead)
    ignore_keys = ["walltime", "ppn", "mem", "type", "name"]
    
    for key, val in full_params.items():
        if key in ignore_keys: continue
        
        if isinstance(val, bool):
            if val:
                # Special case for no-vtk which is a store_false action
                if key == "no-vtk":
                    cmd_parts.append("--no-vtk")
                else:
                    cmd_parts.append(f"--{key}")
        elif val is not None:
            if isinstance(val, (dict, list)):
                # Auto-serialize complex structures to JSON
                cmd_parts.append(f"--{key} '{json.dumps(val)}'")
            else:
                cmd_parts.append(f"--{key} {val}")
    
    return " \\\n    ".join(cmd_parts)
    
def get_active_job_ids(user):
    """
    Fetch current active job IDs for the user.
    """
    try:
        result = subprocess.run(["qstat", "-u", user], capture_output=True, text=True, stdin=subprocess.DEVNULL)
        if result.returncode != 0:
            return []
        lines = result.stdout.splitlines()
        ids = []
        for line in lines:
            parts = line.split()
            # Expecting ID in the first column, usually like '3281.master'
            if len(parts) > 0 and "." in parts[0] and user in line:
                ids.append(parts[0])
        return ids
    except Exception as e:
        try: print(f"Warning: Could not fetch active jobs: {e}")
        except OSError: pass
        return []

def safe_print(msg):
    try:
        print(msg)
    except OSError:
        pass

def submit_jobs(jobs, submit=True, max_concurrent=4, user="guest"):
    """
    Takes a list of job configurations, generates PBS scripts, and optionally submits them.
    Returns a list of generated PBS files and a list of submitted job IDs (if submit=True).
    """
    tails = []
    if submit:
        tails = get_active_job_ids(user)
        # Limit to the most recent max_concurrent jobs if there are already many
        if len(tails) > max_concurrent:
            tails = tails[-max_concurrent:]
        safe_print(f"Current active jobs detected: {len(tails)}. Limit: {max_concurrent}")

    if not os.path.exists("PBS_Output"): os.makedirs("PBS_Output")
    if not os.path.exists("temp/temp_pbs"): os.makedirs("temp/temp_pbs")

    generated_files = []
    submitted_ids = []

    for i, job in enumerate(jobs):
        name = job.get("name", "Job_" + datetime.now().strftime("%H%M%S"))
        walltime = job.get("walltime", DEFAULTS["walltime"])
        ppn = job.get("ppn", DEFAULTS["ppn"])
        mem = job.get("mem", DEFAULTS["mem"])
        
        try:
            command = generate_command(job)
            pbs_content = PBS_TEMPLATE.format(
                job_name=name,
                walltime=walltime,
                ppn=ppn,
                mem=mem,
                command=command
            )
            
            pbs_file = f"temp/temp_pbs/{name}.pbs"
            with open(pbs_file, "w") as f:
                f.write(pbs_content)
            
            generated_files.append(pbs_file)
            
            if submit:
                # 1. Lineage/Manual Dependency (Specific parent job)
                manual_dep = job.get("dependency")
                
                # 2. Concurrency Dependency (Round-robin to stay under max_concurrent)
                concurrency_dep = None
                if len(tails) >= max_concurrent:
                    concurrency_dep = tails[i % max_concurrent]
                
                # Combine dependencies (PBS supports afterany:id1:id2)
                deps = []
                if manual_dep: deps.append(manual_dep)
                if concurrency_dep and concurrency_dep not in deps: deps.append(concurrency_dep)
                
                cmd = ["qsub"]
                if deps:
                    cmd.extend(["-W", f"depend=afterany:{':'.join(deps)}"])

                cmd.append(pbs_file)
                
                msg = f"Submitting: {name}"
                if deps: msg += f" (depends on {', '.join(deps)})"
                safe_print(msg)

                
                res = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL)
                if res.returncode == 0:
                    new_id = res.stdout.strip()
                    submitted_ids.append(new_id)
                    if len(tails) < max_concurrent:
                        tails.append(new_id)
                    else:
                        tails[i % max_concurrent] = new_id
                else:
                    safe_print(f"Error submitting {name}: {res.stderr.strip()}")
            else:
                safe_print(f"Generated: {pbs_file}")
        except Exception as e:
            safe_print(f"Error on {name}: {e}")
            
    return generated_files, submitted_ids

def main():
    parser = argparse.ArgumentParser(description="Flexible PBS Job Submitter")
    parser.add_argument("--mode", choices=["fill", "fill_resume", "flow", "flow_resume"], required=True, help="Mode of simulation study")
    parser.add_argument("--submit", action="store_true", help="Submit jobs to the queue (otherwise only generates files)")
    parser.add_argument("--max-concurrent", type=int, default=4, help="Max number of concurrent jobs allowed (default: 4)")
    parser.add_argument("--user", default="guest", help="User to check for active jobs (default: guest)")
    args = parser.parse_args()

    jobs = get_job_list(args.mode)
    submit_jobs(jobs, submit=args.submit, max_concurrent=args.max_concurrent, user=args.user)

if __name__ == "__main__":
    main()
