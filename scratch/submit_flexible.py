import os
import json
import subprocess
import argparse
import random
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
    # Currently resumed
    # ----------------------------------------------------------------------
    # resume_source_dirs = {
    #     4:  "dumping_yard/Hopper_Fill/Grid_Fill_1H_S808086",
    #     12: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S481365",
    #     24: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S469898",
    #     48: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S437513"
    # }
    # ----------------------------------------------------------------------
    # Next to resume
    # ----------------------------------------------------------------------
    resume_source_dirs = {
        4:  "dumping_yard/Hopper_Fill/Grid_Fill_1H_S374952",
        12: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S935840",
        24: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S240587",
        48: "dumping_yard/Hopper_Fill/Grid_Fill_1H_S358667"
    }

    for n in Ns:
        n_fill = n_atoms // n
        
        match mode:
            case "fill":
                # ----------------------------------------------------------------------
                # TYPE 1: FILL (New Hopper Filling Run)
                # ----------------------------------------------------------------------
                jobs.append({
                    "name": f"Fill_N{n}",
                    "type": "fill",
                    "walltime": "48:00:00",
                    "ppn": 16,
                    "mem": "16gb",
                    "params": {
                        "num_procs": 8,
                        "num_threads": 1,
                        "N": n,
                        "seed": random.randint(100000, 999999),
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
                    jobs.append({
                        "name": f"Resume_Fill_N{n}",
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
                            "seed": random.randint(100000, 999999),
                        }
                    })

            case "flow":
                if n in resume_source_dirs and os.path.exists(resume_source_dirs[n]):
                    for f in freqs:
                        for a in amps:
                            jobs.append({
                                "name": f"Flow_N{n}_F{f}_A{a}",
                                "type": "flow",
                                "walltime": "48:00:00",
                                "params": {
                                    "source_dir": resume_source_dirs[n],
                                    "seed": random.randint(100000, 999999),
                                    "freq": f,
                                    "amp": a,
                                    "run_steps": 2000000,
                                    "osc_dir": "z",
                                    "simulation": f"Flow_Study_N{n}_F{f}_A{a}"
                                }
                            })
            
            case "flow_resume":
                if n in resume_source_dirs and os.path.exists(resume_source_dirs[n]):
                    for f in freqs:
                        for a in amps:
                            jobs.append({
                                "name": f"Res_Flow_N{n}_F{f}_A{a}",
                                "type": "flow_resume",
                                "walltime": "24:00:00",
                                "params": {
                                    "restart_path": resume_source_dirs[n],
                                    "run_steps": 1000000,
                                    "dt": 1e-06,
                                    "seed": random.randint(100000, 999999),
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

def main():
    parser = argparse.ArgumentParser(description="Flexible PBS Job Submitter")
    parser.add_argument("--mode", choices=["fill", "fill_resume", "flow", "flow_resume"], required=True, help="Mode of simulation study")
    parser.add_argument("--submit", action="store_true", help="Submit jobs to the queue (otherwise only generates files)")
    args = parser.parse_args()

    jobs = get_job_list(args.mode)
    
    if not os.path.exists("PBS_Output"): os.makedirs("PBS_Output")
    if not os.path.exists("temp/temp_pbs"): os.makedirs("temp/temp_pbs")

    for job in jobs:
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
            
            if args.submit:
                print(f"Submitting: {name}")
                subprocess.run(["qsub", pbs_file])
            else:
                print(f"Generated: {pbs_file}")
        except Exception as e:
            print(f"Error on {name}: {e}")

if __name__ == "__main__":
    main()
