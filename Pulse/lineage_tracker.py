import os
import sys
import json
import re
from pathlib import Path
import subprocess

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from utilities.env_loader import load_env
LINEAGE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lineage.json")

# Load dynamic environment configuration
env = load_env()
DUMPING_YARD = env.get("DUMPING_YARD", os.path.join(ROOT_DIR, "dumping_yard"))

def get_max_step(run_dir):
    """Finds the maximum timestep reached in the run efficiently using shell streams."""
    chain_dir = os.path.join(run_dir, "chain")
    if not os.path.exists(chain_dir):
        return 0
    try:
        # Use shell streams to find the last file alphabetically (corresponds to max step)
        cmd = f"ls -1f {chain_dir} | grep '^chain_' | sort -V | tail -n 1"
        last_file = subprocess.check_output(cmd, shell=True, text=True).strip()
        if last_file:
            match = re.search(r'chain_(\d+)\.dump', last_file)
            if match:
                return int(match.group(1))
    except Exception as e:
        pass
    return 0

def get_max_step_generalized(run_dir, outputs_spec=None):
    """Finds the maximum timestep reached inside a run directory dynamically based on output specs or fallbacks."""
    if outputs_spec and "dump_files" in outputs_spec:
        dump_patterns = outputs_spec["dump_files"]
    else:
        # Legacy fallback
        dump_patterns = ["chain/chain_*.dump"]
        
    max_step = 0
    for pattern in dump_patterns:
        rel_dir = os.path.dirname(pattern)
        file_glob = os.path.basename(pattern)
        
        target_dir = os.path.join(run_dir, rel_dir) if rel_dir else run_dir
        if not os.path.exists(target_dir):
            continue
            
        prefix = file_glob.split("*")[0]
        suffix = file_glob.split("*")[-1] if "*" in file_glob else ""
        
        try:
            files = os.listdir(target_dir)
            for f in files:
                if f.startswith(prefix) and f.endswith(suffix):
                    core = f[len(prefix):]
                    if suffix:
                        core = core[:-len(suffix)]
                    match = re.search(r'(\d+)', core)
                    if match:
                        step_val = int(match.group(1))
                        if step_val > max_step:
                            max_step = step_val
        except Exception:
            pass
            
    # Fallback to parsing lammps.log
    if max_step == 0:
        log_files = []
        if outputs_spec and "log_file" in outputs_spec:
            log_files.append(outputs_spec["log_file"])
        log_files.extend(["lammps.log", "lammps_resume.log"])
        
        for log_file in log_files:
            log_path = os.path.join(run_dir, log_file)
            if os.path.exists(log_path):
                try:
                    with open(log_path, "r", encoding="utf-8") as lf:
                         lines = lf.readlines()
                    for line in reversed(lines):
                        parts = line.split()
                        if parts and parts[0].isdigit():
                            step_val = int(parts[0])
                            if step_val > max_step:
                                max_step = step_val
                            break
                except Exception:
                    pass
    return max_step

def get_params_from_in_file(run_dir):
    """Parses dt, freq, and amp from the LAMMPS input file."""
    params = {"dt": None, "freq": None, "amp": None}
    try:
        in_files = [f for f in os.listdir(run_dir) if f.startswith("in.")]
        if not in_files:
            return params
        
        in_file_path = os.path.join(run_dir, in_files[0])
        with open(in_file_path, 'r') as f:
            content = f.read()
            
        for key in params.keys():
            match = re.search(fr'variable\s+{key}\s+string\s+([\w.e-]+)', content)
            if match:
                params[key] = match.group(1)
    except Exception:
        pass
    return params

def get_seed_from_name(name):
    """Extracts seed from name like Grid_Fill_..._S123456."""
    match = re.search(r'_S(\d+)$', name)
    if match:
        return match.group(1)
    match = re.search(r'_S(\d+)_', name)
    if match:
        return match.group(1)
    return None

def scan_dumping_yard():
    """Scans the dumping yard and updates the lineage JSON with ancestral inheritance."""
    lineage = {}
    if os.path.exists(LINEAGE_FILE):
        try:
            with open(LINEAGE_FILE, 'r') as f:
                lineage = json.load(f)
        except Exception:
            lineage = {}
            
    # Mark all currently known runs as 'Archived'.
    for run_id in lineage:
        lineage[run_id]["status"] = "Archived"

    # Walk through the dumping yard
    for root, dirs, files in os.walk(DUMPING_YARD, followlinks=True):
        metadata_file = None
        if "grid_metadata.json" in files:
            metadata_file = "grid_metadata.json"
        elif "metadata.json" in files:
            metadata_file = "metadata.json"
            
        if metadata_file:
            dirs[:] = [] 
            run_dir = os.path.realpath(root)
            run_name = os.path.basename(run_dir)
            
            try:
                with open(os.path.join(run_dir, metadata_file), 'r') as f:
                    data = json.load(f)
                
                parent_path = data.get("parent") or data.get("source_dir") or data.get("restart_path")
                if parent_path:
                    parent_path = os.path.realpath(parent_path)
                    if parent_path == run_dir:
                        parent_path = None
                
                n_val = data.get("N", 0)
                if isinstance(n_val, list) and len(n_val) > 0:
                    n_val = n_val[0]
                
                in_params = get_params_from_in_file(run_dir)
                seed = get_seed_from_name(run_name)
                # Normalize rel_path to keep zip names clean
                if run_dir.startswith("/Data/palash_data/dumping_yard"):
                    rel_path = run_dir.replace("/Data/palash_data/", "")
                else:
                    rel_path = os.path.relpath(run_dir, os.getcwd())

                n_fill = data.get("n_fill")
                if isinstance(n_fill, list) and len(n_fill) > 0:
                    n_fill = n_fill[0]

                # Progress parsing
                outputs_spec = data.get("outputs")
                steps_reached = get_max_step_generalized(run_dir, outputs_spec)

                run_params = {
                    "seed": seed or data.get("seed"),
                    "freq": in_params["freq"] or data.get("freq"),
                    "amp": in_params["amp"] or data.get("amp"),
                    "dt": in_params["dt"] or data.get("dt"),
                    "n_hoppers": data.get("n_hoppers", 1),
                    "n_fill": n_fill,
                    "geometry_vars": data.get("geometry_vars", {})
                }
                
                # Merge all generalized overrides
                if "params" in data and isinstance(data["params"], dict):
                    for k, v in data["params"].items():
                        if k not in run_params:
                            run_params[k] = v

                run_info = {
                    "name": run_name,
                    "parent": parent_path,
                    "simulation": data.get("simulation", data.get("simulation_type", "Unknown")),
                    "N": n_val,
                    "steps": steps_reached,
                    "path": run_dir,
                    "rel_path": rel_path,
                    "status": "Active",
                    "params": run_params
                }
                
                if run_dir in lineage:
                    run_info["sync_status"] = lineage[run_dir].get("sync_status", "Local")
                    run_info["sync_time"] = lineage[run_dir].get("sync_time", None)
                else:
                    run_info["sync_status"] = "Local"
                    run_info["sync_time"] = None

                lineage[run_dir] = run_info
                print(f"DEBUG: Parsed {run_name} with N={n_val}, steps={steps_reached}")
                
            except Exception as e:
                print(f"Error parsing metadata for {run_name}: {e}")

    # --- ANCESTRAL INHERITANCE PASS ---
    missing_parents = {}
    for rid, info in lineage.items():
        parent_path = info.get("parent")
        if parent_path:
            # If parent is missing OR exists as a corrupted/archived placeholder with N=0
            if parent_path not in lineage or (lineage[parent_path].get("status") == "Archived" and lineage[parent_path].get("N", 0) == 0):
                parent_name = os.path.basename(parent_path)
                parent_seed = get_seed_from_name(parent_name)
                
                # Normalize rel_path to keep zip names clean
                if parent_path.startswith("/Data/palash_data/dumping_yard"):
                    rel_path = parent_path.replace("/Data/palash_data/", "")
                else:
                    try: rel_path = os.path.relpath(parent_path, os.getcwd())
                    except: rel_path = parent_path
                
                missing_parents[parent_path] = {
                    "name": parent_name,
                    "parent": None, 
                    "simulation": "Hopper_Fill" if "Hopper_Fill" in parent_path else "Unknown", 
                    "N": info.get("N", 0), 
                    "steps": 0, 
                    "path": parent_path,
                    "rel_path": rel_path,
                    "status": "Archived",
                    "params": {"seed": parent_seed},
                    "sync_status": "Synced",
                    "sync_time": None
                }
    
    lineage.update(missing_parents)

    with open(LINEAGE_FILE, 'w') as f:
        json.dump(lineage, f, indent=4)
    
    active_count = len([r for r in lineage.values() if r["status"] == "Active"])
    archived_count = len(lineage) - active_count
    print(f"Lineage scan complete. Total runs: {len(lineage)} ({active_count} active, {archived_count} archived).")

    return lineage

if __name__ == "__main__":
    scan_dumping_yard()
