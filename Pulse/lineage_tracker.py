import os
import json
import re
from pathlib import Path

LINEAGE_FILE = "Pulse/lineage.json"
DUMPING_YARD = "dumping_yard"

def get_max_step(run_dir):
    """Finds the maximum timestep reached in the run."""
    chain_dir = os.path.join(run_dir, "chain")
    if not os.path.exists(chain_dir):
        return 0
    max_step = 0
    try:
        for f in os.listdir(chain_dir):
            match = re.match(r'chain_(\d+)\.dump', f)
            if match:
                max_step = max(max_step, int(match.group(1)))
    except Exception:
        pass
    return max_step

def get_params_from_in_file(run_dir):
    """Parses dt, freq, and amp from the LAMMPS input file."""
    params = {"dt": None, "freq": None, "amp": None}
    try:
        # Find the input file (starts with 'in.')
        in_files = [f for f in os.listdir(run_dir) if f.startswith("in.")]
        if not in_files:
            return params
        
        in_file_path = os.path.join(run_dir, in_files[0])
        with open(in_file_path, 'r') as f:
            content = f.read()
            
        # Parse variables like 'variable dt string 1e-06'
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
    # Also handle intermediate seeds in chains
    match = re.search(r'_S(\d+)_', name)
    if match:
        return match.group(1)
    return None

def scan_dumping_yard():
    """Scans the dumping yard and updates the lineage JSON."""
    lineage = {}
    if os.path.exists(LINEAGE_FILE):
        try:
            with open(LINEAGE_FILE, 'r') as f:
                lineage = json.load(f)
        except Exception:
            lineage = {}
            
    # Mark all currently known runs as 'Archived'. We will unmark them if found on disk.
    for run_id in lineage:
        lineage[run_id]["status"] = "Archived"

    # Walk through the dumping yard
    for root, dirs, files in os.walk(DUMPING_YARD):
        metadata_file = None
        if "grid_metadata.json" in files:
            metadata_file = "grid_metadata.json"
        elif "metadata.json" in files:
            metadata_file = "metadata.json"
            
        if metadata_file:
            run_dir = os.path.abspath(root)
            run_name = os.path.basename(run_dir)
            
            try:
                with open(os.path.join(run_dir, metadata_file), 'r') as f:
                    data = json.load(f)
                
                # Determine parent path (Source of truth for lineage)
                parent_path = data.get("source_dir") or data.get("restart_path")
                if parent_path:
                    parent_path = os.path.abspath(parent_path)
                    # If parent is just the current dir (rare bug), set to None
                    if parent_path == run_dir:
                        parent_path = None
                
                # Extract N (handles different metadata formats)
                n_val = data.get("N", 0)
                if isinstance(n_val, list) and len(n_val) > 0:
                    n_val = n_val[0]
                
                # Parse additional params from .in file and name
                in_params = get_params_from_in_file(run_dir)
                seed = get_seed_from_name(run_name)
                
                # Get relative path for cleaner display
                rel_path = os.path.relpath(run_dir, os.getcwd())

                # Build run info
                run_info = {
                    "name": run_name,
                    "parent": parent_path,
                    "simulation": data.get("simulation", "Unknown"),
                    "N": n_val,
                    "steps": get_max_step(run_dir),
                    "path": run_dir,
                    "rel_path": rel_path,
                    "status": "Active",
                    "params": {
                        "seed": seed or data.get("seed"),
                        "freq": in_params["freq"] or data.get("freq"),
                        "amp": in_params["amp"] or data.get("amp"),
                        "dt": in_params["dt"] or data.get("dt")
                    }
                }
                
                # If we already had this run and it was archived, we update it
                lineage[run_dir] = run_info
                
            except Exception as e:
                print(f"Warning: Failed to parse metadata in {run_dir}: {e}")

    # Ensure the directory exists
    os.makedirs(os.path.dirname(LINEAGE_FILE), exist_ok=True)
    
    with open(LINEAGE_FILE, 'w') as f:
        json.dump(lineage, f, indent=4)
    
    active_count = len([r for r in lineage.values() if r["status"] == "Active"])
    archived_count = len(lineage) - active_count
    print(f"Lineage scan complete. Total runs: {len(lineage)} ({active_count} active, {archived_count} archived).")

if __name__ == "__main__":
    scan_dumping_yard()
