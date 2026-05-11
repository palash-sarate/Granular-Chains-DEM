import json
import os
import re
from pathlib import Path

# Configuration
OLD_BASE = "/home/guest/palash/Granular-Chains-DEM/dumping_yard"
NEW_BASE = "/Data/palash_data/dumping_yard"

def refactor_json(file_path):
    if not os.path.exists(file_path):
        print(f"Skipping {file_path} (not found)")
        return
    
    print(f"Refactoring {file_path}...")
    with open(file_path, "r") as f:
        content = f.read()
    
    # Replace absolute paths
    new_content = content.replace(OLD_BASE, NEW_BASE)
    
    with open(file_path, "w") as f:
        f.write(new_content)

def refactor_metadata_files():
    print(f"Walking through {NEW_BASE} to update internal metadata...")
    for root, dirs, files in os.walk(NEW_BASE):
        for name in ["grid_metadata.json", "metadata.json"]:
            if name in files:
                p = os.path.join(root, name)
                try:
                    with open(p, "r") as f:
                        data = json.load(f)
                    
                    changed = False
                    # Update source_dir and restart_path
                    for key in ["source_dir", "restart_path"]:
                        if key in data and isinstance(data[key], str) and OLD_BASE in data[key]:
                            data[key] = data[key].replace(OLD_BASE, NEW_BASE)
                            changed = True
                    
                    if changed:
                        with open(p, "w") as f:
                            json.dump(data, f, indent=4)
                except Exception as e:
                    print(f"Warning: Could not refactor {p}: {e}")

if __name__ == "__main__":
    # 1. Update lineage files
    refactor_json("Pulse/lineage.json")
    refactor_json("Pulse/lineage_notes.json")
    
    # 2. Update individual simulation metadata
    refactor_metadata_files()
    
    print("Path migration complete.")
