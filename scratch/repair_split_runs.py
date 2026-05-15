import os
import shutil
import json
import re
from pathlib import Path

# Paths
DUMPING_YARD = Path("/Data/palash_data/dumping_yard")
SRC_DIR = DUMPING_YARD / "Grid_Hopper_Filling"
TARGET_ROOTS = [DUMPING_YARD / "Hopper_Fill", DUMPING_YARD / "Hopper_Fill_Resume"]
REPAIR_DIR = DUMPING_YARD / "Repair_merge"

def get_max_step(run_dir):
    """Finds the maximum timestep reached in the run by scanning dump files."""
    chain_dir = run_dir / "chain"
    if not chain_dir.exists():
        return 0
    max_ts = 0
    for f in chain_dir.glob("chain_*.dump"):
        match = re.search(r'chain_(\d+)\.dump', f.name)
        if match:
            max_ts = max(max_ts, int(match.group(1)))
    return max_ts

def repair():
    if not REPAIR_DIR.exists():
        REPAIR_DIR.mkdir(parents=True)
        print(f"Created temporary repair directory: {REPAIR_DIR}")

    # List all accidental fragments
    fragments = [d for d in SRC_DIR.iterdir() if d.is_dir()]
    print(f"Found {len(fragments)} fragments in {SRC_DIR}")

    repaired_count = 0
    for frag in fragments:
        run_name = frag.name
        parent_dir = None
        
        # Find the original parent
        for root in TARGET_ROOTS:
            potential = root / run_name
            if potential.exists():
                parent_dir = potential
                break
        
        if not parent_dir:
            # If no parent found, we just treat it as a new run to be copied to repair dir
            print(f"[NEW] {run_name} - No parent found in Hopper_Fill or Hopper_Fill_Resume. Just copying.")
            dest_dir = REPAIR_DIR / run_name
            if dest_dir.exists(): shutil.rmtree(dest_dir)
            shutil.copytree(frag, dest_dir)
            repaired_count += 1
            continue

        print(f"[MERGE] {run_name} - Merging {frag.relative_to(DUMPING_YARD)} into {parent_dir.relative_to(DUMPING_YARD)}")
        
        # 1. Create a fresh copy of the parent in the repair directory
        dest_dir = REPAIR_DIR / run_name
        if dest_dir.exists():
            shutil.rmtree(dest_dir)
        shutil.copytree(parent_dir, dest_dir)
        
        # 2. Layer the fragment on top (Merge subdirectories)
        for sub in ["chain", "bond", "angle", "restart", "split_states"]:
            s_sub = frag / sub
            d_sub = dest_dir / sub
            if s_sub.exists():
                if sub == "split_states":
                    # For split_states, we usually just want the newest ones
                    if d_sub.exists(): shutil.rmtree(d_sub)
                    shutil.copytree(s_sub, d_sub)
                else:
                    d_sub.mkdir(parents=True, exist_ok=True)
                    for f in s_sub.iterdir():
                        if f.is_file():
                            shutil.copy2(f, d_sub / f.name)
        
        # 3. Copy root files (.data, .log, .inc, .json)
        for f in frag.iterdir():
            if f.is_file():
                # For lammps.log and input scripts, preserve both if possible or just use newest
                if f.name.endswith(".log") or f.name.startswith("in."):
                    shutil.copy2(f, dest_dir / f"{f.stem}_fragment{f.suffix}")
                else:
                    shutil.copy2(f, dest_dir / f.name)

        # 4. Consolidate Metadata
        meta_data = None
        meta_path = None
        for name in ["metadata.json", "grid_metadata.json"]:
            p = dest_dir / name
            if p.exists():
                with open(p, 'r') as f:
                    meta_data = json.load(f)
                    meta_path = p
                break
        
        if meta_data:
            # Update step count
            new_max = get_max_step(dest_dir)
            meta_data["steps"] = new_max
            
            # Ensure simulation name in metadata matches parent context if desired
            # (Though keeping it as is in parent is safer for this stage)
            
            with open(meta_path, 'w') as f:
                json.dump(meta_data, f, indent=4)
            print(f"      -> Updated steps to {new_max} in {meta_path.name}")
        else:
            print(f"      -> Warning: No metadata file found in {dest_dir} to update.")

        repaired_count += 1

    print(f"\nRepair complete. {repaired_count} runs prepared in {REPAIR_DIR}")
    print("Please verify the contents before manually moving them to their final destination.")

if __name__ == "__main__":
    repair()
