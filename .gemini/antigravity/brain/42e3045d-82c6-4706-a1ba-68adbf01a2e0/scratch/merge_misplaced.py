import os
import shutil
import json
import glob
from pathlib import Path

RUNS = [
    "Grid_Fill_1H_S358667_S348132_S521842_S161014_S897208_S769729",
    "Grid_Fill_1H_S437513_S402899",
    "Grid_Fill_1H_S442391_S430735_S744817_S842756_S169565",
    "Grid_Fill_1H_S694955_S950784_S371727_S801731_S381262",
    "Grid_Fill_1H_S904266_S944243_S153340_S131626_S155721"
]

BASE_SRC = Path("/Data/palash_data/dumping_yard/Hopper_Fill")
BASE_DST = Path("/Data/palash_data/dumping_yard/Hopper_Fill_Resume")

def merge_run(run_name):
    src = BASE_SRC / run_name
    dst = BASE_DST / run_name
    
    if not src.exists():
        print(f"Skipping {run_name}: Source directory not found.")
        return

    print(f"Merging {run_name}...")
    
    # 1. Sync everything from SRC to DST
    # rsync -av will copy all files and subdirs, preserving times, and overwriting DST with newer SRC files.
    os.system(f"rsync -av {src}/ {dst}/")

    # 2. Update metadata in DST
    meta_path = dst / "grid_metadata.json"
    if meta_path.exists():
        with open(meta_path, 'r') as f:
            meta = json.load(f)
        
        # Find highest restart step
        restarts = glob.glob(str(dst / "restart" / "restart.*.bin"))
        steps = []
        for r in restarts:
            try:
                # restart.STEP.bin
                stem = Path(r).stem
                if '.' in stem:
                    s = int(stem.split('.')[-1])
                    steps.append(s)
            except:
                continue
        
        if steps:
            max_step = max(steps)
            print(f"  Updating steps to {max_step}")
            meta["steps"] = max_step
            
            # Also update active_seeds if there's a new lammps_S<seed>.log
            logs = glob.glob(str(dst / "lammps_S*.log"))
            for l in logs:
                try:
                    stem = Path(l).stem
                    if '_S' in stem:
                        seed = stem.split('_S')[-1]
                        if "active_seeds" not in meta: meta["active_seeds"] = []
                        if seed not in meta["active_seeds"]:
                            meta["active_seeds"].append(seed)
                except:
                    continue
            
            with open(meta_path, 'w') as f:
                json.dump(meta, f, indent=4)
    else:
        print(f"  Warning: No metadata found in {dst}")

    print(f"  Done merging {run_name}.")

if __name__ == "__main__":
    for run in RUNS:
        merge_run(run)
    
    print("\nAll runs merged. Please verify one manually before deleting the misplaced directories.")
