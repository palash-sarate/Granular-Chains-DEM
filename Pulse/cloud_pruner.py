import json
import os
import subprocess
from pathlib import Path

# Configuration
RCLONE_PATH = "/home/guest/miniconda3/envs/gchain/bin/rclone"
LINEAGE_FILE = "Pulse/lineage.json"
CLOUD_PATH = "gdrive:Granular-Chains-DEM/dumping_yard_zips/"

def load_lineage_names():
    if not os.path.exists(LINEAGE_FILE):
        return set()
    with open(LINEAGE_FILE, "r") as f:
        lineage = json.load(f)
    
    # We collect all names currently in the lineage
    names = set()
    for info in lineage.values():
        names.add(info.get("name"))
    return names

def get_cloud_zips():
    print("Fetching cloud inventory from Google Drive...")
    try:
        result = subprocess.check_output([RCLONE_PATH, "lsf", CLOUD_PATH], text=True)
        return result.strip().split("\n")
    except Exception as e:
        print(f"Error fetching cloud inventory: {e}")
        return []

def prune():
    active_names = load_lineage_names()
    cloud_zips = get_cloud_zips()
    
    orphans = []
    
    print(f"\nAnalyzing {len(cloud_zips)} cloud files...")
    for zip_file in cloud_zips:
        if not zip_file.endswith(".tar.gz"):
            continue
            
        # Try to find the run name in the zip name
        # Convention: dumping_yard_Hopper_Fill_NAME.tar.gz
        # We strip the prefix and suffix
        name_part = zip_file.replace(".tar.gz", "")
        
        found_match = False
        for active_name in active_names:
            if active_name and active_name in name_part:
                found_match = True
                break
        
        if not found_match:
            orphans.append(zip_file)

    if not orphans:
        print("✅ No orphaned zips found. Your cloud storage is lean and clean!")
        return

    print(f"\n🚨 FOUND {len(orphans)} ORPHANED ZIPS ON GOOGLE DRIVE:")
    print("-" * 50)
    for i, orphan in enumerate(orphans):
        print(f"[{i+1}] {orphan}")
    print("-" * 50)
    print("These files exist in the cloud but are no longer in your Lineage History (likely due to merges).")
    
    confirm = input("\nDo you want to delete these orphans? (type 'DELETE' to confirm, anything else to cancel): ")
    
    if confirm == "DELETE":
        print("\nStarting safe deletion...")
        for orphan in orphans:
            print(f"Deleting {orphan}...")
            try:
                subprocess.run([RCLONE_PATH, "deletefile", f"{CLOUD_PATH}{orphan}"], check=True)
                print(f"Done.")
            except Exception as e:
                print(f"Error deleting {orphan}: {e}")
        print("\nCleanup complete.")
    else:
        print("\nDeletion cancelled. No files were removed.")

if __name__ == "__main__":
    prune()
