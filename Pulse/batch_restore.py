import subprocess
import os

# Configuration
RCLONE_PATH = "/home/guest/miniconda3/envs/gchain/bin/rclone"
CLOUD_BASE = "gdrive:Granular-Chains-DEM/dumping_yard_zips/"
LOCAL_TARGET_BASE = "/Data/palash_data/dumping_yard/Hopper_Fill/"
ZIP_TEMP_DIR = "/home/guest/palash/Granular-Chains-DEM/Zips_dir"

RUNS = [
    "Grid_Fill_1H_S468103",
    "Grid_Fill_1H_S935840",
    "Grid_Fill_1H_S437513",
    "Grid_Fill_1H_S481365",
    "Grid_Fill_1H_S240587",
    "Grid_Fill_1H_S904266",
    "Grid_Fill_1H_S358667",
    "Grid_Fill_1H_S808086",
    "Grid_Fill_1H_S442391",
    "Grid_Fill_1H_S469898",
    "Grid_Fill_1H_S374952",
    "Grid_Fill_1H_S694955"
]

def restore():
    os.makedirs(ZIP_TEMP_DIR, exist_ok=True)
    os.makedirs(LOCAL_TARGET_BASE, exist_ok=True)
    
    print(f"Starting batch restoration to {LOCAL_TARGET_BASE}...")
    
    for run_name in RUNS:
        zip_name = f"dumping_yard_Hopper_Fill_{run_name}.tar.gz"
        zip_path = os.path.join(ZIP_TEMP_DIR, zip_name)
        
        print(f"\n[1/3] Downloading {zip_name}...")
        try:
            subprocess.run([RCLONE_PATH, "copy", f"{CLOUD_BASE}{zip_name}", ZIP_TEMP_DIR], check=True)
            
            print(f"[2/3] Extracting to {LOCAL_TARGET_BASE}...")
            # Extracting one level up because the tar usually contains 'dumping_yard/Hopper_Fill/...' 
            # or just the folder. Let's check a tar structure.
            # Assuming it extracts the folder Grid_Fill_... directly.
            subprocess.run(["tar", "-xzf", zip_path, "-C", LOCAL_TARGET_BASE], check=True)
            
            print(f"[3/3] Cleaning up {zip_name}...")
            os.remove(zip_path)
            print(f"DONE: {run_name} is now local.")
            
        except Exception as e:
            print(f"ERROR: Failed to restore {run_name}: {e}")

if __name__ == "__main__":
    restore()
