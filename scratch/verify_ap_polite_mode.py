import sys
import os
import shutil
import json
from unittest.mock import MagicMock, patch

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)

# Path definitions
AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
AUTO_PILOT_BACKUP = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json.bak")

def setup_test_config():
    # 1. Backup original config if exists
    if os.path.exists(AUTO_PILOT_FILE):
        shutil.copyfile(AUTO_PILOT_FILE, AUTO_PILOT_BACKUP)
        print("Backed up auto_pilot.json successfully.")

    # 2. Write test configuration
    test_data = {
        "settings": {
            "max_concurrent": 4,
            "polite_mode": True,
            "cron_interval_minutes": 60,
            "enabled": True,
            "last_heartbeat": None
        },
        "goals": {
            "NewRoot_Polite": {
                "target_steps": 1000000,
                "increment": 100000,
                "in_place": True,
                "last_submitted": None,
                "status": "Idle",
                "mode": "fill",
                "polite_mode": True,
                "params": {
                    "seed": 111111,
                    "num_procs": 8,
                    "num_threads": 1,
                    "walltime": "24:00:00",
                    "ppn": 16,
                    "mem": "16gb"
                }
            },
            "NewRoot_Priority": {
                "target_steps": 1000000,
                "increment": 100000,
                "in_place": True,
                "last_submitted": None,
                "status": "Idle",
                "mode": "fill",
                "polite_mode": False,
                "params": {
                    "seed": 222222,
                    "num_procs": 8,
                    "num_threads": 1,
                    "walltime": "24:00:00",
                    "ppn": 16,
                    "mem": "16gb"
                }
            }
        }
    }

    with open(AUTO_PILOT_FILE, "w") as f:
        json.dump(test_data, f, indent=4)
    print("Created test auto_pilot.json config.")

def restore_backup():
    if os.path.exists(AUTO_PILOT_BACKUP):
        shutil.copyfile(AUTO_PILOT_BACKUP, AUTO_PILOT_FILE)
        os.remove(AUTO_PILOT_BACKUP)
        print("Restored original auto_pilot.json config.")
    elif os.path.exists(AUTO_PILOT_FILE):
        os.remove(AUTO_PILOT_FILE)

def run_verification():
    setup_test_config()
    
    try:
        from Pulse.auto_pilot_manager import run_manager
        
        # 1. Mock active student jobs
        mock_jobs = [
            {"id": "1001.master", "name": "StudentRun", "user": "student_user", "status": "R"},
            {"id": "1002.master", "name": "AnotherStudent", "user": "other_student", "status": "Q"}
        ]
        
        # We will intercept subprocess.run to catch the qsub call and read the generated PBS scripts
        submitted_jobs = []
        
        def mock_subprocess_run(args, **kwargs):
            # Check if it is calling qsub
            if args[0].endswith("qsub") or "qsub" in args[0]:
                script_path = args[1]
                with open(script_path, "r") as f:
                    content = f.read()
                
                # Extract job properties
                job_name = None
                priority = None
                for line in content.splitlines():
                    if "#PBS -N" in line:
                        job_name = line.split("-N")[-1].strip()
                    if "#PBS -p" in line:
                        priority = line.split("-p")[-1].strip()
                
                submitted_jobs.append({
                    "name": job_name,
                    "priority": priority,
                    "content": content
                })
                
                # Return dummy success
                mock_res = MagicMock()
                mock_res.stdout = f"Job_{len(submitted_jobs)}"
                return mock_res
            return MagicMock()

        # Patch PBS querying and submission inside run_manager
        with patch("Pulse.auto_pilot_manager.get_pbs_jobs", return_value=mock_jobs), \
             patch("Pulse.auto_pilot_manager.subprocess.run", side_effect=mock_subprocess_run), \
             patch("Pulse.auto_pilot_manager.fcntl.flock"): # Prevent lock file issues
            
            print("\nExecuting Auto-Pilot Manager under mock active student condition...")
            run_manager()

        # 2. Assertions
        print("\n--- ASSERTION VERIFICATION ---")
        
        # Verify that only the priority job was submitted
        if len(submitted_jobs) != 1:
            print(f"❌ FAILED: Expected exactly 1 job to be submitted, but got {len(submitted_jobs)}.")
            sys.exit(1)
            
        submitted_job = submitted_jobs[0]
        print(f"Submitted Job Name: {submitted_job['name']}")
        print(f"Submitted Job Priority: {submitted_job['priority']}")
        
        if "AP_NewRoot_Priority" not in submitted_job['name']:
            print(f"❌ FAILED: The wrong job was submitted! Expected 'AP_NewRoot_Priority', got '{submitted_job['name']}'.")
            sys.exit(1)
            
        if submitted_job['priority'] != "0":
            print(f"❌ FAILED: Expected priority 0 for priority job, got '{submitted_job['priority']}'.")
            sys.exit(1)
            
        print("✅ SUCCESS: Priority goal bypasses student active stand-down check and is successfully scheduled with normal priority (p = 0)!")
        
        # Verify polite job was yielded
        with open(AUTO_PILOT_FILE, "r") as f:
            updated_data = json.load(f)
            
        polite_goal = updated_data["goals"]["NewRoot_Polite"]
        priority_goal = updated_data["goals"]["NewRoot_Priority"]
        
        if polite_goal.get("status") == "Running" or polite_goal.get("last_submitted") is not None:
            print("❌ FAILED: Polite goal was not skipped / yielded correctly.")
            sys.exit(1)
            
        if priority_goal.get("status") != "Running" or priority_goal.get("last_submitted") is None:
            print("❌ FAILED: Priority goal was not updated to Running/Submitted status in config.")
            sys.exit(1)
            
        print("✅ SUCCESS: Config statuses are updated correctly!")
        
    finally:
        restore_backup()

if __name__ == "__main__":
    run_verification()
