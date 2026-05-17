import sys
import os
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(ROOT_DIR, "scratch"))

import submit_flexible

jobs = [{
    "name": "Test_Polite_Job",
    "type": "fill",
    "walltime": "01:00:00",
    "ppn": 8,
    "mem": "4gb",
    "priority": -1024,
    "params": {
        "num_procs": 4,
        "num_threads": 1,
        "N": 4,
        "seed": 123456,
        "n_fill": 100,
        "n_hoppers": 1,
        "relax_steps": 1000,
        "simulation": "Test_Job",
        "dt": 1e-06,
        "viscosity": 0.000001
    }
}]

generated_files, _ = submit_flexible.submit_jobs(jobs, submit=False)
print(f"Generated file: {generated_files[0]}")

with open(generated_files[0], 'r') as f:
    content = f.read()
    print("\n--- PBS Script Content ---")
    print(content)
    if "#PBS -p -1024" in content:
        print("\n✅ Verification SUCCESS: Priority flag found!")
    else:
        print("\n❌ Verification FAILED: Priority flag NOT found!")
