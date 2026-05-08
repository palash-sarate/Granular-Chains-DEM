import os
import sys
import argparse

# Add project root to sys.path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pulse_core import SyncManager

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Pulse Sync Cycle as a PBS Job")
    parser.add_argument("--user", default="guest", help="User filter for PBS check")
    args = parser.parse_args()
    
    print(f"--- HPC Sync Job Started for user: {args.user} ---")
    SyncManager.run_sync_cycle(user=args.user)
    print("--- HPC Sync Job Finished ---")
