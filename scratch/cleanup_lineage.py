import json
import os
from pathlib import Path

LINEAGE_PATH = Path("/home/guest/palash/Granular-Chains-DEM/Pulse/lineage.json")

def cleanup():
    if not LINEAGE_PATH.exists():
        print("Lineage file not found.")
        return

    with open(LINEAGE_PATH, 'r') as f:
        lineage = json.load(f)

    # Group by name to identify duplicates
    by_name = {}
    for path, info in lineage.items():
        name = info.get("name")
        if name not in by_name:
            by_name[name] = []
        by_name[name].append((path, info))

    new_lineage = {}
    removed_count = 0

    for name, entries in by_name.items():
        if len(entries) == 1:
            # Only one entry, keep it
            new_lineage[entries[0][0]] = entries[0][1]
            continue
        
        # Multiple entries for the same name.
        # Find the best one (Active, or if multiple active, the one in the standard dumping yard)
        active_entries = [e for e in entries if e[1].get("status") == "Active"]
        
        if active_entries:
            # Keep only the active ones
            for path, info in active_entries:
                new_lineage[path] = info
            removed_count += (len(entries) - len(active_entries))
        else:
            # None are active, just keep the first one or keep all (to be safe)
            for path, info in entries:
                new_lineage[path] = info

    # Save the cleaned lineage
    with open(LINEAGE_PATH, 'w') as f:
        json.dump(new_lineage, f, indent=4)

    print(f"Cleanup complete. Removed {removed_count} archived duplicate entries from lineage.json.")

if __name__ == "__main__":
    cleanup()
