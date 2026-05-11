import os
import sys
import shutil
import json
import re
from pathlib import Path

def merge_folders(target_dir: str, source_dir: str):
    """
    Merges source_dir (child) into target_dir (parent).
    """
    target = Path(target_dir)
    source = Path(source_dir)

    # Helper for extracting seeds from names
    def extract_seeds(name):
        return re.findall(r'_S(\d+)', name)

    if not target.exists() or not source.exists():
        print(f"Error: One of the paths does not exist: {target}, {source}")
        return False

    # 1. Load Metadata
    def load_meta(d: Path):
        for name in ["grid_metadata.json", "metadata.json"]:
            p = d / name
            if p.exists():
                with open(p, 'r') as f:
                    return json.load(f), p
        return None, None

    target_meta, target_meta_path = load_meta(target)
    source_meta, source_meta_path = load_meta(source)

    if not target_meta or not source_meta:
        print("Error: Could not find metadata in one of the folders.")
        return False

    print(f"Merging {source.name} into {target.name}...")

    # 2. Merge Dump Files (chain, bond, angle)
    for sub in ["chain", "bond", "angle"]:
        s_sub = source / sub
        t_sub = target / sub
        if not s_sub.exists(): continue
        t_sub.mkdir(parents=True, exist_ok=True)

        for f in s_sub.glob("*.dump"):
            dest = t_sub / f.name
            if dest.exists():
                # Compare sizes or just overwrite with the newer one from child
                if f.stat().st_size >= dest.stat().st_size:
                    shutil.move(str(f), str(dest))
                else:
                    f.unlink() # Child's version is smaller/broken? (Unlikely)
            else:
                shutil.move(str(f), str(dest))

    # 3. Merge Restarts
    s_restart = source / "restart"
    t_restart = target / "restart"
    if s_restart.exists():
        t_restart.mkdir(parents=True, exist_ok=True)
        for f in s_restart.glob("*.bin"):
            dest = t_restart / f.name
            # Keep restart.final.bin from child
            if f.name == "restart.final.bin" or not dest.exists():
                shutil.move(str(f), str(dest))
            else:
                # Compare step numbers if present in filename restart.100000.bin
                match_s = re.search(r'restart\.(\d+)\.bin', f.name)
                match_t = re.search(r'restart\.(\d+)\.bin', dest.name)
                if match_s and match_t:
                    if int(match_s.group(1)) > int(match_t.group(1)):
                        shutil.move(str(f), str(dest))
                    else:
                        f.unlink()
                else:
                    # Non-numbered restart, keep child's if newer
                    if f.stat().st_mtime > dest.stat().st_mtime:
                        shutil.move(str(f), str(dest))
                    else:
                        f.unlink()

    # 4. Replace split_states and .data files
    s_states = source / "split_states"
    t_states = target / "split_states"
    if s_states.exists():
        if t_states.exists():
            shutil.rmtree(t_states)
        shutil.move(str(s_states), str(t_states))
    
    # Move root .data files (e.g. final_grid.data)
    for f in source.glob("*.data"):
        dest = target / f.name
        if dest.exists():
            dest.unlink()
        shutil.move(str(f), str(dest))

    # Preserve Logs and Input Scripts with child seed
    seed_list = extract_seeds(source.name)
    suffix = f"_S{seed_list[-1]}" if seed_list else f"_{source.name}"
    
    for log_file in source.glob("lammps.log*"):
        new_name = f"lammps{suffix}.log"
        shutil.move(str(log_file), str(target / new_name))
    
    for in_file in source.glob("in.*"):
        new_name = f"{in_file.name}{suffix}"
        shutil.move(str(in_file), str(target / new_name))

    # 5. Consolidate Metadata
    # Update target steps by scanning dump files for the true max timestep
    max_ts = 0
    for f in (target / "chain").glob("*.dump"):
        match = re.search(r'chain_(\d+)\.dump', f.name)
        if match: max_ts = max(max_ts, int(match.group(1)))
    
    target_meta["steps"] = max_ts
    
    # Consolidate active_seeds
    t_seeds = target_meta.get("active_seeds", [])
    s_seeds = source_meta.get("active_seeds", [])
    
    all_seeds = set(t_seeds + s_seeds + extract_seeds(target.name) + extract_seeds(source.name))
    target_meta["active_seeds"] = sorted(list(all_seeds), key=lambda x: int(x) if x.isdigit() else 0)

    # Save target meta
    with open(target_meta_path, 'w') as f:
        json.dump(target_meta, f, indent=4)

    # 6. Update Lineage and Notes
    try:
        from Pulse.pulse_core import PBSManager
        lineage = PBSManager.load_lineage()
        source_id = str(source.absolute())
        target_id = str(target.absolute())

        # Load and update notes
        notes_path = Path("Pulse/lineage_notes.json")
        notes = {}
        if notes_path.exists():
            try:
                with open(notes_path, 'r') as f:
                    notes = json.load(f)
            except Exception:
                notes = {}

        if source_id in notes:
            child_note = notes.pop(source_id)
            if target_id in notes:
                notes[target_id] += f"\n\n--- Merged from {source.name} ---\n{child_note}"
            else:
                notes[target_id] = child_note
            
            with open(notes_path, 'w') as f:
                json.dump(notes, f, indent=4)
            print(f"Migrated research notes for {source.name} to {target.name}")

        if source_id in lineage and target_id in lineage:
            # Update target in lineage
            target_info = lineage[target_id]
            source_info = lineage[source_id]
            
            target_info["steps"] = max(target_info.get("steps", 0), source_info.get("steps", 0))
            target_info["status"] = "Active"    # Ensure it's not gray/archived
            target_info["sync_status"] = "Local" # Force re-sync
            
            # Record the merge in a history field for future tracking
            if "merge_history" not in target_info: target_info["merge_history"] = []
            target_info["merge_history"].append({
                "source_name": source.name,
                "source_seed": source_info.get("params", {}).get("seed"),
                "merged_at": target_info["steps"]
            })

            # Re-point children of source to target
            reparented_count = 0
            for rid, info in lineage.items():
                if info.get("parent") == source_id:
                    lineage[rid]["parent"] = target_id
                    reparented_count += 1
                    
                    # PHYSICAL REPAIR: Update the metadata file on disk for the survivor
                    for meta_name in ["grid_metadata.json", "metadata.json"]:
                        child_meta_path = Path(rid) / meta_name
                        if child_meta_path.exists():
                            try:
                                with open(child_meta_path, 'r') as f:
                                    child_meta = json.load(f)
                                
                                # Update parent path pointer
                                if "source_dir" in child_meta: child_meta["source_dir"] = target_id
                                if "restart_path" in child_meta: child_meta["restart_path"] = target_id
                                
                                with open(child_meta_path, 'w') as f:
                                    json.dump(child_meta, f, indent=4)
                                print(f"  -> Physically updated metadata for survivor: {info['name']}")
                            except Exception as e:
                                print(f"  -> Warning: Could not update child metadata file: {e}")
            
            # Remove source from lineage
            del lineage[source_id]
            
            PBSManager.save_lineage(lineage)
            print(f"Updated lineage.json: {source.name} merged into {target.name} ({reparented_count} descendants re-parented)")
    except Exception as e:
        print(f"Warning: Could not update lineage or notes: {e}")

    # 7. Cleanup Source
    try:
        shutil.rmtree(source)
        print(f"Successfully merged and deleted {source.name}")
    except Exception as e:
        print(f"Warning: Could not delete source directory {source}: {e}")

    return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_runs.py <target_parent_dir> <source_child_dir>")
        sys.exit(1)
    
    merge_folders(sys.argv[1], sys.argv[2])
