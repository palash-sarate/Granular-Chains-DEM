import streamlit as st
import os
import json
import glob
import time
import random
import re
import pandas as pd
from pathlib import Path
import streamlit.components.v1 as components
from Pulse.pulse_core import PBSManager, SyncManager
import submit_flexible
from lineage_tracker import scan_dumping_yard

# File paths
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NOTES_FILE = os.path.join(ROOT_DIR, "Pulse/lineage_notes.json")
AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse/auto_pilot.json")

def load_notes():
    if os.path.exists(NOTES_FILE):
        try:
            with open(NOTES_FILE, "r") as f: return json.load(f)
        except: pass
    return {}

def save_note(run_id, note):
    notes = load_notes()
    notes[run_id] = note
    with open(NOTES_FILE, "w") as f: json.dump(notes, f)

def load_auto_pilot():
    if os.path.exists(AUTO_PILOT_FILE):
        try:
            with open(AUTO_PILOT_FILE, "r") as f: return json.load(f)
        except: pass
    return {"settings": {"enabled": False}, "goals": {}}

def save_auto_pilot(data):
    with open(AUTO_PILOT_FILE, "w") as f:
        json.dump(data, f, indent=4)

def wrap_sim_name(name, max_width=25):
    """Wraps simulation names by breaking on underscores to keep nodes compact."""
    if len(name) <= max_width: return name
    parts = name.split("_")
    lines = []
    current_line = ""
    for part in parts:
        if len(current_line) + len(part) + 1 > max_width and current_line:
            lines.append(current_line)
            current_line = part
        else:
            current_line = (current_line + "_" + part) if current_line else part
    if current_line: lines.append(current_line)
    return "<br/>".join(lines)

def render_lineage(user_filter: str):
    st.subheader("🧬 Simulation Genealogy & Lineage")
    st.markdown("Persistent history of all simulation runs and their parent-child relationships.")
    
    c1, c2 = st.columns([1, 4])
    if c1.button("🔄 Sync from Disk", use_container_width=True):
        with st.spinner("Scanning dumping_yard..."):
            scan_dumping_yard()
            st.rerun()
    
    # Auto-initialize lineage once per session
    if "lineage_auto_scanned" not in st.session_state:
        with st.spinner("Scanning simulation lineage..."):
            try:
                scan_dumping_yard()
            except Exception:
                pass
            st.session_state["lineage_auto_scanned"] = True
    
    lineage = PBSManager.load_lineage()
    
    if not lineage:
        st.info("No lineage data found. Click 'Sync from Disk' to scan your simulations.")
    else:
        # Load Auto-Pilot goals for graph badges
        auto_data = load_auto_pilot()
        auto_goals = auto_data.get("goals", {})

        # Parse PBS Jobs for Status Overlay
        pbs_manager = PBSManager()
        try:
            active_jobs = pbs_manager.get_jobs()
        except Exception:
            active_jobs = []
            
        running_seeds = set()
        queued_ghosts = []
        node_to_jobid = {} # Map Mermaid node IDs (short_id or Ghost_id) to PBS Job IDs
        
        # Map of seed -> parent_short_id
        seed_to_shortid = {}
        for run_id, info in lineage.items():
            short_id = info['name'].replace('-', '_').replace('.', '_')
            seed = str(info.get('params', {}).get('seed', ''))
            if seed:
                seed_to_shortid[seed] = short_id
                
        for job in active_jobs:
            job_name = job.get('Job_Name', '')
            job_id = job.get('id', '')
            state = job.get('job_state', '')
            
            # Match convention: [Prefix][ParentSeed]_[ChildSeed]
            match = re.match(r"^[A-Z]+(\d{3,6})_(\d{3,6})$", job_name)
            if match:
                parent_seed, child_seed = match.groups()
                if state == 'R':
                    running_seeds.add(child_seed)
                    # Find corresponding short_id in lineage
                    found_in_lineage = False
                    for rid, info in lineage.items():
                        if str(info.get('params', {}).get('seed', '')) == child_seed:
                            s_id = info['name'].replace('-', '_').replace('.', '_')
                            node_to_jobid[s_id] = job_id
                            found_in_lineage = True
                            break
                    
                    # If running but not yet on disk, treat as a "Starting" ghost
                    if not found_in_lineage:
                        ghost_id = f"Ghost_{child_seed}"
                        ghost_label = f"{job_name}<br/>(Starting...)"
                        parent_short_id = seed_to_shortid.get(parent_seed)
                        queued_ghosts.append((ghost_id, ghost_label, parent_short_id))
                        node_to_jobid[ghost_id] = job_id

                elif state in ['Q', 'H', 'W', 'S']:
                    ghost_id = f"Ghost_{child_seed}"
                    ghost_label = f"{job_name}<br/>(Queued/Hold)"
                    
                    parent_short_id = None
                    if parent_seed in seed_to_shortid:
                        parent_short_id = seed_to_shortid[parent_seed]
                    
                    # Add ghost if it has a known parent OR if it's a root job (seed 000000)
                    if parent_short_id or parent_seed == "000000":
                        queued_ghosts.append((ghost_id, ghost_label, parent_short_id))
                        node_to_jobid[ghost_id] = job_id

        # 1. Prepare Mermaid Diagram
        mermaid_code = "graph LR\n"
        selected_parents = st.session_state.get("lineage_selected_parents", [])

        # Define nodes
        for run_id, info in lineage.items():
            short_id = info['name'].replace('-', '_').replace('.', '_')
            wrapped_name = wrap_sim_name(info['name'], max_width=20)
            node_label = f"{wrapped_name}<br/>(N={info['N']}, {info['steps']:,} steps)"
            
            # Append freq/amp for flow simulations
            p = info.get('params', {})
            f, a = p.get('freq'), p.get('amp')
            if f is not None or a is not None:
                # Handle single values or lists (for grid runs)
                f_val = f[0] if isinstance(f, list) else f
                a_val = a[0] if isinstance(a, list) else a
                if f_val or a_val:
                    node_label += f"<br/>f={f_val}, a={a_val}"
            
            # Append geometry_vars if present
            gv_data = p.get('geometry_vars', [])
            if gv_data:
                gv = gv_data[0] if isinstance(gv_data, list) and len(gv_data) > 0 else gv_data
                if isinstance(gv, dict) and gv:
                    label_map = {"orifice_half": "Orf_H", "hopper_ang": "Ang"}
                    gv_str = ", ".join([f"{label_map.get(k, k)}={v}" for k, v in gv.items()])
                    node_label += f"<br/>{gv_str}"

            # Append lepton_vars if present
            lv = p.get('lepton_vars', {})
            if lv and isinstance(lv, dict):
                label_map_l = {"ktheta": "kθ", "mu_s": "μ_s", "mu_roll": "μ_r", "theta_lim": "θ_lim"}
                lv_str = ", ".join([f"{label_map_l.get(k, k)}={v}" for k, v in lv.items()])
                node_label += f"<br/>Lepton: {lv_str}"

            seed = str(info.get('params', {}).get('seed', ''))
            
            # Base Style
            on_disk = os.path.exists(run_id)
            is_auto = False
            if run_id in auto_goals:
                goal = auto_goals[run_id]
                target = goal.get("target_steps", 0)
                current_steps = info.get("steps", 0)
                if current_steps < target:
                    is_auto = True
            
            if is_auto:
                node_label = f"🤖 {node_label}"

            style = "active"
            
            # Check for in-place running via active_seeds in metadata
            is_inplace_running = False
            if on_disk:
                try:
                    m_path = os.path.join(run_id, "grid_metadata.json")
                    if not os.path.exists(m_path):
                        m_path = os.path.join(run_id, "metadata.json")
                    if os.path.exists(m_path):
                        with open(m_path, 'r') as f:
                            m_data = json.load(f)
                            if "active_seeds" in m_data:
                                for s in m_data["active_seeds"]:
                                    if str(s) in running_seeds:
                                        is_inplace_running = True
                                        break
                except Exception:
                    pass

            if seed in running_seeds or is_inplace_running: 
                style = "running"
            elif info.get("sync_status") == "Synced":
                style = "synced" if on_disk else "cleared"
            elif info["status"] == "Archived": style = "archived"
            elif "Flow" in info["simulation"]: style = "flow"
            
            if is_auto: style = "auto_pilot"

            mermaid_code += f'    {short_id}["{node_label}"]:::{style}\n'
            mermaid_code += f'    click {short_id} call selectNode("{short_id}")\n'
            
            # Apply selection highlight separately
            if run_id in selected_parents:
                mermaid_code += f"    class {short_id} selected\n"

        # Define relationships
        for run_id, info in lineage.items():
            if info["parent"] and info["parent"] in lineage:
                parent_name = lineage[info["parent"]]['name'].replace('-', '_').replace('.', '_')
                child_name = info['name'].replace('-', '_').replace('.', '_')
                mermaid_code += f"    {parent_name} --> {child_name}\n"
                
        # Inject Queued Ghost Nodes
        for ghost_id, ghost_label, parent_short_id in queued_ghosts:
            mermaid_code += f'    {ghost_id}["{ghost_label}"]:::queued\n'
            mermaid_code += f'    click {ghost_id} call selectNode("{ghost_id}")\n'
            if ghost_id in selected_parents:
                mermaid_code += f"    class {ghost_id} selected\n"
            if parent_short_id:
                mermaid_code += f"    {parent_short_id} --> {ghost_id}\n"

        # Add "+" Node for new fill runs
        mermaid_code += '    NewRoot[" + Start New Fill Run "]:::new_node\n'
        mermaid_code += '    click NewRoot call selectNode("NewRoot")\n'
        if "NewRoot" in selected_parents:
            mermaid_code += "    class NewRoot selected\n"

        # Define styles
        mermaid_code += "    classDef active fill:#e1f5fe,stroke:#01579b,stroke-width:2px;\n"
        mermaid_code += "    classDef archived fill:#f5f5f5,stroke:#9e9e9e,stroke-dasharray: 5 5;\n"
        mermaid_code += "    classDef synced fill:#80cbc4,stroke:#00695c,stroke-width:2px;\n"
        mermaid_code += "    classDef cleared fill:#b39ddb,stroke:#512da8,stroke-width:2px;\n"
        mermaid_code += "    classDef flow fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;\n"
        mermaid_code += "    classDef running fill:#a5d6a7,stroke:#2e7d32,stroke-width:3px;\n"
        mermaid_code += "    classDef queued fill:#ffccbc,stroke:#d32f2f,stroke-dasharray: 5 5;\n"
        mermaid_code += "    classDef new_node fill:#ffffff,stroke:#333333,stroke-width:2px,stroke-dasharray: 5 5;\n"
        mermaid_code += "    classDef auto_pilot fill:#b2ebf2,stroke:#00acc1,stroke-width:3px;\n"
        mermaid_code += "    classDef selected stroke:#ff9800,stroke-width:4px;\n"

        # Render Mermaid using custom component
        mermaid_click_component = components.declare_component(
            "mermaid_click",
            path=os.path.join(ROOT_DIR, "Pulse/mermaid_component")
        )
        
        clicked_node = mermaid_click_component(mermaid_code=mermaid_code, default=None)
        
        # Handle click from the custom component natively
        if clicked_node:
            if clicked_node == "NewRoot":
                if "lineage_selected_parents" not in st.session_state:
                    st.session_state["lineage_selected_parents"] = []
                
                if "NewRoot" in st.session_state["lineage_selected_parents"]:
                    st.session_state["lineage_selected_parents"].remove("NewRoot")
                else:
                    st.session_state["lineage_selected_parents"].append("NewRoot")
                st.rerun()

            # 1. Job Deletion Interface (for ongoing/on-hold jobs)
            if clicked_node in node_to_jobid:
                active_jid = node_to_jobid[clicked_node]
                st.divider()
                st.warning(f"⚠️ **Ongoing Job Selected:** `{clicked_node}` (Job ID: `{active_jid}`)")
                
                if st.session_state.get("confirm_delete_id") == active_jid:
                    st.error("Are you absolutely sure you want to terminate this job? This will stop the simulation immediately.")
                    c1, c2 = st.columns([1, 4])
                    if c1.button("🔥 YES, TERMINATE", type="primary", use_container_width=True):
                        try:
                            PBSManager.delete_job(active_jid)
                            st.toast(f"Successfully sent qdel for {active_jid}")
                            del st.session_state["confirm_delete_id"]
                            time.sleep(1.5)
                            st.rerun()
                        except Exception as e:
                            st.error(f"Failed to delete job: {e}")
                    if c2.button("Cancel", use_container_width=True):
                        del st.session_state["confirm_delete_id"]
                        st.rerun()
                else:
                    if st.button(f"🛑 Terminate Ongoing Job ({active_jid})", use_container_width=True, type="secondary"):
                        st.session_state["confirm_delete_id"] = active_jid
                        st.rerun()
            
            # 2. Archived/Metadata-only Deletion
            id_to_path = {info['name'].replace('-', '_').replace('.', '_'): rid for rid, info in lineage.items()}
            if clicked_node in id_to_path:
                rid = id_to_path[clicked_node]
                info = lineage[rid]
                
                if info.get("status") == "Archived":
                    st.divider()
                    st.warning(f"📦 **Archived Run Selected:** `{info['name']}`")
                    
                    if st.session_state.get("confirm_archive_delete") == rid:
                        st.error(f"Remove `{info['name']}` from lineage history? (Files on disk will NOT be deleted)")
                        c1, c2 = st.columns([1, 4])
                        if c1.button("🗑️ REMOVE FROM LINEAGE", type="primary", use_container_width=True):
                            del lineage[rid]
                            PBSManager.save_lineage(lineage)
                            st.toast(f"Removed {info['name']} from lineage.")
                            del st.session_state["confirm_archive_delete"]
                            st.rerun()
                        if c2.button("Cancel", use_container_width=True):
                            del st.session_state["confirm_archive_delete"]
                            st.rerun()
                    else:
                        if st.button(f"🗑️ Remove `{info['name']}` from Lineage Metadata", use_container_width=True):
                            st.session_state["confirm_archive_delete"] = rid
                            st.rerun()

            if clicked_node in id_to_path:
                rid = id_to_path[clicked_node]
                info = lineage[rid]
                
                # 3.1 Display Sync & Restore Options
                sync_status = info.get("sync_status", "Local")
                is_local = os.path.exists(rid)
                
                st.divider()
                st.markdown(f"### 📦 Storage Status: **{sync_status}**")
                
                # Display Run Parameters
                p_info = info.get('params', {})
                if p_info:
                    st.markdown("#### 📋 Run Parameters")
                    main_params = {k: v for k, v in p_info.items() if k not in ["geometry_vars", "lepton_vars"] and v is not None}
                    if main_params:
                        p_cols = st.columns(min(len(main_params), 4))
                        for i, (k, v) in enumerate(main_params.items()):
                            p_cols[i % 4].metric(k.capitalize(), str(v))
                    
                    geo_vars = p_info.get("geometry_vars")
                    if geo_vars and geo_vars != {}:
                        with st.expander("🛠️ Geometry Overrides", expanded=True):
                            st.json(geo_vars)

                    lepton_vars = p_info.get("lepton_vars")
                    if lepton_vars and lepton_vars != {}:
                        with st.expander("🔬 Lepton parameter Overrides", expanded=True):
                            st.json(lepton_vars)
                
                c1, c2 = st.columns([1, 1])
                if sync_status == "Synced" and not is_local:
                    if c1.button(f"📥 Restore {info['name']} to Disk", use_container_width=True):
                        with st.spinner("Downloading from Google Drive..."):
                            success, msg = SyncManager.restore_run(rid)
                            if success:
                                st.success(msg)
                                time.sleep(1)
                                st.rerun()
                            else:
                                st.error(msg)
                elif is_local:
                    c1.success("✅ Files are available offline.")
                
                st.divider()
                st.markdown("#### 🤖 Auto-Pilot Control")
                auto_data = load_auto_pilot()
                if rid in auto_data["goals"]:
                    g = auto_data["goals"][rid]
                    st.success(f"Target: {g['target_steps']} | Increment: {g['increment']}")
                    if st.button("❌ Remove from Auto-Pilot", key=f"rm_auto_{rid}"):
                        del auto_data["goals"][rid]
                        save_auto_pilot(auto_data)
                        st.rerun()
                else:
                    st.info("💡 Use the **Launch** section below to configure standard or Auto-Pilot runs for this node.")

                st.divider()

                # 4. Regular Parent Selection
                new_parent = id_to_path[clicked_node]
                if "lineage_selected_parents" not in st.session_state:
                    st.session_state["lineage_selected_parents"] = []
                
                if new_parent in st.session_state["lineage_selected_parents"]:
                    st.session_state["lineage_selected_parents"].remove(new_parent)
                else:
                    st.session_state["lineage_selected_parents"].append(new_parent)
                st.rerun()

        # 2. Detailed Data View
        with st.expander("📄 View Detailed Parameters"):
            display_lineage = []
            for rid, info in lineage.items():
                n_val = info.get("N", 0)
                if isinstance(n_val, list) and len(n_val) > 0: n_val = n_val[0]

                display_lineage.append({
                    "Name": info["name"],
                    "Type": info["simulation"],
                    "N": n_val,
                    "Steps": info["steps"],
                    "Status": info["status"],
                    "GeoVars": str(info.get("params", {}).get("geometry_vars", {})),
                    "LeptonVars": str(info.get("params", {}).get("lepton_vars", {})),
                    "Path": rid
                })
            st.dataframe(pd.DataFrame(display_lineage), width="stretch", hide_index=True)
            
        st.divider()
        st.subheader("🚀 Launch New Simulation(s) from Lineage")
        
        # Select Parent Runs
        run_options = ["NewRoot"] + list(lineage.keys())
        def format_run(rid):
            if rid == "NewRoot": return "🆕 Start New Fill Run (+)"
            return f"{lineage[rid]['name']} ({lineage[rid]['simulation']})"
        
        if "lineage_selected_parents" not in st.session_state:
            st.session_state["lineage_selected_parents"] = []
        
        selected_parents = st.multiselect(
            "Select Parent Run(s)", 
            run_options, 
            default=st.session_state["lineage_selected_parents"],
            format_func=lambda x: format_run(x)
        )
        
        st.session_state["lineage_selected_parents"] = selected_parents
        
        if selected_parents:
            # Validation and Categorization
            if "NewRoot" in selected_parents and len(selected_parents) > 1:
                st.error("⚠️ **Conflict:** 'New Root' cannot be combined with existing parent runs. Please select one or the other.")
                st.stop()

            if selected_parents == ["NewRoot"]:
                selected_mode = "fill"
                st.info("🆕 **Creating New Root Fill Run** (No Parent)")
                active_parent_jids = {}
                representative_parent = "NewRoot"
            else:
                parent_infos = {rid: lineage[rid] for rid in selected_parents}
                parent_types = set(info.get("simulation", "") for info in parent_infos.values())
                
                if len(parent_types) > 1:
                    st.error(f"⚠️ **Incompatible Types:** You have selected runs of different categories ({parent_types}). All selected runs must be either all 'Fill' or all 'Flow'.")
                    st.stop()
                
                parent_type = list(parent_types)[0]
                if "Flow" in parent_type:
                    allowed_modes = ["flow_resume"]
                else:
                    allowed_modes = ["fill_resume", "flow", "calibration"]
                    
                selected_mode = st.radio("Select Simulation Mode", allowed_modes, horizontal=True)
                representative_parent = selected_parents[0]
                
                # Resolve dependencies
                active_parent_jids = {}
                for rid, info in parent_infos.items():
                    short_id = info['name'].replace('-', '_').replace('.', '_')
                    jid = node_to_jobid.get(short_id)
                    if jid:
                        active_parent_jids[rid] = jid
                
                if active_parent_jids:
                    st.info(f"🔗 **Lineage Dependencies:** {len(active_parent_jids)} of the selected runs are still active. Dependent jobs will wait for them.")
                
                # Consolidation / Merge
                if len(selected_parents) >= 2:
                    def get_depth(node_id, current_lineage):
                        depth = 0
                        curr = node_id
                        while curr and curr in current_lineage and current_lineage[curr].get("parent"):
                            curr = current_lineage[curr]["parent"]
                            depth += 1
                        return depth
                    
                    sorted_nodes = sorted(selected_parents, key=lambda x: get_depth(x, lineage))
                    ultimate_parent = sorted_nodes[0]
                    children_to_merge = sorted_nodes[1:]
                    
                    def is_descendant(child_id, ancestor_id, current_lineage):
                        curr = child_id
                        while curr and curr in current_lineage:
                            p = current_lineage[curr].get("parent")
                            if p == ancestor_id: return True
                            curr = p
                        return False

                    all_valid = all(is_descendant(c, ultimate_parent, lineage) for c in children_to_merge)
                    
                    if all_valid:
                        st.sidebar.markdown("---")
                        st.sidebar.subheader("🛠️ Batch Consolidation")
                        st.sidebar.info(f"Merge **{len(children_to_merge)}** runs into **{lineage[ultimate_parent]['name']}**.")
                        st.sidebar.warning("⚠️ **Storage Tip**: Shortening the lineage to < 3 nodes may cause the Sync Manager to keep this data local rather than archiving it to the cloud, as it considers the branch 'actively being worked on'.")
                        
                        if st.sidebar.button(f"🚀 Merge {len(children_to_merge)} Runs", use_container_width=True):
                            try:
                                from Pulse.merge_runs import merge_folders
                                with st.sidebar.status("Batch consolidating...", expanded=True) as status:
                                    for child_id in children_to_merge:
                                        child_name = lineage[child_id]['name']
                                        status.write(f"Processing {child_name}...")
                                        
                                        if not os.path.exists(ultimate_parent):
                                            status.write(f"Restoring parent: {lineage[ultimate_parent]['name']}...")
                                            SyncManager.restore_run(ultimate_parent)
                                        if not os.path.exists(child_id):
                                            status.write(f"Restoring child: {child_name}...")
                                            SyncManager.restore_run(child_id)
                                        
                                        success = merge_folders(ultimate_parent, child_id)
                                        if not success: raise Exception(f"Merge failed for {child_name}")
                                        
                                    status.update(label="✅ Batch Consolidation Complete!", state="complete")
                                    st.session_state["lineage_selected_parents"] = [ultimate_parent]
                                    time.sleep(1)
                                    st.rerun()
                            except Exception as e:
                                st.sidebar.error(f"Batch failed: {e}")
                    else:
                        st.sidebar.warning("⚠️ **Selection Mismatch:** Selected nodes must belong to the same linear branch for batch merge.")
            
            if selected_parents != ["NewRoot"]:
                primary_parent = selected_parents[0]
                parent_info = lineage[primary_parent]
                
                st.markdown(f"### 🔍 Inspecting: `{parent_info['name']}`")
                note_col, snap_col = st.columns([1, 1])
                
                with note_col:
                    notes_db = load_notes()
                    current_note = notes_db.get(primary_parent, "")
                    new_note = st.text_area("🗒️ Run Notes", value=current_note, height=100, help="Save observations or metadata for this run.")
                    if st.button("💾 Save Notes", key=f"save_note_{primary_parent}"):
                        save_note(primary_parent, new_note)
                        st.toast("Notes saved!")
                
                with snap_col:
                    with st.expander("🖼️ View Simulation Preview"):
                        with st.spinner("Generating snapshot..."):
                            try:
                                from Pulse.snapshot_helper import generate_snapshot
                                force_refresh = st.button("🔄 Refresh Snapshot", key=f"refresh_snap_{primary_parent}")
                                snap_path = generate_snapshot(primary_parent, force=force_refresh)
                                
                                if snap_path:
                                    st.image(snap_path, caption=f"Last Snapshot of {parent_info['name']} (Y-Z Plane)", use_container_width=True)
                                else:
                                    st.info("No snapshot available (no dump files found).")
                            except Exception as e:
                                st.warning(f"Snapshot preview unavailable: {e}")

                st.divider()
            
            # REACTIVE LAUNCH PARAMETERS
            st.markdown("### 🌍 Global Parameters")
            gc1, gc2, gc3, gc4 = st.columns(4)
            walltime = gc1.text_input("Walltime", value="24:00:00")
            ppn = gc2.number_input("PPN", value=16, step=1)
            mem = gc3.text_input("Memory", value="16gb")
            max_concurrent = gc4.number_input("Max Concurrent Jobs", value=4, min_value=1, step=1)
            
            polite_mode = st.checkbox("😇 Polite Mode (Low Priority)", value=True)
            
            gc5, gc6, gc7, gc8 = st.columns(4)
            num_procs = gc5.number_input("Num Procs", value=8, step=1)
            num_threads = gc6.number_input("Num Threads", value=1, step=1)
            dt = gc7.number_input("dt", value=1e-06, format="%e")
            viscosity = gc8.number_input("Viscosity", value=0.000001, format="%f")
            
            seed = st.number_input("Seed", value=random.randint(100000, 999999), step=1)
            
            st.markdown(f"### ⚙️ Mode-Specific Parameters ({selected_mode})")
            mode_params = {}
            
            if selected_mode == "fill":
                fc1, fc2, fc3 = st.columns(3)
                mode_params["N"] = fc1.number_input("Chain Length (N)", value=4, step=1)
                mode_params["n_fill"] = fc2.number_input("Number of Chains", value=3600, step=100)
                mode_params["relax_steps"] = fc3.number_input("Relax Steps", value=1000000, step=100000)
                
                fc4, fc5, fc6 = st.columns(3)
                mode_params["source_dir"] = fc4.text_input("Source Relaxed Chains", value="chain_data/relaxed_2D_x")
                mode_params["spacing"] = fc5.number_input("Hopper Spacing", value=0.5)
                mode_params["n_hoppers"] = fc6.number_input("Number of Hoppers", value=1, step=1)

                fc7, fc8, fc9 = st.columns(3)
                mode_params["mode"] = fc7.selectbox("Pouring Mode", ["2D_stacked", "2D_worst_case"], index=0)
                mode_params["dump_file"] = fc8.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                mode_params["hopper_template_data"] = fc9.text_input("Hopper Template", value="simulation_geometries/2D_hopper_with_orifice_cover.inc")
                mode_params["simulation"] = "Hopper_Fill"
                mode_params["no-vtk"] = True

                with st.expander("🛠️ Geometry Overrides", expanded=False):
                    st.info("Provide custom geometry parameters for the hopper(s) as a JSON object.")
                    geo_json = st.text_area("Geometry Variables (JSON)", value="{}", help='Example: {"orifice_half": 0.04, "hopper_ang": 55.0}')
                    try:
                        if geo_json.strip():
                            parsed_geo = json.loads(geo_json)
                            if not isinstance(parsed_geo, dict):
                                st.error("JSON must be a dictionary/object.")
                                mode_params["_invalid_geo"] = True
                            else:
                                mode_params["geometry_vars"] = parsed_geo
                                if parsed_geo != {}:
                                    st.success("Configuration loaded.")
                        else:
                            mode_params["geometry_vars"] = {}
                    except Exception as e:
                        st.error(f"Invalid JSON: {e}")
                        mode_params["_invalid_geo"] = True

            elif selected_mode == "fill_resume":
                rc1, rc2 = st.columns(2)
                mode_params["relax_steps"] = rc1.number_input("Relax Steps", value=1000000, step=100000)
                mode_params["dump_file"] = rc2.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                mode_params["simulation"] = "Hopper_Fill_Resume"
                
            elif selected_mode == "flow":
                ffc1, ffc2, ffc3, ffc4 = st.columns(4)
                freq_input = ffc1.text_input("Frequency(s)", value="5.0")
                amp_input = ffc2.text_input("Amplitude(s)", value="0.001")
                mode_params["run_steps"] = ffc3.number_input("Run Steps", value=2000000, step=100000)
                mode_params["osc_dir"] = ffc4.text_input("Oscillation Dir", value="z")
                mode_params["_freq_list"] = [f.strip() for f in freq_input.split(",") if f.strip()]
                mode_params["_amp_list"] = [a.strip() for a in amp_input.split(",") if a.strip()]
                
            elif selected_mode == "flow_resume":
                frc1, frc2 = st.columns(2)
                mode_params["run_steps"] = frc1.number_input("Run Steps", value=1000000, step=100000)
                
            elif selected_mode == "calibration":
                cc1, cc2, cc3 = st.columns(3)
                ktheta = cc1.number_input("Angular Stiffness (ktheta)", value=10.0, step=10.0)
                mu_s = cc2.number_input("Sliding Friction (mu_s)", value=0.3, step=0.1)
                mu_roll = cc3.number_input("Rolling Friction (mu_roll)", value=0.3, step=0.1)
                
                cc4, cc5 = st.columns(2)
                mode_params["relax_steps"] = cc4.number_input("Relaxation Steps", value=500000, step=100000)
                mode_params["dump_file"] = cc5.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                mode_params["simulation"] = "Hopper_Fill_Resume"
                
                l_vars = {"ktheta": ktheta, "mu_s": mu_s, "mu_roll": mu_roll}
                
                with st.expander("🔬 Other Lepton parameter Overrides (JSON)", expanded=False):
                    st.info("Provide any other custom lepton variables as a JSON object.")
                    lepton_json = st.text_area("Other Lepton Variables (JSON)", value="{}", help='Example: {"k_n": 1e6, "nu_n0": 50.0}')
                    try:
                        if lepton_json.strip():
                            parsed_lepton = json.loads(lepton_json)
                            if not isinstance(parsed_lepton, dict):
                                st.error("JSON must be a dictionary/object.")
                                mode_params["_invalid_lepton"] = True
                            else:
                                l_vars.update(parsed_lepton)
                    except Exception as e:
                        st.error(f"Invalid JSON: {e}")
                        mode_params["_invalid_lepton"] = True
                
                mode_params["lepton_vars"] = l_vars

            st.markdown("---")
            use_auto_pilot = st.checkbox("🤖 **Handover to Auto-Pilot**", value=False, help="Automatically queue and resume simulations until a target step count is reached.")
            
            ap_target, ap_inc = 10000000, 100000
            if use_auto_pilot:
                ap_col1, ap_col2 = st.columns(2)
                ap_target = ap_col1.number_input("Target Total Steps", value=10000000, step=1000000)
                ap_inc = ap_col2.number_input("Step Increment", value=100000, step=10000)
                inplace_mode = True
            else:
                inplace_mode = False
                if selected_mode != "fill":
                    inplace_mode = st.checkbox("😇 In-place Resumption", value=True, help="Run simulation in the same directory as parent.")
            
            submit_btn = st.button("🚀 Launch Auto-Pilot" if use_auto_pilot else "🚀 Submit to PBS", use_container_width=True)
                
            if submit_btn:
                if mode_params.get("_invalid_geo"):
                    st.error("Submission blocked: Please fix the invalid JSON in Geometry Overrides.")
                    st.stop()
                if mode_params.get("_invalid_lepton"):
                    st.error("Submission blocked: Please fix the invalid JSON in Lepton Overrides.")
                    st.stop()
                jobs_to_submit = []
                
                for p_rid in selected_parents:
                    if p_rid == "NewRoot":
                        p_info = {"name": "S000000", "params": {"seed": "000000"}}
                        p_seed = "000000"
                    else:
                        p_info = lineage[p_rid]
                        p_seed = str(p_info.get("params", {}).get("seed", "000000"))[-6:]
                    
                    p_active_jid = active_parent_jids.get(p_rid)
                    
                    # Handle Parametric Sweep for Flow
                    if selected_mode == "flow" and len(mode_params.get("_freq_list", [])) * len(mode_params.get("_amp_list", [])) > 1:
                        for f in mode_params["_freq_list"]:
                            for a in mode_params["_amp_list"]:
                                child_seed = random.randint(100000, 999999)
                                job_name = f"F{p_seed}_{str(child_seed)[-6:]}"
                                
                                job_params = {
                                    "num_procs": int(num_procs),
                                    "num_threads": int(num_threads),
                                    "dt": dt,
                                    "viscosity": viscosity,
                                    "seed": child_seed,
                                    "source_dir": p_rid,
                                    **mode_params
                                }

                                job_params["freq"] = float(f)
                                job_params["amp"] = float(a)
                                job_params["simulation"] = f"Flow_Study_N{p_info.get('N', 0)}_F{f}_A{a}"
                                job_params.pop("_freq_list", None)
                                job_params.pop("_amp_list", None)
                                
                                jobs_to_submit.append({
                                    "name": job_name,
                                    "type": selected_mode,
                                    "walltime": walltime,
                                    "ppn": int(ppn),
                                    "mem": mem,
                                    "priority": -1024 if polite_mode else 0,
                                    "dependency": p_active_jid,
                                    "params": {**job_params, "inplace": inplace_mode}
                                })
                    else:
                        if len(selected_parents) > 1:
                            current_job_seed = random.randint(100000, 999999)
                        else:
                            current_job_seed = int(seed)
                            
                        new_seed_suffix = str(current_job_seed)[-6:]
                        if selected_mode == "fill_resume": prefix = "R"
                        elif selected_mode == "flow": prefix = "F"
                        elif selected_mode == "flow_resume": prefix = "FR"
                        elif selected_mode == "calibration": prefix = "CAL"
                        else: prefix = "S"
                        
                        job_name = f"{prefix}{p_seed}_{new_seed_suffix}"
                        
                        final_params = {
                            "num_procs": int(num_procs),
                            "num_threads": int(num_threads),
                            "dt": dt,
                            "viscosity": viscosity,
                            "seed": current_job_seed,
                            **mode_params
                        }
                        
                        if selected_mode in ["fill_resume", "calibration"]:
                            final_params["restart_path"] = p_rid
                        elif selected_mode == "flow":
                            final_params["source_dir"] = p_rid
                            final_params["freq"] = float(mode_params.get("_freq_list", [0])[0])
                            final_params["amp"] = float(mode_params.get("_amp_list", [0])[0])
                        elif selected_mode == "flow_resume":
                            final_params["restart_path"] = p_rid
                        
                        final_params.pop("_freq_list", None)
                        final_params.pop("_amp_list", None)

                        job_data_params = {**final_params}
                        if selected_mode != "fill":
                            job_data_params["inplace"] = inplace_mode

                        jobs_to_submit.append({
                            "name": job_name,
                            "type": selected_mode,
                            "walltime": walltime,
                            "ppn": int(ppn),
                            "mem": mem,
                            "priority": -1024 if polite_mode else 0,
                            "dependency": p_active_jid,
                            "params": job_data_params
                        })
                
                if use_auto_pilot:
                    with st.spinner("Enrolling goals..."):
                        auto_data = load_auto_pilot()
                        enrolled_count = 0
                        for p_rid in selected_parents:
                            path_key = p_rid
                            if p_rid == "NewRoot":
                                path_key = f"NewRoot_{selected_mode}_{int(time.time())}"
                            
                            auto_data["goals"][path_key] = {
                                "target_steps": ap_target,
                                "increment": ap_inc,
                                "in_place": inplace_mode,
                                "last_submitted": None,
                                "status": "Idle",
                                "mode": selected_mode,
                                "polite_mode": polite_mode,
                                "paused": False,
                                "params": {
                                    "dt": dt,
                                    "viscosity": viscosity,
                                    "num_procs": int(num_procs),
                                    "num_threads": int(num_threads),
                                    "walltime": walltime,
                                    "ppn": int(ppn),
                                    "mem": mem,
                                    **{k:v for k,v in mode_params.items() if not k.startswith("_") and k not in ["relax_steps", "run_steps"]}
                                }
                            }
                            enrolled_count += 1
                        
                        save_auto_pilot(auto_data)
                        st.success(f"🤖 Successfully enrolled {enrolled_count} run(s) into Auto-Pilot!")
                        time.sleep(2)
                        st.rerun()
                else:
                    try:
                        importlib.reload(submit_flexible)
                        with st.spinner(f"Submitting {len(jobs_to_submit)} job(s)..."):
                            generated, submitted = submit_flexible.submit_jobs(jobs_to_submit, submit=True, max_concurrent=int(max_concurrent), user=user_filter)
                        if submitted:
                            st.success(f"Successfully submitted {len(submitted)} job(s)!")
                            if len(submitted) > 1:
                                st.info(f"First Job ID: {submitted[0]} | Last Job ID: {submitted[-1]}")
                        else:
                            st.warning("Jobs were generated but not submitted or submission failed.")
                    except Exception as e:
                        import traceback
                        with open("submit_error.log", "a") as errf:
                            errf.write(f"\n--- Error at {time.strftime('%Y-%m-%d %H:%M:%S')} ---\n")
                            errf.write(traceback.format_exc())
                        st.error(f"Error submitting job: {e}. Check submit_error.log for details.")
