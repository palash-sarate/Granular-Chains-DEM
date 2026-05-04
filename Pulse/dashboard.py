import streamlit as st
import pandas as pd
import sys
import os
import json

# Prevent "Errno 5" I/O errors by redirecting stdout/stderr to devnull
# This is required for libraries that try to flush standard streams in headless environments
try:
    if sys.stdout is not None:
        sys.stdout.flush()
except OSError:
    sys.stdout = open(os.devnull, 'w')
try:
    if sys.stderr is not None:
        sys.stderr.flush()
except OSError:
    sys.stderr = open(os.devnull, 'w')

import time
import subprocess
import psutil
import tempfile

# Add project root to sys.path to allow importing from analysis module
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pulse_core import PBSManager, SimulationMonitor
import streamlit.components.v1 as components
from analysis.controllers.sim_data import SimDataController
from analysis.controllers.renderer import SimulationRenderer
from analysis.controllers.highlighter import HighlightController
import socket
import random

scratch_dir = os.path.join(ROOT_DIR, "scratch")
if scratch_dir not in sys.path:
    sys.path.append(scratch_dir)
import submit_flexible
from lineage_tracker import scan_dumping_yard


# Persistent Notes Storage
NOTES_FILE = os.path.join(ROOT_DIR, "Pulse/lineage_notes.json")
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

def st_directory_picker(label, key, base_path):
    """A simple directory picker for Streamlit."""
    if key not in st.session_state:
        st.session_state[key] = base_path
        
    curr = st.session_state[key]
    
    # Ensure path exists
    if not os.path.exists(curr):
        curr = base_path
        st.session_state[key] = base_path

    st.markdown(f"**{label}**")
    st.code(curr, language="bash")
    
    c1, c2, c3 = st.columns([1, 1, 3])
    if c1.button("⬆️ Up", key=f"{key}_up"):
        st.session_state[key] = os.path.dirname(curr)
        st.rerun()
    if c2.button("🏠 Home", key=f"{key}_home"):
        st.session_state[key] = base_path
        st.rerun()
        
    try:
        subdirs = sorted([d for d in os.listdir(curr) if os.path.isdir(os.path.join(curr, d)) and not d.startswith(".")])
        if subdirs:
            chosen = st.selectbox("Browse subdirectories:", ["-- Select to enter --"] + subdirs, key=f"{key}_browse")
            if chosen != "-- Select to enter --":
                st.session_state[key] = os.path.join(curr, chosen)
                st.rerun()
    except Exception as e:
        st.error(f"Access error: {e}")
        
    return st.session_state[key]

if st.runtime.exists():
    st.set_page_config(page_title="Pulse Dashboard", page_icon="⚡", layout="wide")

    st.title("⚡ Pulse Job Monitor")
    st.markdown("Real-time PBS Cluster Dashboard")

    # Sidebar for controls
    st.sidebar.header("Controls")
    refresh_rate = st.sidebar.slider("Refresh Rate (seconds)", 5, 60, 10)
    user_filter = st.sidebar.text_input("User Filter", value="guest")

    if st.sidebar.button("Refresh Now", use_container_width=False):
        st.rerun()

    # Create Tabs for Active vs History vs ETA vs Visualizer vs Lineage
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["📊 Active Queue", "🕰️ Job History", "⏱️ Simulation ETA", "🎥 Visualizer", "🧬 Lineage"])

    with tab1:
        @st.fragment(run_every=refresh_rate)
        def render_active_queue():
            # Fetch active jobs
            jobs = PBSManager.get_jobs(user=user_filter)
            
            # 1. Live Tracker for the Latest Job
            if jobs:
                # Find the latest Job ID
                latest_job = max(jobs, key=lambda x: int(x.get("id", "0").split(".")[0]) if x.get("id", "0").split(".")[0].isdigit() else 0)
                
                st.subheader(f"📡 Live Tracker: Job {latest_job.get('id')}")
                lt_col1, lt_col2, lt_col3, lt_col4 = st.columns(4)
                
                state = latest_job.get("job_state", "?")
                lt_col1.metric("Current State", state, help="R: Running, Q: Queued, H: Held")
                lt_col2.metric("Walltime Used", latest_job.get("resources_used.walltime", "00:00:00"))
                lt_col3.metric("Current RAM", latest_job.get("resources_used.mem", "0 MB"))
                lt_col4.metric("CPU Load", f"{latest_job.get('resources_used.cpupercent', '0')}%")
                
                st.divider()

            if not jobs:
                st.info(f"No active jobs found for user: {user_filter}")
            else:
                # Prepare data for display
                display_data = []
                for job in jobs:
                    display_data.append({
                        "Job ID": job.get("id"),
                        "Name": job.get("Job_Name"),
                        "State": job.get("job_state"),
                        "Queue": job.get("queue"),
                        "Walltime": job.get("resources_used.walltime", "00:00:00"),
                        "CPU %": job.get("resources_used.cpupercent", "0"),
                        "RAM": job.get("resources_used.mem", "0kb"),
                        "Node": job.get("exec_vnode", "N/A"),
                        "Comment": job.get("comment") or job.get("Comment") or job.get("depend", "")
                    })
                
                df = pd.DataFrame(display_data)
                
                # Highlight states
                def color_state(val):
                    color = 'white'
                    if val == 'R': color = '#28a745' # Green
                    elif val == 'Q': color = '#ffc107' # Yellow
                    elif val == 'H': color = '#dc3545' # Red
                    return f'background-color: {color}; color: black; font-weight: bold'

                st.subheader("Active Queue")
                st.dataframe(df.style.map(color_state, subset=['State']), width="stretch")

                # Detailed Hold analysis
                held_jobs = [j for j in jobs if j.get("job_state") == "H"]
                if held_jobs:
                    st.warning(f"Found {len(held_jobs)} jobs in HOLD state.")
                    for hj in held_jobs:
                        with st.expander(f"Hold Details: {hj.get('id')} ({hj.get('Job_Name')})"):
                            reason = hj.get('comment') or hj.get('Comment')
                            if not reason and hj.get('depend'):
                                reason = f"Dependency: {hj.get('depend')}"
                            st.error(f"Reason: {reason or 'None'}")
                            st.code(f"Error Path: {hj.get('Error_Path')}")

            # Metrics summary
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Jobs", len(jobs))
            col2.metric("Running", len([j for j in jobs if j.get("job_state") == "R"]))
            col3.metric("Held/Queued", len([j for j in jobs if j.get("job_state") in ["H", "Q"]]))
            
        render_active_queue()

    with tab2:
        st.subheader("🕰️ Persistent Job History")
        
        # 1. Bulk Scan Section
        with st.expander("🚀 Bulk Scan Job Range", expanded=False):
            col_start, col_end = st.columns(2)
            start_id = col_start.number_input("Start Job ID", value=3000, step=1)
            end_id = col_end.number_input("End Job ID", value=3050, step=1)
            
            if st.button("Start Bulk Scan", width="stretch"):
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                def update_progress(current, total):
                    progress = current / total
                    progress_bar.progress(progress)
                    status_text.text(f"Scanning Job {start_id + current - 1}... ({current}/{total})")

                PBSManager.scan_range(int(start_id), int(end_id), user_filter, progress_callback=update_progress)
                st.success(f"Scan complete! Metadata updated.")
                st.rerun()

        # 2. Trace Specific Job
        with st.expander("🔍 Quick Trace ID", expanded=False):
            trace_id = st.text_input("Enter Job ID (e.g. 3033)")
            if st.button("Trace", width="stretch"):
                if trace_id:
                    stats = PBSManager.trace_job(trace_id)
                    if stats:
                        st.success(f"Stats for {trace_id}")
                        st.json(stats)
                        # Also update metadata automatically
                        cache = PBSManager.load_metadata()
                        cache[trace_id] = stats
                        PBSManager.save_metadata(cache)
                    else:
                        st.error("Could not find job trace info.")

        # 3. History Discovery & Display
        st.divider()
        history = PBSManager.get_job_history(user=user_filter)
        
        # Also load all from metadata for this user
        all_metadata = PBSManager.load_metadata()
        user_jobs = []
        for jid, data in all_metadata.items():
            if data.get("owner") == user_filter:
                user_jobs.append({"Job ID": jid, **data})
        
        if user_jobs:
            # Sort by Job ID descending
            user_jobs.sort(key=lambda x: int(x["Job ID"]) if x["Job ID"].isdigit() else 0, reverse=True)
            df_hist = pd.DataFrame(user_jobs)
            
            # Display with data editor to allow selection/deletion
            st.write("Showing all discovered/scanned jobs for your user:")
            
            # Data editor for "Management"
            edited_df = st.data_editor(
                df_hist, 
                hide_index=True, 
                width="stretch",
                num_rows="dynamic", # Allows deleting rows
                disabled=df_hist.columns.tolist(), # Disable editing for all columns
                key="history_editor"
            )
            
            # Check for deletions
            if len(edited_df) < len(df_hist):
                # Find which ones were removed
                remaining_ids = set(edited_df["Job ID"].tolist())
                all_ids = set(df_hist["Job ID"].tolist())
                deleted_ids = all_ids - remaining_ids
                
                for did in deleted_ids:
                    PBSManager.delete_cached_job(did)
                
                st.toast(f"Deleted {len(deleted_ids)} jobs from history.")
                st.rerun()
        else:
            st.info("No historical jobs found in metadata. Try running a Bulk Scan above!")

    with tab3:
        @st.fragment(run_every=refresh_rate)
        def render_eta_estimator():
            st.subheader("⏱️ Simulation ETA Estimator")
            st.markdown("Analyze dump timestamps to estimate completion time, accounting for simulation slowdown.")
            
            base_yard = "/home/guest/palash/Granular-Chains-DEM/dumping_yard"
            
            # 1. Advanced Picker vs Quick Picker
            pick_mode = st.radio("Selection Mode", ["🔍 Auto-Detect", "📂 Manual Browser"], horizontal=True)
            
            dump_dir = ""
            
            if pick_mode == "🔍 Auto-Detect":
                # Find all directories containing chain_*.dump
                with st.spinner("Scanning dumping_yard..."):
                    detected = []
                    if os.path.exists(base_yard):
                        # We limit depth for speed
                        for root, dirs, files in os.walk(base_yard):
                            if any(f.startswith("chain_") and f.endswith(".dump") for f in files):
                                detected.append(root)
                            if len(detected) > 20: break # Safety limit
                    
                    if detected:
                        dump_dir = st.selectbox("Select an active simulation dump folder:", detected)
                    else:
                        st.warning("No active dump folders (chain_*.dump) found in dumping_yard.")
                        st.info("Try switching to 'Manual Browser' mode.")
            else:
                dump_dir = st_directory_picker("Select Dump Directory", "eta_browser_path", base_yard)

            # 2. Target Duration
            target_steps = st.number_input(
                "Target Duration (Steps to Run)", 
                value=int(st.session_state.get("eta_target_steps", 1000000)),
                step=100000,
                help="The number of steps you want this specific job to complete."
            )
            st.session_state["eta_target_steps"] = target_steps
            
            if dump_dir:
                # Save for persistence
                st.session_state["eta_dump_dir"] = dump_dir
                
                with st.spinner("Analyzing simulation progress..."):
                    eta_data = SimulationMonitor.estimate_eta(dump_dir, target_steps)
                    
                if "error" in eta_data:
                    st.error(eta_data["error"])
                else:
                    # 2. Key Metrics
                    st.divider()
                    m_col1, m_col2, m_col3 = st.columns(3)
                    
                    # Progress bar
                    progress = min(1.0, eta_data['progress_percent'] / 100.0)
                    st.progress(progress, text=f"Simulation Progress: {eta_data['progress_percent']:.1f}%")
                    
                    # Format Time Elapsed
                    e_hrs = int(eta_data['time_elapsed_hr'])
                    e_mins = int((eta_data['time_elapsed_hr'] - e_hrs) * 60)
                    m_col1.metric("Time Elapsed", f"{e_hrs}h {e_mins}m")

                    # Format ETA
                    hrs = int(eta_data['time_remaining_hr'])
                    mins = int((eta_data['time_remaining_hr'] - hrs) * 60)
                    m_col2.metric("ETA Remaining", f"{hrs}h {mins}m")
                    
                    m_col3.metric("Real-world Completion", eta_data['completion_time'].split(" ")[1], help=eta_data['completion_time'])

                    # 3. Status Card
                    st.divider()
                    s1, s2, s3 = st.columns(3)
                    
                    steps_left = eta_data['target_relative_steps'] - eta_data['current_relative_step']
                    s1.metric("Steps Remaining", f"{max(0, steps_left):,}")
                    
                    # Format Cost per Dump (using dynamic interval)
                    interval = eta_data.get('step_interval', 1000)
                    cost_m = eta_data.get('cost_per_dump', 0)
                    
                    if cost_m >= 1:
                        cost_str = f"{int(cost_m)}m"
                    else:
                        cost_str = f"{int(cost_m * 60)}s"
                    
                    s2.markdown(f"**Cost per {interval:,} steps**")
                    sc1, sc2, sc3 = s2.columns([1, 1.5, 0.5], vertical_alignment="bottom")
                    sc1.metric("Speed", cost_str, label_visibility="collapsed")
                    if eta_data.get("cost_history"):
                        sc2.line_chart(eta_data["cost_history"], height=60, width="stretch")
                    # sc3 acts as a spacer
                    
                    s3.metric("Data Points", f"{eta_data['data_points']}")
                    if eta_data.get("status") == "Target Reached":
                        st.success(f"✅ {eta_data['message']}")
                    else:
                        st.info(f"🎯 **Target Reached By:** {eta_data['completion_time']}")
                    
                    # 4. Performance Insights
                    with st.expander("📈 Performance Details"):
                        p_col1, p_col2 = st.columns(2)
                        p_col1.write(f"**Absolute Timestep:** {eta_data['current_timestep']:,}")
                        p_col1.write(f"**Steps in this Job:** {eta_data['current_relative_step']:,}")
                        p_col1.write(f"**Total Files Found:** {eta_data.get('total_files', 'N/A')}")
                        p_col2.write(f"**Fitted Points:** {eta_data['data_points']}")
                        p_col2.write(f"**Current Speed:** {eta_data['cost_per_10k_steps']:.1f} min / 10k steps")
                        p_col2.write(f"**Last File Sync:** {eta_data['last_updated']}")
                        
                        st.caption("Note: Estimation uses a quadratic fit to account for simulation slowdown as more particles enter the system.")
                        
        render_eta_estimator()

    with tab4:
        from visualiser import render_visualiser
        render_visualiser()

    with tab5:
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
            # Parse PBS Jobs for Status Overlay
            import re
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
            # Define nodes with styles
            for run_id, info in lineage.items():
                short_id = info['name'].replace('-', '_').replace('.', '_')
                node_label = f"{info['name']}<br/>(N={info['N']}, {info['steps']:,} steps)"
                seed = str(info.get('params', {}).get('seed', ''))
                
                # Style based on type
                style = ""
                if seed in running_seeds:
                    style = ":::running"
                elif info["status"] == "Archived":
                    style = ":::archived"
                elif "Flow" in info["simulation"]:
                    style = ":::flow"
                else:
                    style = ":::active"
                
                mermaid_code += f'    {short_id}["{node_label}"]{style}\n'
                mermaid_code += f'    click {short_id} call selectNode("{short_id}")\n'
            
            # Define relationships
            for run_id, info in lineage.items():
                if info["parent"] and info["parent"] in lineage:
                    parent_name = lineage[info["parent"]]['name'].replace('-', '_').replace('.', '_')
                    child_name = info['name'].replace('-', '_').replace('.', '_')
                    mermaid_code += f"    {parent_name} --> {child_name}\n"
                    
            # Inject Queued Ghost Nodes
            for ghost_id, ghost_label, parent_short_id in queued_ghosts:
                mermaid_code += f'    {ghost_id}["{ghost_label}"]:::queued\n'
                if parent_short_id:
                    mermaid_code += f"    {parent_short_id} --> {ghost_id}\n"
                mermaid_code += f'    click {ghost_id} call selectNode("{ghost_id}")\n'

            
            # Add "+" Node for new fill runs
            mermaid_code += '    NewRoot[" + Start New Fill Run "]:::new_node\n'
            mermaid_code += '    click NewRoot call selectNode("NewRoot")\n'

            
            # Define styles
            mermaid_code += "    classDef active fill:#e1f5fe,stroke:#01579b,stroke-width:2px;\n"
            mermaid_code += "    classDef archived fill:#f5f5f5,stroke:#9e9e9e,stroke-dasharray: 5 5;\n"
            mermaid_code += "    classDef flow fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;\n"
            mermaid_code += "    classDef running fill:#c8e6c9,stroke:#388e3c,stroke-width:3px;\n"
            mermaid_code += "    classDef queued fill:#ffccbc,stroke:#d32f2f,stroke-dasharray: 5 5;\n"
            mermaid_code += "    classDef new_node fill:#ffffff,stroke:#333333,stroke-width:2px,stroke-dasharray: 5 5;\n"


            # Render Mermaid using custom component
            mermaid_click_component = components.declare_component(
                "mermaid_click",
                path=os.path.join(os.path.dirname(__file__), "mermaid_component")
            )
            
            clicked_node = mermaid_click_component(mermaid_code=mermaid_code, default=None)
            
            # Handle click from the custom component natively
            if clicked_node:
                # 0. Handle New Root Click
                if clicked_node == "NewRoot":
                    if st.session_state.get("lineage_selected_parent") != "NewRoot":
                        st.session_state["lineage_selected_parent"] = "NewRoot"
                        st.rerun()

                # 1. Job Deletion Interface (for ongoing/on-hold jobs)
                if clicked_node in node_to_jobid:
                    active_jid = node_to_jobid[clicked_node]
                    st.divider()
                    st.warning(f"⚠️ **Ongoing Job Selected:** `{clicked_node}` (Job ID: `{active_jid}`)")
                    
                    # Double Confirmation Logic
                    if st.session_state.get("confirm_delete_id") == active_jid:
                        st.error("Are you absolutely sure you want to terminate this job? This will stop the simulation immediately.")
                        c1, c2 = st.columns([1, 4])
                        if c1.button("🔥 YES, TERMINATE", type="primary", use_container_width=True):
                            try:
                                PBSManager.delete_job(active_jid)
                                st.toast(f"Successfully sent qdel for {active_jid}")
                                del st.session_state["confirm_delete_id"]
                                time.sleep(1.5) # Give cluster a moment
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

                # 3. Regular Parent Selection for Launching
                if clicked_node in id_to_path:
                    new_parent = id_to_path[clicked_node]
                    # Update if changed to avoid infinite rerun loops
                    if st.session_state.get("lineage_selected_parent") != new_parent:
                        st.session_state["lineage_selected_parent"] = new_parent
                        st.rerun()

            
            # 2. Detailed Data View
            with st.expander("📄 View Detailed Parameters"):
                # Clean up for display
                display_lineage = []
                for rid, info in lineage.items():
                    display_lineage.append({
                        "Name": info["name"],
                        "Type": info["simulation"],
                        "N": info["N"],
                        "Steps": info["steps"],
                        "Status": info["status"],
                        "Path": rid
                    })
                st.dataframe(pd.DataFrame(display_lineage), width="stretch", hide_index=True)
                
            st.divider()
            st.subheader("🚀 Launch New Simulation from Lineage")
            
            # Select Parent Run
            run_options = ["NewRoot"] + list(lineage.keys())
            def format_run(rid):
                if rid == "NewRoot": return "🆕 Start New Fill Run (+)"
                return f"{lineage[rid]['name']} ({lineage[rid]['simulation']})"
            
            # Synchronize dropdown with session state (which may be set by clicking graph)
            default_index = 0
            if "lineage_selected_parent" in st.session_state:
                if st.session_state["lineage_selected_parent"] in run_options:
                    default_index = run_options.index(st.session_state["lineage_selected_parent"])
            
            selected_parent = st.selectbox(
                "Select Parent Run", 
                ["-- Select --"] + run_options, 
                index=default_index + 1 if default_index >= 0 else 0,
                format_func=lambda x: format_run(x) if x != "-- Select --" else x
            )
            
            if selected_parent != "-- Select --":
                # Ensure session state is updated if manually changed in dropdown
                st.session_state["lineage_selected_parent"] = selected_parent
                
                if selected_parent != "NewRoot":
                    parent_info = lineage[selected_parent]
                    parent_type = parent_info.get("simulation", "")
                    
                    # --- ANNOTATED LINEAGE & SNAPSHOT ---
                    note_col, snap_col = st.columns([1, 1])
                    with note_col:
                        notes_db = load_notes()
                        current_note = notes_db.get(selected_parent, "")
                        new_note = st.text_area("🗒️ Run Notes", value=current_note, height=150, help="Save observations or metadata for this run.")
                        if st.button("💾 Save Notes", key=f"save_note_{selected_parent}"):
                            save_note(selected_parent, new_note)
                            st.toast("Notes saved!")
                    
                    with snap_col:
                        with st.spinner("Generating snapshot..."):
                            try:
                                from Pulse.snapshot_helper import generate_snapshot
                                # Support manual refresh
                                force_refresh = st.button("🔄 Refresh Snapshot", key=f"refresh_snap_{selected_parent}")
                                snap_path = generate_snapshot(selected_parent, force=force_refresh)
                                
                                if snap_path:
                                    st.image(snap_path, caption=f"Last Snapshot of {parent_info['name']} (Y-Z Plane)", use_container_width=True)
                                else:
                                    st.info("No snapshot available (no dump files found).")
                            except Exception as e:
                                st.warning(f"Snapshot preview unavailable: {e}")



                    
                    st.divider()
                    
                    # Check for ongoing parent job to enforce lineage dependency
                    parent_short_id = parent_info['name'].replace('-', '_').replace('.', '_')
                    active_parent_jid = node_to_jobid.get(parent_short_id)
                    
                    if active_parent_jid:
                        st.info(f"🔗 **Lineage Dependency Detected:** This job will automatically wait for its parent (`{parent_short_id}`) to finish (PBS ID: `{active_parent_jid}`).")

                    if "Flow" in parent_type:
                        allowed_modes = ["flow_resume"]
                    else:
                        allowed_modes = ["fill_resume", "flow"]
                        
                    selected_mode = st.radio("Select Simulation Mode", allowed_modes, horizontal=True)
                else:
                    parent_info = {"name": "S000000", "simulation": "None", "N": 4, "params": {"seed": "000000"}}
                    selected_mode = "fill"
                    st.info("🆕 **Creating New Root Fill Run** (No Parent)")
                    active_parent_jid = None
                
                with st.form("launch_sim_form"):
                    st.markdown("### 🌍 Global Parameters")
                    gc1, gc2, gc3, gc4 = st.columns(4)
                    walltime = gc1.text_input("Walltime", value="24:00:00")
                    ppn = gc2.number_input("PPN", value=16, step=1)
                    mem = gc3.text_input("Memory", value="16gb")
                    max_concurrent = gc4.number_input("Max Concurrent Jobs", value=4, min_value=1, step=1, help="If you submit more jobs than this limit, they will be queued with depend=afterany.")
                    
                    gc5, gc6, gc7, gc8 = st.columns(4)
                    num_procs = gc5.number_input("Num Procs", value=8, step=1)
                    num_threads = gc6.number_input("Num Threads", value=1, step=1)
                    dt = gc7.number_input("dt", value=1e-06, format="%e")
                    viscosity = gc8.number_input("Viscosity", value=0.001, format="%f")
                    
                    seed = st.number_input("Seed", value=random.randint(100000, 999999), step=1)
                    
                    st.markdown(f"### ⚙️ Mode-Specific Parameters ({selected_mode})")
                    
                    # Store mode specific inputs
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

                        fc7, fc8 = st.columns(2)
                        mode_params["mode"] = fc7.selectbox("Pouring Mode", ["2D_stacked", "2D_worst_case"], index=0)
                        mode_params["dump_file"] = fc8.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                        
                        mode_params["hopper_template_data"] = st.text_input("Hopper Template Data", value="simulation_geometries/2D_hopper_with_orifice_cover.inc")
                        
                        mode_params["simulation"] = "Hopper_Fill"
                        mode_params["no-vtk"] = True

                    elif selected_mode == "fill_resume":
                        rc1, rc2 = st.columns(2)
                        mode_params["relax_steps"] = rc1.number_input("Relax Steps", value=1000000, step=100000)
                        mode_params["dump_file"] = rc2.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                        mode_params["simulation"] = "Hopper_Fill_Resume"
                        mode_params["restart_path"] = selected_parent
                        
                    elif selected_mode == "flow":
                        fc1, fc2, fc3, fc4 = st.columns(4)
                        # Parametric Sweep Support
                        freq_input = fc1.text_input("Frequency(s)", value="5.0", help="Single value (5.0) or comma separated (5.0, 10.0, 15.0)")
                        amp_input = fc2.text_input("Amplitude(s)", value="0.01", help="Single value (0.01) or comma separated (0.01, 0.02)")
                        
                        mode_params["run_steps"] = fc3.number_input("Run Steps", value=2000000, step=100000)
                        mode_params["osc_dir"] = fc4.text_input("Oscillation Dir", value="z")
                        mode_params["source_dir"] = selected_parent
                        
                        # Parse inputs for display preview
                        freqs = [f.strip() for f in freq_input.split(",") if f.strip()]
                        amps = [a.strip() for a in amp_input.split(",") if a.strip()]
                        num_sweep = len(freqs) * len(amps)
                        
                        if num_sweep > 1:
                            st.warning(f"🎰 **Parametric Sweep:** This will submit **{num_sweep}** separate jobs.")
                        
                        mode_params["_freq_list"] = freqs
                        mode_params["_amp_list"] = amps
                        
                    elif selected_mode == "flow_resume":
                        frc1, frc2 = st.columns(2)
                        mode_params["run_steps"] = frc1.number_input("Run Steps", value=1000000, step=100000)
                        mode_params["restart_path"] = selected_parent
                        
                    submit_btn = st.form_submit_button("🚀 Submit to PBS")
                    
                    if submit_btn:
                        parent_seed = str(parent_info.get("params", {}).get("seed", "000000"))[-6:]
                        
                        # Handle Parametric Sweep for Flow
                        jobs_to_submit = []
                        if selected_mode == "flow" and len(mode_params.get("_freq_list", [])) * len(mode_params.get("_amp_list", [])) > 1:
                            for f in mode_params["_freq_list"]:
                                for a in mode_params["_amp_list"]:
                                    child_seed = random.randint(100000, 999999)
                                    job_name = f"F{parent_seed}_{str(child_seed)[-6:]}"
                                    
                                    # Create deep copy/clone of params for this job
                                    job_params = {
                                        "num_procs": int(num_procs),
                                        "num_threads": int(num_threads),
                                        "dt": dt,
                                        "viscosity": viscosity,
                                        "seed": child_seed,
                                        **mode_params
                                    }

                                    # Override specific sweep values
                                    job_params["freq"] = float(f)
                                    job_params["amp"] = float(a)
                                    job_params["simulation"] = f"Flow_Study_N{parent_info.get('N', 0)}_F{f}_A{a}"
                                    # Cleanup internal keys
                                    job_params.pop("_freq_list", None)
                                    job_params.pop("_amp_list", None)
                                    
                                    jobs_to_submit.append({
                                        "name": job_name,
                                        "type": selected_mode,
                                        "walltime": walltime,
                                        "ppn": int(ppn),
                                        "mem": mem,
                                        "dependency": active_parent_jid,
                                        "params": job_params
                                    })
                        else:
                            # Single job submission
                            new_seed = str(seed)[-6:]
                            if selected_mode == "fill_resume": prefix = "R"
                            elif selected_mode == "flow": prefix = "F"
                            elif selected_mode == "flow_resume": prefix = "FR"
                            else: prefix = "S"
                            
                            job_name = f"{prefix}{parent_seed}_{new_seed}"
                            
                            # Finalize params
                            final_params = {
                                "num_procs": int(num_procs),
                                "num_threads": int(num_threads),
                                "dt": dt,
                                "viscosity": viscosity,
                                "seed": int(seed),
                                **mode_params
                            }
                            if selected_mode == "flow":
                                final_params["freq"] = float(mode_params.get("_freq_list", [0])[0])
                                final_params["amp"] = float(mode_params.get("_amp_list", [0])[0])
                                final_params["simulation"] = f"Flow_Study_N{parent_info.get('N', 0)}_F{final_params['freq']}_A{final_params['amp']}"
                            
                            final_params.pop("_freq_list", None)
                            final_params.pop("_amp_list", None)

                            jobs_to_submit.append({
                                "name": job_name,
                                "type": selected_mode,
                                "walltime": walltime,
                                "ppn": int(ppn),
                                "mem": mem,
                                "dependency": active_parent_jid,
                                "params": final_params
                            })
                        
                        try:
                            import importlib
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
                                errf.write(f"\\n--- Error at {time.strftime('%Y-%m-%d %H:%M:%S')} ---\\n")
                                errf.write(traceback.format_exc())
                            st.error(f"Error submitting job: {e}. Check submit_error.log for details.")

    # 5. System Status (Master Node only)
    if "master" in socket.gethostname():

        @st.fragment(run_every=refresh_rate)
        def render_system_status():
            st.divider()
            st.subheader("🖥️ Master Node System Status")
            
            col_cpu, col_mem, col_gpu1, col_gpu2 = st.columns(4)
            
            # CPU Usage
            cpu_usage = psutil.cpu_percent()
            col_cpu.metric("Total CPU Load", f"{cpu_usage}%")
            
            # Memory Usage
            mem = psutil.virtual_memory()
            col_mem.metric("Total RAM Usage", f"{mem.percent}%", f"{mem.used // (1024**3)} GB used")
            
            # GPU Usage
            try:
                gpu_info = subprocess.getoutput("nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits").split(",")
                if len(gpu_info) >= 2:
                    col_gpu1.metric("GPU Utilization", f"{gpu_info[0].strip()}%")
                    col_gpu2.metric("GPU Memory", f"{gpu_info[1].strip()} MiB")
            except:
                col_gpu1.write("GPU monitoring unavailable")

            # 6. Thermal Status
            st.markdown("### 🌡️ Thermal Status")
            col_t_cpu, col_t_gpu, col_t_empty1, col_t_empty2 = st.columns(4)
            
            temps = PBSManager.get_node_temperatures()
            if "cpu" in temps:
                col_t_cpu.metric("CPU Temp", f"{temps['cpu']:.1f} °C", delta=f"{temps['cpu']-60:.1f} °C" if temps['cpu'] > 60 else None, delta_color="inverse")
            if "gpu" in temps:
                col_t_gpu.metric("GPU Temp", f"{temps['gpu']:.1f} °C", delta=f"{temps['gpu']-70:.1f} °C" if temps['gpu'] > 70 else None, delta_color="inverse")

        render_system_status()
