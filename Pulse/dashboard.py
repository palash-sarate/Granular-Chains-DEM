import streamlit as st
import pandas as pd
import sys
import os
import json
import glob
from pathlib import Path

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
import datetime

# Add project root to sys.path to allow importing from analysis module
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pulse_core import PBSManager, SimulationMonitor, SyncManager, AutoPilotManager
import streamlit.components.v1 as components
from analysis.controllers.sim_data import SimDataController
from analysis.controllers.renderer import SimulationRenderer
from analysis.controllers.highlighter import HighlightController
from stpyvista import stpyvista
import socket
import random

scratch_dir = os.path.join(ROOT_DIR, "scratch")
if scratch_dir not in sys.path:
    sys.path.append(scratch_dir)
import submit_flexible
from lineage_tracker import scan_dumping_yard
import importlib
import analysis.controllers.unified_renderer
importlib.reload(analysis.controllers.unified_renderer)
from analysis.controllers.unified_renderer import UnifiedRenderer
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

VIZ_SETTINGS_FILE = os.path.join(ROOT_DIR, "Pulse/viz_settings.json")
def load_viz_settings():
    if os.path.exists(VIZ_SETTINGS_FILE):
        try:
            with open(VIZ_SETTINGS_FILE, "r") as f: return json.load(f)
        except: pass
    return {"offset": [0.0, 0.0, 0.0], "zoom": 1.0, "camera": None, "width": 1280, "height": 720}

def save_viz_settings(offset, zoom, camera=None, width=1280, height=720):
    try:
        with open(VIZ_SETTINGS_FILE, "w") as f:
            json.dump({"offset": offset, "zoom": zoom, "camera": camera, "width": width, "height": height}, f)
    except: pass

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

    # --- Navigation ---
    # We use a horizontal radio button to act as 'True Navigation'
    # Standard st.tabs execute all tabs' code on every rerun, which causes the 1-minute lag for the 3D visualiser.
    # This radio selector ensures ONLY the active page's code is executed.
    st.write('<style>div.row-widget.stRadio > div{flex-direction:row; justify-content: center; gap: 20px;} div.row-widget.stRadio label{background: #f0f2f6; padding: 10px 20px; border-radius: 5px; cursor: pointer;} div.row-widget.stRadio div[role="radiogroup"] > label[data-baseweb="radio"]{background: #f0f2f6; border: 1px solid #ddd;}</style>', unsafe_allow_html=True)
    
    # --- Navigation Persistence ---
    pages = ["📊 Active Queue", "🧬 Lineage", "🤖 Auto-Pilot", "🎬 Visualization", "🎥 Visualizer", "⏱️ ETA", "🕰️ History", "🔄 Sync"]
    
    # Initialize from URL or default
    query_nav = st.query_params.get("tab", pages[0])
    if query_nav not in pages:
        query_nav = pages[0]
    
    current_idx = pages.index(query_nav)
    
    nav = st.radio(
        "Select Page",
        pages,
        index=current_idx,
        horizontal=True,
        label_visibility="collapsed",
        key="nav_selector"
    )

    # Sync URL with selection
    if nav != query_nav:
        st.query_params["tab"] = nav

    st.divider()

    # --- AUTO-PILOT SETTINGS LOADER ---
    AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
    def load_auto_pilot():
        if os.path.exists(AUTO_PILOT_FILE):
            try:
                with open(AUTO_PILOT_FILE, "r") as f: return json.load(f)
            except: pass
        return {"settings": {"enabled": False}, "goals": {}}

    def save_auto_pilot(data):
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)

    if nav == "📊 Active Queue":
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
                        "Polite": "😇" if int(job.get("Priority", 0)) < 0 else "⚡",
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

    elif nav == "🕰️ History":
        st.subheader("🕰️ Comprehensive Job History")
        
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
                # ... deletion logic ...
                remaining_ids = set(edited_df["Job ID"].tolist())
                all_ids = set(df_hist["Job ID"].tolist())
                deleted_ids = all_ids - remaining_ids
                for did in deleted_ids:
                    PBSManager.delete_cached_job(did)
                st.toast(f"Deleted {len(deleted_ids)} jobs from history.")
                st.rerun()
        else:
            st.info("No historical jobs found in metadata. Try running a Bulk Scan above!")

    elif nav == "⏱️ ETA":
        def render_eta_section():
            st.subheader("⏱️ Simulation ETA Estimator")
            st.markdown("Analyze dump timestamps to estimate completion time, accounting for simulation slowdown.")
            base_yard = "/home/guest/palash/Granular-Chains-DEM/dumping_yard"
            pick_mode = st.radio("Selection Mode", ["🔍 Auto-Detect", "📂 Manual Browser"], horizontal=True, key="eta_pick_mode")
            dump_dir = ""
            if pick_mode == "🔍 Auto-Detect":
                with st.spinner("Scanning dumping_yard..."):
                    detected = []
                    if os.path.exists(base_yard):
                        for root, dirs, files in os.walk(base_yard):
                            if any(f.startswith("chain_") and f.endswith(".dump") for f in files):
                                detected.append(root)
                            if len(detected) > 20: break
                    if detected:
                        dump_dir = st.selectbox("Select an active simulation dump folder:", detected)
            else:
                dump_dir = st_directory_picker("Select Dump Directory", "eta_browser_path", base_yard)

            target_steps = st.number_input("Target Duration (Steps to Run)", value=1000000, step=100000)
            if dump_dir and st.button("Calculate ETA", use_container_width=True):
                eta_data = SimulationMonitor.estimate_eta(dump_dir, target_steps)
                if "error" in eta_data:
                    st.error(eta_data["error"])
                else:
                    # 2. Key Metrics
                    st.divider()
                    m_col1, m_col2, m_col3 = st.columns(3)

                    # Format Time Elapsed
                    e_hrs = int(eta_data['time_elapsed_hr'])
                    e_mins = int((eta_data['time_elapsed_hr'] - e_hrs) * 60)
                    m_col1.metric("Time Elapsed", f"{e_hrs}h {e_mins}m")

                    # Format ETA
                    hrs = int(eta_data['time_remaining_hr'])
                    mins = int((eta_data['time_remaining_hr'] - hrs) * 60)
                    m_col2.metric("ETA Remaining", f"{hrs}h {mins}m")
                    
                    m_col3.metric("Completion Time", eta_data['completion_time'].split(" ")[1], help=eta_data['completion_time'])

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
                        sc2.line_chart(eta_data["cost_history"], height=60, use_container_width=True)
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
                        
        render_eta_section()

    elif nav == "🎥 Visualizer":
        from Pulse.visualiser import render_visualiser
        render_visualiser()

    elif nav == "🎬 Visualization":
        st.subheader("🎬 Lineage Visualization & Movie Generation")
        st.markdown("Create high-quality movies across multiple simulation stages with automatic cloud-restoration.")
        
        from Pulse.viz_manager import VizManager
        lineage = PBSManager.load_lineage()
        
        if not lineage:
            st.info("No lineage data found. Please scan your simulations in the 'Lineage' tab first.")
        else:
            # 1. Chain Selection
            st.write("### 1. Select Lineage Chain")
            all_paths = sorted(list(lineage.keys()))
            names = {p: lineage[p]['name'] for p in all_paths}
            
            c1, c2 = st.columns(2)
            start_node = c1.selectbox("Start Node", all_paths, format_func=lambda x: names[x], index=0)
            
            # Filter end_node options based on descendants of start_node (including start_node itself)
            possible_ends = [start_node] + VizManager.get_descendants(start_node, lineage)
            end_node = c2.selectbox("End Node", possible_ends, format_func=lambda x: names[x], index=len(possible_ends)-1)

            if end_node:
                chain = VizManager.get_lineage_chain(start_node, end_node)
            else:
                chain = []
            
            if not chain:
                st.error("No valid direct lineage chain found. Select a different Start/End node.")
            else:
                st.success(f"Resolved Chain: {' → '.join([names[p] for p in chain])}")
                
                # Check for missing synced runs in the chain
                missing_synced = []
                for p in chain:
                    if not os.path.exists(p) and lineage.get(p, {}).get("sync_status") == "Synced":
                        missing_synced.append(p)
                
                if missing_synced:
                    st.warning(f"⚠️ **{len(missing_synced)} runs** in this chain are synced but not local. They must be restored for scrubbing and preview.")
                    if st.button("📥 Restore Missing Data in Chain", use_container_width=True, type="primary"):
                        progress_bar = st.progress(0)
                        status_text = st.empty()
                        for i, p in enumerate(missing_synced):
                            status_text.text(f"Restoring {os.path.basename(p)}... ({i+1}/{len(missing_synced)})")
                            SyncManager.restore_run(p)
                            progress_bar.progress((i + 1) / len(missing_synced))
                        st.success("Restoration complete!")
                        st.session_state.pop("viz_chain_frames", None)
                        st.session_state.pop("last_viz_chain", None) # Force re-scan
                        st.rerun()
                
                # Frame Scanning & Scrubbing
                if "last_viz_chain" not in st.session_state or st.session_state["last_viz_chain"] != chain:
                    st.session_state["last_viz_chain"] = chain
                    with st.spinner("Scanning chain frames..."):
                        st.session_state["viz_chain_frames"] = VizManager.get_chain_frames(chain)
                
                all_frames = st.session_state.get("viz_chain_frames", [])
                
                if not all_frames:
                    st.warning("No dump frames found in this chain. Ensure simulations have data files.")
                else:
                    st.divider()
                    st.write("### ⏱️ Timeline Scrubber")
                    scrub_idx = st.slider("Scrub through all frames in chain", 0, len(all_frames)-1, 0, 
                                         format="Frame %d", help="Drag to select a specific moment from any run in the chain.")
                    target_frame = all_frames[scrub_idx]
                    viz_dt = st.session_state.get("viz_dt", 1e-6)
                    st.caption(f"📍 **Selected**: `{os.path.basename(target_frame['path'])}` | Timestep: `{target_frame['ts']}` | `{target_frame['ts']*viz_dt:.3f}s`")

                # 2. Configuration
                st.divider()
                st.write("### 2. Rendering Configuration")
                
                col1, col2 = st.columns(2)
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write("**Canvas & Resolution**")
                    c_w, c_h = st.columns(2)
                    
                    saved_settings = load_viz_settings()
                    canvas_width = c_w.number_input("Width (px)", 320, 3840, saved_settings.get("width", 1280))
                    canvas_height = c_h.number_input("Height (px)", 240, 2160, saved_settings.get("height", 720))
                    fps = st.number_input("FPS", 1, 60, 30)
                    dt = st.number_input("Time Step (dt)", 1e-8, 1e-3, 1e-6, format="%.2e", key="viz_dt")

                    st.write("**Scene Translation (Offset)**")
                    off_x, off_y, off_z = st.columns(3)
                    
                    # Ensure persistence using file-based settings and session state
                    if "viz_scene_offset" not in st.session_state:
                        saved_settings = load_viz_settings()
                        st.session_state["viz_scene_offset"] = saved_settings.get("offset", [0.0, 0.0, 0.0])
                        st.session_state["viz_zoom"] = saved_settings.get("zoom", 1.0)
                        if saved_settings.get("camera"):
                            st.session_state["viz_last_cam_pos"] = saved_settings["camera"]

                    trans_x = off_x.number_input("Move X", value=st.session_state["viz_scene_offset"][0], format="%.3f")
                    trans_y = off_y.number_input("Move Y", value=st.session_state["viz_scene_offset"][1], format="%.3f")
                    trans_z = off_z.number_input("Move Z", value=st.session_state["viz_scene_offset"][2], format="%.3f")
                    
                    # Sync to session state
                    new_offset = [trans_x, trans_y, trans_z]
                    
                    # Store these in session state for UI responsiveness
                    st.session_state["viz_canvas_res"] = [canvas_width, canvas_height]
                    st.session_state["viz_scene_offset"] = new_offset

                    zoom = st.number_input("Base Zoom", min_value=0.01, max_value=100.0, value=st.session_state["viz_zoom"], format="%.2f", help="Set the camera magnification level.")
                    
                    # If anything changed, save to disk (excluding camera, which is saved via button)
                    if zoom != st.session_state["viz_zoom"] or \
                       new_offset != saved_settings["offset"] or \
                       canvas_width != saved_settings.get("width") or \
                       canvas_height != saved_settings.get("height"):
                        
                        st.session_state["viz_zoom"] = zoom
                        save_viz_settings(new_offset, zoom, camera=st.session_state.get("viz_last_cam_pos"), width=canvas_width, height=canvas_height)
                
                with col2:
                    st.write("**Visual Elements**")
                    
                    # Discover all available VTKs in the chain
                    all_chain_vtks = []
                    for p in chain:
                        geo_dir = os.path.join(p, "Geometry_vtk")
                        if os.path.exists(geo_dir):
                            all_chain_vtks.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))
                    
                    vtk_basenames = sorted(list(set([os.path.basename(v) for v in all_chain_vtks])))
                    selected_vtk_names = st.multiselect("Visible Geometry Layers", vtk_basenames, default=vtk_basenames)
                    
                    show_geo = st.checkbox("Show Geometry (Global)", value=True)
                    show_interactive = st.checkbox("Enable Interactive 3D View", value=False)
                    show_axes = st.checkbox("Show Corner Axes", value=True)
                    show_grid = st.checkbox("Show 3D Grid", value=False)

                # Interactive 3D Configuration Preview
                st.divider()
                st.write("#### 🕹️ Interactive Configuration Preview")
                st.caption("Rotate, zoom, and pan to set the perfect camera view for your movie.")
                
                # Use the target_frame from the scrubber above
                target_dump = None
                if target_frame:
                    f_path, f_ts = target_frame['path'], target_frame['ts']
                    s_dirs = [os.path.join(f_path, "chain"), f_path]
                    for sd in s_dirs:
                        if os.path.isdir(sd):
                            matches = glob.glob(os.path.join(sd, f"*{f_ts}.dump"))
                            if matches:
                                target_dump = matches[0]
                                break
                
                # Collect geometry files
                vtk_files = []
                for p in chain:
                    geo_dir = os.path.join(p, "Geometry_vtk")
                    if os.path.exists(geo_dir):
                        vtk_files.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))

                # Check for missing data in the interactive preview
                geo_vtk_dir = os.path.join(target_frame['path'], "Geometry_vtk") if target_frame else None
                missing_vtks = geo_vtk_dir and not os.path.exists(geo_vtk_dir)

                if target_frame and not os.path.exists(target_frame['path']):
                    st.error(f"Data not local for {os.path.basename(target_frame['path'])}. Restore it first for a preview.")
                    if st.button(f"📥 Restore {os.path.basename(target_frame['path'])} Now", use_container_width=True, type="primary"):
                        with st.spinner("Restoring..."):
                            success, msg = SyncManager.restore_run(target_frame['path'])
                            if success:
                                st.success("Restore complete!")
                                st.rerun()
                            else:
                                st.error(msg)
                elif target_frame:
                    # --- Advanced Geometry Tools Expander ---
                    with st.expander("🛠️ Advanced Geometry Tools", expanded=missing_vtks):
                        st.write("#### 🏗️ VTK Mesh (Re)generation")
                        st.caption("Reconstruct the hopper geometry from simulation metadata. Useful if files are missing or you want higher resolution.")
                        
                        v_col1, v_col2 = st.columns(2)
                        lattice_spacing = v_col1.number_input("Lattice Spacing (m)", 0.0001, 0.02, 0.002, format="%.4f", help="Resolution of the sampling grid. Lower = Sharper but Slower.")
                        recon_radius = v_col2.number_input("Recon. Radius (m)", 0.0, 0.1, 0.0, format="%.4f", help="Radius for surface reconstruction. Leave 0.0 for auto-calculation (1.2 * spacing).")
                        
                        if st.button("🔄 (Re)generate Geometry Meshes", use_container_width=True, type="primary" if missing_vtks else "secondary"):
                            with st.spinner("Regenerating geometry meshes..."):
                                try:
                                    import warnings
                                    warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
                                    
                                    from simulation.grid_hopper_manager import GridHopperManager
                                    from simulation.runner import SimulationRunner
                                    from analysis.geometry_extractor import GeometryExtractor
                                    
                                    # 1. Robust Metadata Search
                                    meta_path = None
                                    curr_p = Path(target_frame['path'])
                                    for _ in range(4):
                                        for m_name in ["grid_metadata.json", "metadata.json"]:
                                            test_p = curr_p / m_name
                                            if test_p.exists():
                                                meta_path = test_p
                                                break
                                        if meta_path: break
                                        curr_p = curr_p.parent
                                        if curr_p.name == "dumping_yard": break

                                    if not meta_path:
                                        st.error("Metadata not found. Cannot regenerate geometry.")
                                    else:
                                        with open(meta_path, 'r') as f:
                                            meta = json.load(f)
                                        
                                        runner = SimulationRunner(lammps_executable="lmp") 
                                        mgr = GridHopperManager(runner)
                                        
                                        inc_name = meta.get("geometry_inc", "replicated_geometry.inc")
                                        inc_path = os.path.join(target_frame['path'], os.path.basename(inc_name))
                                        
                                        if not os.path.exists(inc_path):
                                            mgr._generate_replicated_geometry(
                                                setup_path=Path(meta["hopper_template_data"]),
                                                n_hoppers=len(meta["metadata"]),
                                                spacing=meta["spacing"],
                                                job_dir=Path(target_frame['path']),
                                                normalized_geo_vars=meta["geometry_vars"],
                                                inc_name=os.path.basename(inc_name)
                                            )
                                        
                                        extractor = GeometryExtractor(lammps_cmd="lmp")
                                        env = meta.get("envelope", {})
                                        bounds = [env['total_bounds']['x'][0], env['total_bounds']['x'][1],
                                                  env['total_bounds']['y'][0], env['total_bounds']['y'][1],
                                                  env['total_bounds']['z'][0], env['total_bounds']['z'][1]] if 'total_bounds' in env else None
                                        
                                        extractor.extract(
                                            inc_file=Path(inc_path),
                                            outdir=Path(geo_vtk_dir),
                                            auto_vis=True,
                                            combined=False,
                                            bounds=bounds,
                                            spacing=lattice_spacing,
                                            radius=recon_radius if recon_radius > 0 else None
                                        )
                                        st.success("Geometry VTKs generated successfully!")
                                        st.rerun()
                                except Exception as e:
                                    st.error(f"Failed to generate VTKs: {e}")
                                    import traceback
                                    st.code(traceback.format_exc())

                    if missing_vtks:
                        st.warning("⚠️ **Geometry VTKs Missing**: Hopper boundaries will not be visible. Use the tools above to generate them.")

                    # Filter VTKs based on selection
                    active_vtk_files = []
                    for v in vtk_files:
                        if os.path.basename(v) in selected_vtk_names:
                            active_vtk_files.append(v)

                    # Setup Plotter with custom resolution
                    plotter = UnifiedRenderer.setup_plotter(window_size=st.session_state["viz_canvas_res"])

                    UnifiedRenderer.apply_scene(
                        plotter,
                        vtk_files=active_vtk_files,
                        dump_path=target_dump,
                        show_geometry=show_geo,
                        show_particles=True,
                        show_axes=show_axes,
                        show_grid=show_grid,
                        camera_state=st.session_state.get("viz_last_cam_pos"),
                        offset=st.session_state["viz_scene_offset"],
                        zoom=zoom
                    )
                    
                    # Viewfinder Border CSS
                    vw, vh = st.session_state["viz_canvas_res"]
                    st.markdown(f"""
                        <style>
                        .viz-viewfinder {{
                            border: 2px solid #333;
                            border-radius: 5px;
                            padding: 10px;
                            background: #000;
                            width: 100%;
                            display: flex;
                            justify-content: center;
                            align-items: center;
                            margin-bottom: 10px;
                        }}
                        .viz-info-overlay {{
                            color: #666;
                            font-family: monospace;
                            font-size: 0.8rem;
                            text-align: center;
                            margin-top: 5px;
                        }}
                        </style>
                        <div class="viz-info-overlay">Recording Frame: {vw}x{vh}px (Aspect Ratio: {round(vw/vh, 2)})</div>
                    """, unsafe_allow_html=True)

                    # Display interactive plotter in a bordered 'Viewfinder' container (if enabled)
                    if show_interactive:
                        with st.container(border=True):
                            st_state = stpyvista(plotter, key="viz_preview_plot")
                        
                        # Persist camera state on interaction
                        if st_state is not None:
                            # Extract camera position safely
                            new_pos = st_state.get("camera_position")
                            if new_pos:
                                # Subtract current offset before saving to maintain 'subject-relative' view
                                import numpy as np
                                off = np.array(st.session_state["viz_scene_offset"])
                                clean_pos = list(np.array(new_pos[0]) - off)
                                clean_fp = list(np.array(new_pos[1]) - off)
                                st.session_state["viz_last_cam_pos"] = [clean_pos, clean_fp, new_pos[2]]
                    else:
                        st.info("💡 **Interactive Preview Disabled**: Use 'High-Res Snapshot' below to see your current scene setup.")
                    
                    c_cam1, c_cam2, c_cam3 = st.columns(3)
                    if c_cam1.button("💾 Save View", use_container_width=True, help="Save the current camera angle as the default for this session and future restarts."):
                        save_viz_settings(
                            st.session_state["viz_scene_offset"], 
                            st.session_state["viz_zoom"], 
                            camera=st.session_state.get("viz_last_cam_pos"),
                            width=canvas_width,
                            height=canvas_height
                        )
                        st.success("Camera orientation and resolution saved to disk!")
                    
                    if c_cam2.button("🔄 Reset View", use_container_width=True):
                        st.session_state.pop("viz_last_cam_pos", None)
                        st.rerun()
                    
                    if c_cam3.button("📸 High-Res Snapshot", use_container_width=True):
                        with st.spinner("Generating High-Res PNG..."):
                            params = {
                                "chain_paths": chain,
                                "fps": fps,
                                "dt": dt,
                                "show_geometry": show_geo,
                                "selected_vtks": selected_vtk_names,
                                "show_axes": show_axes,
                                "show_grid": show_grid,
                                "camera_position": st.session_state.get("viz_last_cam_pos"),
                                "offset": st.session_state["viz_scene_offset"],
                                "resolution": st.session_state["viz_canvas_res"],
                                "zoom": zoom
                            }
                            shot_path = VizManager.generate_high_res_snapshot(params, target_frame)
                            st.session_state["last_viz_snapshot"] = shot_path

                    if "last_viz_snapshot" in st.session_state:
                        st.write("---")
                        st.write("🖼️ **Latest High-Res Snapshot**")
                        st.image(st.session_state["last_viz_snapshot"], use_container_width=True)
                        if st.button("🗑️ Clear Snapshot"):
                            st.session_state.pop("last_viz_snapshot")
                            st.rerun()

                    st.info("📌 **Capture Active**: Your current orientation, zoom, resolution, and offset will be used for the movie.")
                else:
                    st.info("No frames found in the selected lineage. Please check if simulation data exists.")

                # 3. Job Submission
                st.divider()
                st.write("### 3. Execution")
                
                # Check for existing job
                lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
                is_running = False
                job_id = ""
                if os.path.exists(lock_path):
                    try:
                        with open(lock_path, "r") as f:
                            content = f.read().strip()
                            if content.startswith("PBS:"):
                                job_id = content.replace("PBS:", "")
                                is_running = True
                    except: pass

                if is_running:
                    st.warning(f"⚠️ **Visualization job is currently running** (PBS ID: `{job_id}`)")
                    if st.button("🛑 Cancel Viz Job", type="secondary"):
                        subprocess.run(["qdel", job_id])
                        if os.path.exists(lock_path): os.remove(lock_path)
                        st.rerun()
                else:
                    output_name = st.text_input("Movie Name", value=f"{names[start_node]}_to_{names[end_node]}.mp4")
                    
                    if st.button("🎬 Submit Visualization Job", type="primary", use_container_width=True):
                        params = {
                            "chain_paths": chain,
                            "output_name": output_name,
                            "fps": fps,
                            "dt": dt,
                            "show_geometry": show_geo,
                            "selected_vtks": selected_vtk_names,
                            "show_axes": show_axes,
                            "show_grid": show_grid,
                            "camera_position": st.session_state.get("viz_last_cam_pos"),
                            "offset": st.session_state["viz_scene_offset"],
                            "resolution": st.session_state["viz_canvas_res"],
                            "zoom": zoom
                        }
                        
                        success, msg = VizManager.submit_viz_job(params)
                        if success:
                            st.success(msg)
                            st.rerun()
                        else:
                            st.error(msg)
                
                # Show Restore Button if needed
                if st.session_state.get("viz_restore_error"):
                    st.error(st.session_state["viz_restore_error"])
                    t_path = st.session_state.get("viz_restore_path")
                    if t_path and st.button(f"📥 Restore {os.path.basename(t_path)} Now", use_container_width=True, type="primary"):
                        with st.spinner("Restoring from Google Drive..."):
                            success, msg = SyncManager.restore_run(t_path)
                            if success:
                                st.success("Restore complete! You can now preview the snapshot.")
                                st.session_state.pop("viz_restore_error", None)
                                st.rerun()
                            else:
                                st.error(msg)

                st.info("💡 Note: Missing simulation data will be automatically restored from Google Drive and cleaned up after the job finishes.")

    elif nav == "🧬 Lineage":
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
                seed = str(info.get('params', {}).get('seed', ''))
                
                # Base Style
                on_disk = os.path.exists(run_id)
                is_auto = run_id in auto_goals
                
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
            mermaid_code += "    classDef synced fill:#80cbc4,stroke:#00695c,stroke-width:2px;\n" # Deep Teal
            mermaid_code += "    classDef cleared fill:#b39ddb,stroke:#512da8,stroke-width:2px;\n" # Deep Purple
            mermaid_code += "    classDef flow fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;\n"
            mermaid_code += "    classDef running fill:#a5d6a7,stroke:#2e7d32,stroke-width:3px;\n"
            mermaid_code += "    classDef queued fill:#ffccbc,stroke:#d32f2f,stroke-dasharray: 5 5;\n"
            mermaid_code += "    classDef new_node fill:#ffffff,stroke:#333333,stroke-width:2px,stroke-dasharray: 5 5;\n"
            mermaid_code += "    classDef auto_pilot fill:#b2ebf2,stroke:#00acc1,stroke-width:3px;\n"
            mermaid_code += "    classDef selected stroke:#ff9800,stroke-width:4px;\n"


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

                if clicked_node in id_to_path:
                    rid = id_to_path[clicked_node]
                    info = lineage[rid]
                    
                    # 3.1 Display Sync & Restore Options
                    sync_status = info.get("sync_status", "Local")
                    is_local = os.path.exists(rid)
                    
                    st.divider()
                    st.markdown(f"### 📦 Storage Status: **{sync_status}**")
                    
                    # Display Run Parameters for the selected node
                    p_info = info.get('params', {})
                    if p_info:
                        st.markdown("#### 📋 Run Parameters")
                        # Display main params in columns
                        main_params = {k: v for k, v in p_info.items() if k != "geometry_vars" and v is not None}
                        if main_params:
                            p_cols = st.columns(min(len(main_params), 4))
                            for i, (k, v) in enumerate(main_params.items()):
                                p_cols[i % 4].metric(k.capitalize(), str(v))
                        
                        # Dedicated section for geometry variables
                        geo_vars = p_info.get("geometry_vars")
                        if geo_vars and geo_vars != {}:
                            with st.expander("🛠️ Geometry Overrides", expanded=True):
                                st.json(geo_vars)
                    
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

                    # 4. Regular Parent Selection for Launching
                    new_parent = id_to_path[clicked_node]
                    
                    # Initialize multiselect state if not present
                    if "lineage_selected_parents" not in st.session_state:
                        st.session_state["lineage_selected_parents"] = []
                    
                    # Toggle selection
                    if new_parent in st.session_state["lineage_selected_parents"]:
                        st.session_state["lineage_selected_parents"].remove(new_parent)
                    else:
                        st.session_state["lineage_selected_parents"].append(new_parent)
                    st.rerun()

            
            # 2. Detailed Data View
            with st.expander("📄 View Detailed Parameters"):
                # Clean up for display
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
            
            # Initialize session state for multiselect
            if "lineage_selected_parents" not in st.session_state:
                st.session_state["lineage_selected_parents"] = []
            
            selected_parents = st.multiselect(
                "Select Parent Run(s)", 
                run_options, 
                default=st.session_state["lineage_selected_parents"],
                format_func=lambda x: format_run(x)
            )
            
            # Update session state source of truth
            st.session_state["lineage_selected_parents"] = selected_parents
            
            if selected_parents:
                # Validation and Categorization
                if "NewRoot" in selected_parents and len(selected_parents) > 1:
                    st.error("⚠️ **Conflict:** 'New Root' cannot be combined with existing parent runs. Please select one or the other.")
                    st.stop()

                if selected_parents == ["NewRoot"]:
                    selected_mode = "fill"
                    st.info("🆕 **Creating New Root Fill Run** (No Parent)")
                    active_parent_jids = {} # No dependencies
                    representative_parent = "NewRoot"
                else:
                    # Check compatibility of multiple parents
                    parent_infos = {rid: lineage[rid] for rid in selected_parents}
                    parent_types = set(info.get("simulation", "") for info in parent_infos.values())
                    
                    if len(parent_types) > 1:
                        st.error(f"⚠️ **Incompatible Types:** You have selected runs of different categories ({parent_types}). All selected runs must be either all 'Fill' or all 'Flow'.")
                        st.stop()
                    
                    parent_type = list(parent_types)[0]
                    if "Flow" in parent_type:
                        allowed_modes = ["flow_resume"]
                    else:
                        allowed_modes = ["fill_resume", "flow"]
                        
                    selected_mode = st.radio("Select Simulation Mode", allowed_modes, horizontal=True)
                    representative_parent = selected_parents[0]
                    
                    # Resolve dependencies for each selected parent
                    active_parent_jids = {}
                    for rid, info in parent_infos.items():
                        short_id = info['name'].replace('-', '_').replace('.', '_')
                        jid = node_to_jobid.get(short_id)
                        if jid:
                            active_parent_jids[rid] = jid
                    
                    if active_parent_jids:
                        st.info(f"🔗 **Lineage Dependencies:** {len(active_parent_jids)} of the selected runs are still active. Dependent jobs will wait for them.")
                    
                    # --- CONSOLIDATION / MERGE SECTION ---
                    if len(selected_parents) >= 2:
                        # 1. Identify Ultimate Parent (closest to root)
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
                        
                        # 2. Verify Lineage Continuity (all must descend from ultimate_parent)
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
                            
                            # Storage implication warning
                            st.sidebar.warning("⚠️ **Storage Tip**: Shortening the lineage to < 3 nodes may cause the Sync Manager to keep this data local rather than archiving it to the cloud, as it considers the branch 'actively being worked on'.")
                            
                            if st.sidebar.button(f"🚀 Merge {len(children_to_merge)} Runs", use_container_width=True):
                                try:
                                    from Pulse.pulse_core import SyncManager
                                    from Pulse.merge_runs import merge_folders
                                    
                                    with st.sidebar.status("Batch consolidating...", expanded=True) as status:
                                        for child_id in children_to_merge:
                                            child_name = lineage[child_id]['name']
                                            status.write(f"Processing {child_name}...")
                                            
                                            # Restore if needed
                                            if not os.path.exists(ultimate_parent):
                                                status.write(f"Restoring parent: {lineage[ultimate_parent]['name']}...")
                                                SyncManager.restore_run(ultimate_parent)
                                            if not os.path.exists(child_id):
                                                status.write(f"Restoring child: {child_name}...")
                                                SyncManager.restore_run(child_id)
                                            
                                            # Execute Merge
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
                    # Notes and Snapshot (Only for the first selected parent to avoid clutter)
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
                                    # Support manual refresh
                                    force_refresh = st.button("🔄 Refresh Snapshot", key=f"refresh_snap_{primary_parent}")
                                    snap_path = generate_snapshot(primary_parent, force=force_refresh)
                                    
                                    if snap_path:
                                        st.image(snap_path, caption=f"Last Snapshot of {parent_info['name']} (Y-Z Plane)", use_container_width=True)
                                    else:
                                        st.info("No snapshot available (no dump files found).")
                                except Exception as e:
                                    st.warning(f"Snapshot preview unavailable: {e}")



                    
                    st.divider()
                    
                
                # REACTIVE LAUNCH SECTION
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
                viscosity = gc8.number_input("Viscosity", value=0.001, format="%f")
                
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

                    fc7, fc8 = st.columns(2)
                    mode_params["mode"] = fc7.selectbox("Pouring Mode", ["2D_stacked", "2D_worst_case"], index=0)
                    mode_params["dump_file"] = fc8.text_input("Dump File Inc", value="simulation_templates/default_dump.inc")
                    mode_params["simulation"] = "Hopper_Fill"
                    mode_params["no-vtk"] = True

                    with st.expander("🛠️ Geometry Overrides", expanded=False):
                        st.info("Provide custom geometry parameters for the hopper(s) as a JSON object.")
                        geo_json = st.text_area("Geometry Variables (JSON)", value="{}", help='Example: {"orifice_half": 0.04, "hopper_ang": 55.0}. You can also provide lists for per-hopper settings: {"orifice_half": [0.05, 0.04]}')
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
                    amp_input = ffc2.text_input("Amplitude(s)", value="0.01")
                    mode_params["run_steps"] = ffc3.number_input("Run Steps", value=2000000, step=100000)
                    mode_params["osc_dir"] = ffc4.text_input("Oscillation Dir", value="z")
                    mode_params["_freq_list"] = [f.strip() for f in freq_input.split(",") if f.strip()]
                    mode_params["_amp_list"] = [a.strip() for a in amp_input.split(",") if a.strip()]
                    
                elif selected_mode == "flow_resume":
                    frc1, frc2 = st.columns(2)
                    mode_params["run_steps"] = frc1.number_input("Run Steps", value=1000000, step=100000)

                st.markdown("---")
                use_auto_pilot = st.checkbox("🤖 **Handover to Auto-Pilot**", value=False, help="Automatically queue and resume simulations until a target step count is reached.")
                
                ap_target, ap_inc = 10000000, 100000
                if use_auto_pilot:
                    ap_col1, ap_col2 = st.columns(2)
                    ap_target = ap_col1.number_input("Target Total Steps", value=10000000, step=1000000)
                    ap_inc = ap_col2.number_input("Step Increment", value=100000, step=10000)
                    st.info("💡 **Auto-Pilot Note**: Simulations will always run 'In-place' to ensure that each increment resumes from the absolute latest state of the previous run.")
                    inplace_mode = True
                else:
                    inplace_mode = False
                    if selected_mode != "fill":
                        inplace_mode = st.checkbox("😇 In-place Resumption", value=True, help="Run simulation in the same directory as parent. Recommended for long relaxations.")
                
                submit_btn = st.button("🚀 Launch Auto-Pilot" if use_auto_pilot else "🚀 Submit to PBS", use_container_width=True)
                    
                if submit_btn:
                    if mode_params.get("_invalid_geo"):
                        st.error("Submission blocked: Please fix the invalid JSON in Geometry Overrides.")
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
                            # Single job submission per parent
                            # Use unique random seeds if batching multiple parents
                            if len(selected_parents) > 1:
                                current_job_seed = random.randint(100000, 999999)
                            else:
                                current_job_seed = int(seed)
                                
                            new_seed_suffix = str(current_job_seed)[-6:]
                            if selected_mode == "fill_resume": prefix = "R"
                            elif selected_mode == "flow": prefix = "F"
                            elif selected_mode == "flow_resume": prefix = "FR"
                            else: prefix = "S"
                            
                            job_name = f"{prefix}{p_seed}_{new_seed_suffix}"
                            
                            # Finalize params
                            final_params = {
                                "num_procs": int(num_procs),
                                "num_threads": int(num_threads),
                                "dt": dt,
                                "viscosity": viscosity,
                                "seed": current_job_seed,
                                **mode_params
                            }
                            
                            # Inject parent-specific path
                            if selected_mode == "fill_resume":
                                final_params["restart_path"] = p_rid
                            elif selected_mode == "flow":
                                final_params["source_dir"] = p_rid
                                final_params["freq"] = float(mode_params.get("_freq_list", [0])[0])
                                final_params["amp"] = float(mode_params.get("_amp_list", [0])[0])
                            elif selected_mode == "flow_resume":
                                final_params["restart_path"] = p_rid
                            
                            final_params.pop("_freq_list", None)
                            final_params.pop("_amp_list", None)

                            # Final job parameters
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
                        # --- ENROLL IN AUTO-PILOT ---
                        with st.spinner("Enrolling goals..."):
                            auto_data = load_auto_pilot()
                            enrolled_count = 0
                            for p_rid in selected_parents:
                                path_key = p_rid
                                if p_rid == "NewRoot":
                                    import time
                                    path_key = f"NewRoot_{selected_mode}_{int(time.time())}"
                                
                                auto_data["goals"][path_key] = {
                                    "target_steps": ap_target,
                                    "increment": ap_inc,
                                    "in_place": inplace_mode,
                                    "last_submitted": None,
                                    "status": "Idle",
                                    "mode": selected_mode,
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
                            st.info("The system will start them during the next hourly cycle (or when you trigger a manual heartbeat).")
                            time.sleep(2)
                            st.rerun()
                    else:
                        # --- STANDARD PBS SUBMISSION ---
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
                                errf.write(f"\n--- Error at {time.strftime('%Y-%m-%d %H:%M:%S')} ---\n")
                                errf.write(traceback.format_exc())
                            st.error(f"Error submitting job: {e}. Check submit_error.log for details.")


    elif nav == "🤖 Auto-Pilot":
        st.header("🤖 Auto-Pilot Fleet Management")
        auto_data = load_auto_pilot()
        
        # --- Heartbeat Indicator ---
        hb_str = auto_data["settings"].get("last_heartbeat")
        if hb_str:
            last_hb = datetime.datetime.fromisoformat(hb_str)
            diff = datetime.datetime.now() - last_hb
            if diff.total_seconds() < 3900: # 65 minutes
                st.success(f"🟢 **SYSTEM ACTIVE** (Last heartbeat: {diff.total_seconds()/60:.1f}m ago)")
            else:
                st.error(f"🔴 **SYSTEM STALE** (Last heartbeat: {diff.total_seconds()/60:.1f}m ago)")
        else:
            st.warning("⚪ **SYSTEM INACTIVE** (No heartbeat recorded yet)")

        st.divider()
        
        # --- Manual Trigger ---
        ap_col1, ap_col2 = st.columns([2, 1])
        if AutoPilotManager.is_running():
            ap_col1.info("⚙️ Auto-Pilot Manager is currently running...")
            if ap_col2.button("🔄 Refresh Status", use_container_width=True):
                st.rerun()
        else:
            if ap_col1.button("🚀 Run Auto-Pilot Now", type="primary", use_container_width=True):
                success, msg = AutoPilotManager.trigger_manual()
                if success:
                    st.success(msg)
                    time.sleep(1)
                    st.rerun()
                else:
                    st.error(msg)
            
            if ap_col2.button("📄 View Log", use_container_width=True):
                log_path = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.log")
                if os.path.exists(log_path):
                    with open(log_path, "r") as f:
                        log_lines = f.readlines()
                        st.code("".join(log_lines[-50:]), language="text")
                else:
                    st.info("No log file found.")

        st.divider()
        
        # --- Global Settings ---
        with st.expander("⚙️ Auto-Pilot Logic Settings", expanded=False):
            col1, col2, col3 = st.columns(3)
            auto_data["settings"]["enabled"] = col1.toggle("Auto-Pilot Enabled", value=auto_data["settings"].get("enabled", True))
            auto_data["settings"]["max_concurrent"] = col2.number_input("Max Concurrent Jobs", value=auto_data["settings"].get("max_concurrent", 4), min_value=1)
            auto_data["settings"]["polite_mode"] = col3.toggle("Polite Mode (Yield to Students)", value=auto_data["settings"].get("polite_mode", True))
            if st.button("Save Logic Settings"):
                save_auto_pilot(auto_data)
                st.success("Logic settings saved.")
                st.rerun()

        # --- Scheduler Settings ---
        with st.expander("⏰ Scheduler Settings (Cron)", expanded=False):
            st.markdown("##### System Crontab Status")
            try:
                curr_cron = subprocess.check_output(["crontab", "-l"], text=True).strip()
                st.code(curr_cron, language="bash")
            except:
                st.info("No active crontab found for this user.")
            
            st.divider()
            st.markdown("##### Update Schedule")
            freq = st.selectbox("Select Frequency", 
                               ["Hourly (Recommended)", "Every 30 Minutes", "Every 2 Hours", "Every 6 Hours", "Daily (Midnight)"],
                               index=0)
            
            cron_map = {
                "Hourly (Recommended)": "0 * * * *",
                "Every 30 Minutes": "*/30 * * * *",
                "Every 2 Hours": "0 */2 * * *",
                "Every 6 Hours": "0 */6 * * *",
                "Daily (Midnight)": "0 0 * * *"
            }
            
            if st.button("Update System Schedule"):
                new_schedule = cron_map[freq]
                python_path = "/home/guest/miniconda3/envs/gchain/bin/python"
                manager_path = os.path.join(ROOT_DIR, "Pulse/auto_pilot_manager.py")
                log_path = os.path.join(ROOT_DIR, "Pulse/auto_pilot.log")
                
                cron_line = f"{new_schedule} cd {ROOT_DIR} && {python_path} {manager_path} >> {log_path} 2>&1"
                
                try:
                    # Clear old and add new
                    subprocess.run(f"(crontab -l 2>/dev/null | grep -v 'auto_pilot_manager.py'; echo '{cron_line}') | crontab -", shell=True, check=True)
                    st.success(f"Schedule updated to: {freq}")
                    time.sleep(1)
                    st.rerun()
                except Exception as e:
                    st.error(f"Failed to update crontab: {e}")

        # --- Fleet Progress ---
        st.subheader("📋 Active Fleet Progress")
        goals = auto_data.get("goals", {})
        if not goals:
            st.info("Your fleet is currently empty. Add nodes from the '🧬 Lineage' tab to start automated runs.")
        else:
            # Prepare table data
            # We need to read the latest lineage to show progress
            lineage = scan_dumping_yard() or {} # Refresh lineage for accuracy
            
            table_data = []
            for run_path, goal in goals.items():
                name = os.path.basename(run_path)
                current_info = lineage.get(run_path, {})
                current_steps = current_info.get("steps", 0)
                target = goal["target_steps"]
                progress = min(100, (current_steps / target * 100)) if target > 0 else 0
                
                # Get params from goal or lineage fallback
                params = goal.get("params", {})
                n_val = params.get("N")
                n_fill_val = params.get("n_fill")
                geo_vars = params.get("geometry_vars")
                
                # Fallback to lineage if not in goal params
                if n_val is None: n_val = current_info.get("N", "-")
                if n_fill_val is None: n_fill_val = current_info.get("params", {}).get("n_fill", "-")

                # Normalize to scalar for table display (prevents Arrow mixed-type errors)
                if isinstance(n_val, list) and len(n_val) > 0: n_val = n_val[0]
                if isinstance(n_fill_val, list) and len(n_fill_val) > 0: n_fill_val = n_fill_val[0]
                if not geo_vars: geo_vars = current_info.get("params", {}).get("geometry_vars", {})
                
                if isinstance(geo_vars, dict) and geo_vars:
                    geo_str = ", ".join([f"{k}:{v}" for k, v in geo_vars.items()])
                else:
                    geo_str = "-"

                table_data.append({
                    "Simulation": name,
                    "N": n_val,
                    "n_fill": n_fill_val,
                    "Geo": geo_str,
                    "Progress": f"{progress:.1f}%",
                    "Current": f"{current_steps:,}",
                    "Target": f"{target:,}",
                    "Increment": f"{goal['increment']:,}",
                    "Status": goal.get("status", "Idle"),
                    "Last Submit": goal.get("last_submitted", "Never")
                })
            
            st.table(table_data)
            
            if st.button("Clear Completed Goals"):
                new_goals = {k: v for k, v in goals.items() if lineage.get(k, {}).get("steps", 0) < v["target_steps"]}
                auto_data["goals"] = new_goals
                save_auto_pilot(auto_data)
                st.rerun()

    elif nav == "🔄 Sync":
        st.subheader("🔄 Sophisticated Google Drive Sync")
        st.markdown("""
        This module synchronizes your simulation data to Google Drive and manages your local storage.
        - **Safe Sync**: Automatically skips ongoing PBS jobs to avoid partial uploads.
        - **Smart Cleanup**: Deletes local files for synced runs that are NOT leaf nodes (parents).
        - **Lineage Integration**: Sync status is stored in `lineage.json` and visible in the graph.
        """)
        
        is_running = SyncManager.is_running()
        
        col1, col2 = st.columns([1, 1])
        
        if not is_running:
            c1, c2 = st.columns(2)
            if c1.button("🚀 Run in Dashboard (Thread)", use_container_width=True):
                success, msg = SyncManager.start_sync(user=user_filter)
                if success:
                    st.toast(msg)
                    st.rerun()
                else:
                    st.error(msg)
            
            if c2.button("💾 Submit as HPC Job (ppn=8)", use_container_width=True, type="primary"):
                success, msg = SyncManager.start_sync_pbs(user=user_filter)
                if success:
                    st.success(msg)
                    st.rerun()
                else:
                    st.error(msg)
        else:
            run_info = SyncManager.get_running_info()
            mode = run_info["mode"] if run_info else "Unknown"
            jid = run_info["id"] if run_info else "???"
            
            if mode == "HPC Job":
                st.warning(f"⚠️ **Sync is running as an HPC Job** (PBS ID: `{jid}`)")
            else:
                st.info(f"ℹ️ **Sync is running in Dashboard Thread** (PID: `{jid}`)")

            if col1.button("🛑 Stop Sync", use_container_width=True, type="secondary"):
                success, msg = SyncManager.stop_sync()
                if success:
                    st.toast(msg)
                    st.rerun()
                else:
                    st.error(msg)
        
        if col2.button("🔄 Refresh Logs", use_container_width=True):
            st.rerun()

        st.divider()
        
        # Live Log Section
        st.markdown("### 📝 Sync Logs & Progress")
        
        @st.fragment(run_every=2 if is_running else None)
        def render_sync_logs():
            logs = SyncManager.get_logs(max_lines=50)
            st.code(logs, language="bash")
            
            if SyncManager.is_running():
                st.info("🔄 Sync cycle is currently running in the background...")
                
                # 1. Detect Rclone Progress
                if "Transferred:" in logs:
                    try:
                        progress_match = re.search(r'(\d+)%,', logs)
                        if progress_match:
                            progress_val = int(progress_match.group(1))
                            st.progress(progress_val / 100.0, text=f"Syncing to Google Drive: {progress_val}%")
                    except:
                        pass
                # 2. Detect Compression Progress
                elif "PROGRESS_COMPRESS:" in logs:
                    try:
                        comp_match = re.findall(r'PROGRESS_COMPRESS: (\d+) / (\d+)', logs)
                        if comp_match:
                            current, total = map(int, comp_match[-1])
                            st.progress(current / total, text=f"Preparing Archives: {current}/{total} folders")
                    except:
                        pass
            else:
                st.success("✅ Sync is idle or completed.")
        
        render_sync_logs()

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
