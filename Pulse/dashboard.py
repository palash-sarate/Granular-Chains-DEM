import streamlit as st
import pandas as pd
import sys
import os
import time
import subprocess
import psutil
import tempfile

# Add project root to sys.path to allow importing from analysis module
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from pulse_core import PBSManager, SimulationMonitor
import streamlit.components.v1 as components
from analysis.controllers.sim_data import SimDataController
from analysis.controllers.renderer import SimulationRenderer
from analysis.controllers.highlighter import HighlightController
import vedo

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
                subprocess.run(["python", "Pulse/lineage_tracker.py"])
                st.rerun()
        
        lineage = PBSManager.load_lineage()
        if not lineage:
            st.info("No lineage data found. Click 'Sync from Disk' to scan your simulations.")
        else:
            # 1. Prepare Mermaid Diagram
            mermaid_code = "graph LR\n"
            # Define nodes with styles
            for run_id, info in lineage.items():
                short_id = info['name'].replace('-', '_')
                node_label = f"{info['name']}<br/>(N={info['N']}, {info['steps']:,} steps)"
                
                # Style based on type
                style = ""
                if info["status"] == "Archived":
                    style = ":::archived"
                elif "Flow" in info["simulation"]:
                    style = ":::flow"
                else:
                    style = ":::active"
                
                mermaid_code += f'    {short_id}["{node_label}"]{style}\n'
            
            # Define relationships
            for run_id, info in lineage.items():
                if info["parent"] and info["parent"] in lineage:
                    parent_name = lineage[info["parent"]]['name'].replace('-', '_')
                    child_name = info['name'].replace('-', '_')
                    mermaid_code += f"    {parent_name} --> {child_name}\n"
            
            # Define styles
            mermaid_code += "    classDef active fill:#e1f5fe,stroke:#01579b,stroke-width:2px;\n"
            mermaid_code += "    classDef archived fill:#f5f5f5,stroke:#9e9e9e,stroke-dasharray: 5 5;\n"
            mermaid_code += "    classDef flow fill:#fff9c4,stroke:#fbc02d,stroke-width:2px;\n"

            # Render Mermaid
            components.html(
                f"""
                <div class="mermaid" style="background-color: white; padding: 20px; border-radius: 10px;">
                    {mermaid_code}
                </div>
                <script type="module">
                    import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
                    mermaid.initialize({{ startOnLoad: true, theme: 'neutral' }});
                </script>
                """,
                height=600,
                scrolling=True
            )
            
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

    # 5. System Status (Master Node only)
    if "master" in subprocess.getoutput("hostname"):
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
