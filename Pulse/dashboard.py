import streamlit as st
import pandas as pd
from pulse_core import PBSManager
import time
import subprocess
import psutil

st.set_page_config(page_title="Pulse Dashboard", page_icon="⚡", layout="wide")

st.title("⚡ Pulse Job Monitor")
st.markdown("Real-time PBS Cluster Dashboard")

# Sidebar for controls
st.sidebar.header("Controls")
refresh_rate = st.sidebar.slider("Refresh Rate (seconds)", 5, 60, 10)
user_filter = st.sidebar.text_input("User Filter", value="guest")

if st.sidebar.button("Refresh Now"):
    st.rerun()

# Create Tabs for Active vs History
tab1, tab2 = st.tabs(["📊 Active Queue", "history 🕰️ Job History"])

with tab1:
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
                "Comment": job.get("comment", "")
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
        st.dataframe(df.style.map(color_state, subset=['State']), use_container_width=True)

        # Detailed Hold analysis
        held_jobs = [j for j in jobs if j.get("job_state") == "H"]
        if held_jobs:
            st.warning(f"Found {len(held_jobs)} jobs in HOLD state.")
            for hj in held_jobs:
                with st.expander(f"Hold Details: {hj.get('id')} ({hj.get('Job_Name')})"):
                    st.error(f"Reason: {hj.get('comment')}")
                    st.code(f"Error Path: {hj.get('Error_Path')}")

    # Metrics summary
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Jobs", len(jobs))
    col2.metric("Running", len([j for j in jobs if j.get("job_state") == "R"]))
    col3.metric("Held/Queued", len([j for j in jobs if j.get("job_state") in ["H", "Q"]]))

with tab2:
    st.subheader("🕰️ Persistent Job History")
    
    # 1. Bulk Scan Section
    with st.expander("🚀 Bulk Scan Job Range", expanded=False):
        col_start, col_end = st.columns(2)
        start_id = col_start.number_input("Start Job ID", value=3000, step=1)
        end_id = col_end.number_input("End Job ID", value=3050, step=1)
        
        if st.button("Start Bulk Scan"):
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
        if st.button("Trace"):
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
            use_container_width=True,
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

# 5. System Status (Master Node only)
if "master" in subprocess.getoutput("hostname"):
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

# Auto-refresh
time.sleep(refresh_rate)
st.rerun()
