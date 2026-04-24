import streamlit as st
import pandas as pd
from pulse_core import PBSManager
import time
import subprocess

st.set_page_config(page_title="Pulse Dashboard", page_icon="⚡", layout="wide")

st.title("⚡ Pulse Job Monitor")
st.markdown("Real-time PBS Cluster Dashboard")

# Sidebar for controls
st.sidebar.header("Controls")
refresh_rate = st.sidebar.slider("Refresh Rate (seconds)", 5, 60, 10)
user_filter = st.sidebar.text_input("User Filter", value="guest")

if st.sidebar.button("Refresh Now"):
    st.rerun()

# Fetch jobs
jobs = PBSManager.get_jobs(user=user_filter)

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

# 5. System Status (Master Node only)
if "master" in subprocess.getoutput("hostname"):
    st.divider()
    st.subheader("🖥️ Master Node System Status")
    
    col_gpu1, col_gpu2 = st.columns(2)
    
    try:
        gpu_info = subprocess.getoutput("nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits").split(",")
        if len(gpu_info) >= 2:
            col_gpu1.metric("GPU Utilization", f"{gpu_info[0].strip()}%")
            col_gpu2.metric("GPU Memory", f"{gpu_info[1].strip()} MiB")
    except:
        st.write("GPU monitoring unavailable")

# Auto-refresh
time.sleep(refresh_rate)
st.rerun()
