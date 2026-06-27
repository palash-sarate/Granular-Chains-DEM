import streamlit as st
import pandas as pd
from Pulse.pulse_core import PBSManager

def render_active_queue(refresh_rate: int, user_filter: str):
    @st.fragment(run_every=refresh_rate)
    def active_queue_fragment():
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

            # Active Job Inspector (details and live log preview)
            import os
            st.write("---")
            st.subheader("🔍 Active Job Inspector")
            st.caption("Select a running/queued job to view resolved simulation/pipeline details and live output logs.")
            
            job_opts = [f"{j.get('id')} ({j.get('Job_Name')})" for j in jobs]
            selected_job_opt = st.selectbox("Select Job to Inspect", job_opts, key="active_q_inspect_job")
            
            if selected_job_opt:
                selected_job_id = selected_job_opt.split(" ")[0]
                sel_job = next((j for j in jobs if j.get("id") == selected_job_id), None)
                if sel_job:
                    lineage = PBSManager.load_lineage()
                    resolved = PBSManager.resolve_job_details(sel_job, lineage)
                    
                    st.markdown("##### Resolved Details")
                    col_j1, col_j2 = st.columns(2)
                    col_j1.markdown(f"**Run Name:** `{resolved['run_name']}`")
                    if resolved['run_path']:
                        col_j1.markdown(f"**Run Path:** `{resolved['run_path']}`")
                    col_j2.markdown(f"**Resolved Task/Stage:** `{resolved['stage']}`")
                    if resolved['log_path']:
                        col_j2.markdown(f"**Log File:** `{os.path.basename(resolved['log_path'])}`")
                    
                    if resolved['log_path']:
                        if os.path.exists(resolved['log_path']):
                            with st.expander("📄 View Live Output Log", expanded=True):
                                try:
                                    with open(resolved['log_path'], 'r', encoding='utf-8', errors='replace') as lf:
                                        log_content = lf.read()
                                        if log_content.strip():
                                            st.code(log_content[-20000:], language="text") # Show last 20k characters
                                        else:
                                            st.info("Log file is currently empty.")
                                except Exception as e:
                                    st.error(f"Error reading log file: {e}")
                        else:
                            st.info(f"Log file not created yet on disk at: `{resolved['log_path']}`")

        # Metrics summary
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Jobs", len(jobs))
        col2.metric("Running", len([j for j in jobs if j.get("job_state") == "R"]))
        col3.metric("Held/Queued", len([j for j in jobs if j.get("job_state") in ["H", "Q"]]))

    active_queue_fragment()
