import streamlit as st
import os
from Pulse.pulse_core import SimulationMonitor

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

def render_eta():
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
