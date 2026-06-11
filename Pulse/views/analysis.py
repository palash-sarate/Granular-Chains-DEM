import streamlit as st
import pandas as pd
import numpy as np
import os
import glob
from Pulse.pulse_core import PBSManager
from analysis.controllers.sim_data import SimDataController

def render_analysis():
    st.subheader("📈 Bond Angles & Lengths Analysis")
    st.markdown("Analyze the time-evolution of bond lengths and bond angles for all particle pairs in a run.")
    
    lineage = PBSManager.load_lineage()
    
    # Selection of run
    use_manual = st.checkbox("🧩 Use Manual Path (e.g. for Archived Runs)", key="anal_use_manual")
    
    run_path = None
    if use_manual:
        manual_path = st.text_input("Absolute Path to Run Directory", 
                                   placeholder="/Data/palash_data/dumping_yard/Archived_Runs/...", key="anal_manual_path")
        if manual_path:
            manual_path = manual_path.strip()
            if os.path.exists(manual_path):
                run_path = manual_path
                st.success(f"Manual Path Resolved: `{os.path.basename(manual_path)}`")
            else:
                st.error("❌ Path does not exist. Please provide a valid absolute path.")
    else:
        if not lineage:
            st.info("No lineage data found. Please scan your simulations in the 'Lineage' tab first, or use a manual path.")
        else:
            all_paths = sorted(list(lineage.keys()))
            names = {p: lineage[p]['name'] for p in all_paths}
            run_path = st.selectbox("Select Simulation Run", all_paths, format_func=lambda x: names.get(x, x), key="anal_run_path")

    if run_path:
        # Check files in run_path/bond and run_path/angle
        bond_dir = os.path.join(run_path, "bond")
        angle_dir = os.path.join(run_path, "angle")
        
        bond_files = glob.glob(os.path.join(bond_dir, "*.dump")) if os.path.isdir(bond_dir) else []
        angle_files = glob.glob(os.path.join(angle_dir, "*.dump")) if os.path.isdir(angle_dir) else []
        
        if not bond_files and not angle_files:
            st.warning("⚠️ No bond or angle dump files found in this run. Verify that the simulation configuration enabled local dumps.")
        else:
            # Initialize session state for analysis loading
            if "analysis_load_state" not in st.session_state:
                st.session_state.analysis_load_state = {
                    "run_path": None,
                    "active": False,
                    "current_batch": 0,
                    "total_batches": 0,
                    "decimation": 1,
                    "ctrl": None,
                    "df_bonds": pd.DataFrame(),
                    "df_angles": pd.DataFrame(),
                    "time_map": {},
                    "stopped": False
                }
            
            state = st.session_state.analysis_load_state
            
            # Slicing/Decimation factor selection
            decimation = st.sidebar.slider("Decimate/Subsample Data (Load every Nth step)", 1, 100, 1, 
                                           help="Higher values speed up parsing for large simulations by skipping frames.")
            
            # Reset state if run path or decimation changes
            if state["run_path"] != run_path or state["decimation"] != decimation:
                state["run_path"] = run_path
                state["decimation"] = decimation
                state["active"] = False
                state["current_batch"] = 0
                state["total_batches"] = 0
                state["ctrl"] = None
                state["df_bonds"] = pd.DataFrame()
                state["df_angles"] = pd.DataFrame()
                state["time_map"] = {}
                state["stopped"] = False
            
            # Load triggering button
            if not state["active"] and state["df_bonds"].empty and state["df_angles"].empty and not state["stopped"]:
                if st.button("🔄 Load Simulation Data", type="primary", use_container_width=True):
                    state["active"] = True
                    state["current_batch"] = 0
                    state["stopped"] = False
                    
                    ctrl = SimDataController()
                    success, err, queued = ctrl.load_folder(run_path, enable_preloading=False)
                    if not success:
                        st.error(f"Failed to scan folder: {err}")
                        state["active"] = False
                    else:
                        # Apply decimation if requested
                        if decimation > 1 and ctrl.sim_source:
                            ctrl.sim_source.timesteps = ctrl.sim_source.timesteps[::decimation]
                            ctrl.sim_source.atom_files = ctrl.sim_source.atom_files[::decimation]
                            ctrl.sim_source.bond_files = ctrl.sim_source.bond_files[::decimation]
                            ctrl.sim_source.angle_files = ctrl.sim_source.angle_files[::decimation]
                            ctrl.sim_source.loaded_batches.clear()
                            ctrl.timesteps = ctrl.sim_source.timesteps
                        
                        state["ctrl"] = ctrl
                        state["total_batches"] = ctrl.sim_source.get_batch_count() if ctrl.sim_source else 0
                        st.rerun()
            
            # Stepwise batch loading loop
            if state["active"]:
                if st.button("🛑 Stop Loading", key="btn_stop_loading", type="primary", use_container_width=True):
                    state["active"] = False
                    state["stopped"] = True
                    st.rerun()
                
                ctrl = state["ctrl"]
                total = state["total_batches"]
                
                # Dynamic WebSocket progress bar & status text container
                status_text = st.empty()
                progress_bar = st.progress(0)
                
                completed_successfully = True
                for curr in range(state["current_batch"], total):
                    status_text.info(f"⏳ Loading and parsing batch {curr + 1} of {total}...")
                    progress_bar.progress(curr / total)
                    
                    try:
                        ctrl.sim_source.load_batch(curr)
                        state["current_batch"] = curr + 1
                    except Exception as parse_err:
                        st.error(f"Error parsing batch {curr + 1}: {parse_err}")
                        completed_successfully = False
                        state["active"] = False
                        break
                
                if completed_successfully:
                    state["active"] = False
                    
                    df_b = ctrl.df_bonds
                    df_a = ctrl.df_angles
                    time_map = {ts: ctrl.sim_source.step_to_time.get(ts, ts) for ts in ctrl.timesteps} if ctrl.sim_source else {}
                    
                    if not df_b.empty:
                        df_b = df_b.reset_index()
                        df_b['time'] = df_b['timestep'].map(lambda x: time_map.get(x, x))
                    if not df_a.empty:
                        df_a = df_a.reset_index()
                        df_a['time'] = df_a['timestep'].map(lambda x: time_map.get(x, x))
                        
                    state["df_bonds"] = df_b
                    state["df_angles"] = df_a
                    state["time_map"] = time_map
                    st.success("✅ Data loaded successfully!")
                    st.rerun()
                else:
                    st.rerun()
            
            # If stopped, load whatever partial data has been parsed so far
            if state["stopped"] and state["df_bonds"].empty and state["df_angles"].empty and state["ctrl"] is not None:
                ctrl = state["ctrl"]
                df_b = ctrl.df_bonds
                df_a = ctrl.df_angles
                time_map = {ts: ctrl.sim_source.step_to_time.get(ts, ts) for ts in ctrl.timesteps} if ctrl.sim_source else {}
                
                if not df_b.empty:
                    df_b = df_b.reset_index()
                    df_b['time'] = df_b['timestep'].map(lambda x: time_map.get(x, x))
                if not df_a.empty:
                    df_a = df_a.reset_index()
                    df_a['time'] = df_a['timestep'].map(lambda x: time_map.get(x, x))
                    
                state["df_bonds"] = df_b
                state["df_angles"] = df_a
                state["time_map"] = time_map
            
            # Render control panel for resetting loaded data
            if not state["active"] and (not state["df_bonds"].empty or not state["df_angles"].empty or state["stopped"]):
                cols = st.columns([2, 1])
                with cols[0]:
                    if state["stopped"]:
                        st.warning("⚠️ Showing partial dataset (stopped by user).")
                    else:
                        st.success(f"📊 Dataset loaded successfully ({state['current_batch']} batches).")
                with cols[1]:
                    if st.button("🧹 Clear & Reset Data", use_container_width=True):
                        state["df_bonds"] = pd.DataFrame()
                        state["df_angles"] = pd.DataFrame()
                        state["time_map"] = {}
                        state["stopped"] = False
                        state["active"] = False
                        state["current_batch"] = 0
                        st.rerun()
            
            df_b = state["df_bonds"]
            df_a = state["df_angles"]
            
            if df_b.empty and df_a.empty:
                if not state["active"]:
                    st.info("ℹ️ Please click 'Load Simulation Data' above to retrieve and parse metrics.")
            else:
                tab_b, tab_a = st.tabs(["🔗 Bond Lengths", "📐 Bond Angles"])
                
                with tab_b:
                    st.subheader("🔗 Bond Lengths (Distance) Analysis")
                    if df_b is None or df_b.empty:
                        st.info("No bond length data available.")
                    else:
                        # 1. Global statistics over time
                        df_b_stats = df_b.groupby('time')['dist'].agg(['mean', 'std', 'min', 'max']).reset_index()
                        
                        st.markdown("#### Global Bond Lengths Over Time")
                        st.caption("Average bond length (with min/max boundaries) across all particle pairs in the system.")
                        st.line_chart(df_b_stats.set_index('time')[['mean', 'min', 'max']])
                        
                        # 2. Individual bond selection
                        st.markdown("#### Individual Bond Analysis")
                        bond_ids = sorted(df_b['id'].unique())
                        selected_bonds = st.multiselect("Select Bond IDs to trace", bond_ids, default=bond_ids[:min(5, len(bond_ids))])
                        
                        if selected_bonds:
                            df_b_sel = df_b[df_b['id'].isin(selected_bonds)]
                            df_b_pivot = df_b_sel.pivot(index='time', columns='id', values='dist')
                            st.line_chart(df_b_pivot)
                        else:
                            st.info("Please select one or more bond IDs above to view their individual trajectories.")
                            
                        # 3. Distribution at a specific time
                        st.markdown("#### Bond Length Distribution")
                        times = sorted(df_b['time'].unique())
                        selected_time = st.select_slider("Select Time for Distribution Hist", options=times, value=times[-1])
                        
                        df_b_t = df_b[df_b['time'] == selected_time]
                        counts, bin_edges = np.histogram(df_b_t['dist'].dropna(), bins=20)
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        hist_df = pd.DataFrame({'Bond Length': bin_centers, 'Count': counts})
                        st.bar_chart(hist_df.set_index('Bond Length'))
                
                with tab_a:
                    st.subheader("📐 Bond Angles (Theta) Analysis")
                    if df_a is None or df_a.empty:
                        st.info("No angle data available.")
                    else:
                        use_deg = st.checkbox("Show angles in Degrees (converted from Radians)", value=True, key="anal_use_deg")
                        angle_col = 'theta_deg' if use_deg else 'theta'
                        angle_unit = "Degrees" if use_deg else "Radians"
                        
                        if use_deg and 'theta_deg' not in df_a.columns:
                            df_a['theta_deg'] = df_a['theta'] * 180.0 / np.pi
                            
                        # 1. Global statistics over time
                        df_a_stats = df_a.groupby('time')[angle_col].agg(['mean', 'std', 'min', 'max']).reset_index()
                        
                        st.markdown(f"#### Global Bond Angles Over Time ({angle_unit})")
                        st.caption("Average bond angle (with min/max boundaries) across all particle triplets in the system.")
                        st.line_chart(df_a_stats.set_index('time')[['mean', 'min', 'max']])
                        
                        # 2. Individual angle selection
                        st.markdown("#### Individual Angle Analysis")
                        angle_ids = sorted(df_a['id'].unique())
                        selected_angles = st.multiselect("Select Angle IDs to trace", angle_ids, default=angle_ids[:min(5, len(angle_ids))])
                        
                        if selected_angles:
                            df_a_sel = df_a[df_a['id'].isin(selected_angles)]
                            df_a_pivot = df_a_sel.pivot(index='time', columns='id', values=angle_col)
                            st.line_chart(df_a_pivot)
                        else:
                            st.info("Please select one or more angle IDs above to view their individual trajectories.")
                            
                        # 3. Distribution at a specific time
                        st.markdown("#### Bond Angle Distribution")
                        times_a = sorted(df_a['time'].unique())
                        selected_time_a = st.select_slider("Select Time for Distribution Hist", options=times_a, value=times_a[-1], key="anal_time_a")
                        
                        df_a_t = df_a[df_a['time'] == selected_time_a]
                        counts, bin_edges = np.histogram(df_a_t[angle_col].dropna(), bins=20)
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        hist_df_a = pd.DataFrame({f'Bond Angle ({angle_unit})': bin_centers, 'Count': counts})
                        st.bar_chart(hist_df_a.set_index(f'Bond Angle ({angle_unit})'))
