import streamlit as st
import pandas as pd
import numpy as np
import os
import glob
import json
import plotly.express as px
from typing import List, Dict
from Pulse.pulse_core import PBSManager
from Pulse.pipeline.pipeline_engine import PipelineEngine

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def render_results(user_filter: str):
    st.subheader("🧪 Results Management & Pipeline Visualizer")
    st.markdown("Run run-level data processing pipelines, pool parameters, and visualize simulation results.")
    
    lineage = PBSManager.load_lineage()
    if not lineage:
        st.info("No lineage data found. Please scan your simulations in the 'Lineage' tab first.")
        return
        
    # Get active PBS jobs once for the status check
    active_jobs = PBSManager.get_jobs(user=user_filter)
    active_pbs_ids = {str(j["id"]) for j in active_jobs}
    
    tab_pipeline, tab_viz = st.tabs(["📦 Data Processing Pipelines", "📊 Results Visualization"])
    
    with tab_pipeline:
        st.write("### 🔍 Run-Level Pipeline Status")
        st.caption("Select a run to view its data pipeline stages, redo computations, or inspect PBS job logs.")
        
        all_paths = sorted(list(lineage.keys()))
        names = {p: lineage[p]['name'] for p in all_paths}
        selected_run_path = st.selectbox("Select Simulation Run", all_paths, format_func=lambda x: names.get(x, x), key="res_run_select")
        
        if selected_run_path:
            run_name = names[selected_run_path]
            st.markdown(f"#### 🏷️ Pipeline for: `{run_name}`")
            
            # Fetch pipeline status
            status_data = PipelineEngine.get_run_pipeline_status(selected_run_path, active_pbs_ids)
            
            # Draw pipeline stages dynamically based on RUN_PIPELINE registry
            stages = list(PipelineEngine.RUN_PIPELINE.keys())
            cols = st.columns(len(stages))
            
            for idx, stage in enumerate(stages):
                config = PipelineEngine.RUN_PIPELINE[stage]
                with cols[idx]:
                    with st.container(border=True):
                        st.markdown(f"##### Stage: `{stage}`")
                        st.caption(config["description"])
                        st.write(f"**Output:** `{config['output_file']}`")
                        
                        stg_status = status_data[stage]["status"]
                        job_id = status_data[stage]["job_id"]
                        
                        if stg_status == "completed":
                            st.success("🟢 Completed")
                        elif stg_status == "running":
                            st.info(f"🔵 Running (PBS: {job_id})")
                        elif stg_status == "pending":
                            st.warning(f"🟡 Pending (PBS: {job_id})")
                        elif stg_status == "failed":
                            st.error("🔴 Failed")
                            if status_data[stage].get("error"):
                                st.caption(f"⚠️ {status_data[stage]['error']}")
                        else:
                            st.markdown("⚪ **Not Started**")
                            
                        if st.button(f"🔄 Redo {stage}", key=f"redo_{stage}", use_container_width=True):
                            try:
                                PipelineEngine.trigger_pipeline(selected_run_path, stage)
                                st.toast(f"Submitted '{stage}' and dependent downstream jobs!")
                                st.rerun()
                            except Exception as e:
                                st.error(f"Failed to submit: {e}")
            
            # Job logs section
            pipeline_dir = os.path.join(selected_run_path, "results_pipeline")
            log_files = glob.glob(os.path.join(pipeline_dir, "job_*.log"))
            
            if log_files:
                st.write("---")
                st.markdown("##### 📄 Pipeline Job Logs")
                for log in sorted(log_files):
                    s_name = os.path.basename(log).replace("job_", "").replace(".log", "")
                    with st.expander(f"View Logs for {s_name}"):
                        try:
                            with open(log, 'r', encoding='utf-8', errors='replace') as f:
                                st.code(f.read(), language="text")
                        except Exception as log_err:
                            st.error(f"Could not read log file: {log_err}")
                            
        st.write("---")
        st.write("### 🌐 Global Pipeline Status")
        st.caption("Pool the computed results from all completed run-level dependencies.")
        
        global_status = PipelineEngine.get_global_pipeline_status(active_pbs_ids)
        stages_g = list(PipelineEngine.GLOBAL_PIPELINE.keys())
        
        for stage_g in stages_g:
            config_g = PipelineEngine.GLOBAL_PIPELINE[stage_g]
            g_st = global_status[stage_g]["status"]
            g_job_id = global_status[stage_g]["job_id"]
            
            col_g1, col_g2 = st.columns([3, 1])
            with col_g1:
                with st.container(border=True):
                    st.markdown(f"##### Global Stage: `{stage_g}`")
                    st.caption(config_g["description"])
                    st.write(f"**Output:** `{config_g['output_file']}`")
                    
                    if g_st == "completed":
                        st.success("🟢 Completed")
                    elif g_st == "running":
                        st.info(f"🔵 Running (PBS: {g_job_id})")
                    elif g_st == "pending":
                        st.warning(f"🟡 Pending (PBS: {g_job_id})")
                    elif g_st == "failed":
                        st.error("🔴 Failed")
                        if global_status[stage_g].get("error"):
                            st.caption(f"⚠️ {global_status[stage_g]['error']}")
                    else:
                        st.markdown("⚪ **Not Started**")
                        
            with col_g2:
                # Resolve completed runs matching dependencies
                completed_runs = []
                dep_stage = config_g["dependencies"][0]
                dep_output = PipelineEngine.RUN_PIPELINE[dep_stage]["output_file"]
                
                for path in all_paths:
                    dep_path = os.path.join(path, "results_pipeline", dep_output)
                    if os.path.exists(dep_path):
                        completed_runs.append(path)
                        
                if st.button(f"🔄 Redo {stage_g}", key=f"redo_{stage_g}", use_container_width=True, disabled=not completed_runs):
                    try:
                        PipelineEngine.submit_global_stage_job(stage_g, completed_runs)
                        st.toast(f"Submitted Global Pooling Job for {stage_g}!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Failed to submit: {e}")
                        
        global_log_path = os.path.join(ROOT_DIR, "Results_Pipeline", "job_pool_time_durations.log")
        if os.path.exists(global_log_path):
            with st.expander("View Global Pooling Job Logs"):
                try:
                    with open(global_log_path, 'r', encoding='utf-8', errors='replace') as f:
                        st.code(f.read(), language="text")
                except:
                    pass

    with tab_viz:
        st.write("### 📊 Parameter Sweeps & Mass Flow Rates")
        
        st.markdown("#### 1️⃣ Hopper Emptying Time Parameter Sweep")
        pooled_path = os.path.join(ROOT_DIR, "Results_Pipeline", "time_durations.parquet")
        
        if os.path.exists(pooled_path):
            try:
                df_pool = pd.read_parquet(pooled_path)
                if not df_pool.empty:
                    st.caption("Configure Plotly 2D Scatter parameters. Markers map different parameters dynamically.")
                    
                    param_options = ["N", "orifice_width", "frequency"]
                    
                    c1, c2, c3 = st.columns(3)
                    x_var = c1.selectbox("X-Axis Variable", param_options, index=1, key="sweep_x")
                    color_var = c2.selectbox("Color Category", param_options, index=2, key="sweep_color")
                    symbol_var = c3.selectbox("Marker Style Category", param_options, index=0, key="sweep_symbol")
                    
                    df_plot = df_pool.copy()
                    for col in param_options:
                        if col in df_plot.columns:
                            df_plot[col] = df_plot[col].astype(str)
                            
                    fig = px.scatter(
                        df_plot,
                        x=x_var,
                        y="time_taken",
                        color=color_var,
                        symbol=symbol_var,
                        title="Time Duration taken to Empty Hopper vs parameters",
                        labels={"time_taken": "Emptying Time (s)"},
                        hover_data=["run_name", "empty_timestep"]
                    )
                    fig.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))
                    fig.update_layout(xaxis_title=x_var.replace('_', ' ').title(), yaxis_title="Time Duration to Empty (s)")
                    
                    st.plotly_chart(fig, use_container_width=True)
                    st.dataframe(df_pool, use_container_width=True, hide_index=True)
                else:
                    st.info("The pooled results parquet is empty. Check pipeline jobs.")
            except Exception as viz_err:
                st.error(f"Error rendering sweep plot: {viz_err}")
        else:
            st.info("💡 **Global durations not pooled yet**: Run the global pooling stage inside 'Data Processing Pipelines' first.")
            
        st.divider()
        st.markdown("#### 2️⃣ Mass Flow Rate Curve Comparison")
        st.caption("Select multiple runs to compare their mass decay over time.")
        
        selected_compare_runs = st.multiselect("Select Runs to Compare", all_paths, format_func=lambda x: names.get(x, x), key="compare_runs_multi")
        
        if selected_compare_runs:
            comparison_records = []
            for p in selected_compare_runs:
                flow_path = os.path.join(p, "results_pipeline", "mass_flow_rate.parquet")
                if os.path.exists(flow_path):
                    try:
                        df_flow = pd.read_parquet(flow_path)
                        if not df_flow.empty:
                            run_info = lineage.get(p, {})
                            dt = float(run_info.get("params", {}).get("dt", 1e-6))
                            df_flow["time"] = df_flow["timestep"] * dt
                            df_flow["Run Name"] = names[p]
                            comparison_records.append(df_flow)
                    except Exception as e:
                        st.error(f"Error loading mass flow rate for {names[p]}: {e}")
                        
            if comparison_records:
                df_compare = pd.concat(comparison_records)
                fig_flow = px.line(
                    df_compare,
                    x="time",
                    y="mass_in_hopper",
                    color="Run Name",
                    title="Hopper Mass Decay Comparison",
                    labels={"time": "Simulation Time (s)", "mass_in_hopper": "Mass in Hopper (kg)"}
                )
                fig_flow.update_layout(xaxis_title="Simulation Time (s)", yaxis_title="Mass in Hopper (kg)")
                st.plotly_chart(fig_flow, use_container_width=True)
            else:
                st.warning("⚠️ None of the selected runs have completed the 'mass_flow_rate' stage yet.")
        else:
            st.info("💡 Select one or more runs from the dropdown above to render mass flow rate comparisons.")
