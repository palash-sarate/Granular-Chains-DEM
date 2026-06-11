import streamlit as st
import pandas as pd
from Pulse.pulse_core import PBSManager

def render_history(user_filter: str):
    st.subheader("🕰️ Comprehensive Job History")
    
    # 1. Bulk Scan Section
    with st.expander("🚀 Bulk Scan Job Range", expanded=False):
        col_start, col_end = st.columns(2)
        start_id = col_start.number_input("Start Job ID", value=3000, step=1)
        end_id = col_end.number_input("End Job ID", value=3050, step=1)
        
        if st.button("Start Bulk Scan", use_container_width=True):
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
        if st.button("Trace", use_container_width=True):
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
    
    # Load all from metadata for this user
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
            remaining_ids = set(edited_df["Job ID"].tolist())
            all_ids = set(df_hist["Job ID"].tolist())
            deleted_ids = all_ids - remaining_ids
            for did in deleted_ids:
                PBSManager.delete_cached_job(did)
            st.toast(f"Deleted {len(deleted_ids)} jobs from history.")
            st.rerun()
    else:
        st.info("No historical jobs found in metadata. Try running a Bulk Scan above!")
