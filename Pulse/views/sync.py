import streamlit as st
import re
from Pulse.pulse_core import SyncManager

def render_sync(user_filter: str):
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
