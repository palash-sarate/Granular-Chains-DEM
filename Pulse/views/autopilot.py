import streamlit as st
import os
import datetime
import time
import subprocess
from Pulse.pulse_core import PBSManager, AutoPilotManager
from lineage_tracker import scan_dumping_yard

# File paths and helpers
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
NOTES_FILE = os.path.join(ROOT_DIR, "Pulse/lineage_notes.json")

def load_auto_pilot():
    if os.path.exists(AUTO_PILOT_FILE):
        try:
            import json
            with open(AUTO_PILOT_FILE, "r") as f: return json.load(f)
        except: pass
    return {"settings": {"enabled": False}, "goals": {}}

def save_auto_pilot(data):
    try:
        import json
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)
    except: pass

def load_notes():
    if os.path.exists(NOTES_FILE):
        try:
            import json
            with open(NOTES_FILE, "r") as f: return json.load(f)
        except: pass
    return {}

def save_note(run_id, note):
    notes = load_notes()
    notes[run_id] = note
    try:
        import json
        with open(NOTES_FILE, "w") as f: json.dump(notes, f)
    except: pass

def render_autopilot(user_filter: str):
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
                with open(log_path, "r", encoding="utf-8", errors="replace") as f:
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
                           ["Hourly (Recommended)", "Every 30 Minutes", "Every 15 Minutes", "Every 5 Minutes", "Every 2 Hours", "Every 6 Hours", "Daily (Midnight)"],
                           index=0)
        
        cron_map = {
            "Hourly (Recommended)": "0 * * * *",
            "Every 30 Minutes": "*/30 * * * *",
            "Every 15 Minutes": "*/15 * * * *",
            "Every 5 Minutes": "*/5 * * * *",
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
        # Prepare interactive layout
        lineage = scan_dumping_yard() or {} # Refresh lineage for accuracy
        
        # Fetch active jobs once for real-time status updates
        active_jobs = PBSManager.get_jobs(user=user_filter)
        
        # Interactive columns grid
        hdr_cols = st.columns([2.0, 1.5, 0.8, 1.0, 1.0, 1.2, 0.5])
        hdr_cols[0].markdown("**Simulation Run**")
        hdr_cols[1].markdown("**Step Progress**")
        hdr_cols[2].markdown("**Status**")
        hdr_cols[3].markdown("**Polite Mode**")
        hdr_cols[4].markdown("**Paused**")
        hdr_cols[5].markdown("**Last Submitted**")
        hdr_cols[6].markdown("**Actions**")
        st.markdown("<hr style='margin: 0px 0px 10px 0px; border-color: rgba(49, 51, 63, 0.2);'>", unsafe_allow_html=True)
        
        for run_path, goal in list(goals.items()):
            name = os.path.basename(run_path)
            current_info = lineage.get(run_path, {})
            current_steps = current_info.get("steps", 0)
            target = goal["target_steps"]
            progress = min(100.0, (current_steps / target * 100.0)) if target > 0 else 0.0
            
            # Retrieve parameters
            params = goal.get("params", {})
            n_val = params.get("N")
            n_fill_val = params.get("n_fill")
            geo_vars = params.get("geometry_vars")
            lepton_vars = params.get("lepton_vars")
            
            if n_val is None: n_val = current_info.get("N", "-")
            if n_fill_val is None: n_fill_val = current_info.get("params", {}).get("n_fill", "-")
            if isinstance(n_val, list) and len(n_val) > 0: n_val = n_val[0]
            if isinstance(n_fill_val, list) and len(n_fill_val) > 0: n_fill_val = n_fill_val[0]
            if not geo_vars: geo_vars = current_info.get("params", {}).get("geometry_vars", {})
            if not lepton_vars: lepton_vars = current_info.get("params", {}).get("lepton_vars", {})
            
            if isinstance(geo_vars, dict) and geo_vars:
                geo_str = ", ".join([f"{k}:{v}" for k, v in geo_vars.items()])
            else:
                geo_str = "-"

            if isinstance(lepton_vars, dict) and lepton_vars:
                lepton_str = ", ".join([f"{k}:{v}" for k, v in lepton_vars.items()])
            else:
                lepton_str = "-"
            
            # Determine mode/type icon
            mode = goal.get("mode", "fill_resume")
            if "flow" in mode:
                mode_icon = "🌊"
                mode_label = "Flow"
            elif mode == "calibration":
                mode_icon = "🔬"
                mode_label = "Calibration"
            else:
                mode_icon = "🧬"
                mode_label = "Fill"
            
            # Create row
            row_cols = st.columns([2.0, 1.5, 0.8, 1.0, 1.0, 1.2, 0.5])
            
            # Column 0: Simulation name and details
            detail_text = f"Type: {mode_label} | N: {n_val} | geo: {geo_str}"
            if lepton_str != "-":
                detail_text += f" | lepton: {lepton_str}"
            row_cols[0].markdown(
                f"**{mode_icon} {name}**<br>"
                f"<small style='color: grey;'>{detail_text}</small>", 
                unsafe_allow_html=True
            )
            
            # Column 1: Progress visual + details
            row_cols[1].markdown(f"**{progress:.1f}%** ({current_steps:,} / {target:,})")
            row_cols[1].progress(progress / 100.0)

            # Editable Target Steps
            try:
                safe_key = f"target_input_{abs(hash(run_path))}"
                new_target = row_cols[1].number_input(
                    "Target Steps",
                    value=int(target),
                    step=100000,
                    key=safe_key,
                    label_visibility="collapsed",
                )
            except Exception:
                new_target = target

            if int(new_target) != int(target):
                auto_data["goals"][run_path]["target_steps"] = int(new_target)
                save_auto_pilot(auto_data)
                st.toast(f"Updated target for {name} to {int(new_target):,}")
                time.sleep(0.5)
                st.rerun()
            
            # Determine live status dynamically
            target_seed = str(goal.get("params", {}).get("seed", ""))
            if not target_seed:
                import re
                match = re.search(r'S(\d+)', name)
                if match:
                    target_seed = match.group(1)
            
            if target_seed:
                expected_job_name = f"AP_{target_seed}_{name}"[:15]
            else:
                expected_job_name = f"AP_{name}"[:15]
            
            legacy_name = f"AP_{name}"
            
            is_job_running = False
            for j in active_jobs:
                j_name = j.get("Job_Name", "")
                if j_name == expected_job_name:
                    is_job_running = True
                    break
                if legacy_name == j_name or (len(j_name) == 15 and legacy_name.startswith(j_name)):
                    is_job_running = True
                    break

            is_paused = goal.get("paused", False)
            if is_paused:
                status = "Paused"
            elif current_steps >= target:
                status = "Completed"
            elif is_job_running:
                status = "Running"
            else:
                status = "Idle"
            
            # Column 2: Status
            if status == "Running":
                status_badge = "⚡ <span style='color: #00e676; font-weight: bold;'>Running</span>"
            elif status == "Completed":
                status_badge = "✅ <span style='color: #aeea00; font-weight: bold;'>Completed</span>"
            elif status == "Idle":
                status_badge = "💤 <span style='color: #80d8ff; font-weight: bold;'>Idle</span>"
            elif status == "Paused":
                status_badge = "⏸️ <span style='color: #ffb300; font-weight: bold;'>Paused</span>"
            else:
                status_badge = f"📋 <span>{status}</span>"
            row_cols[2].markdown(status_badge, unsafe_allow_html=True)
            
            # Column 3: Polite Mode toggle
            is_polite = goal.get("polite_mode", True)
            new_polite = row_cols[3].toggle(
                "Polite", 
                value=is_polite, 
                key=f"polite_toggle_{run_path}", 
                label_visibility="collapsed"
            )
            if new_polite != is_polite:
                auto_data["goals"][run_path]["polite_mode"] = new_polite
                save_auto_pilot(auto_data)
                st.toast(f"Updated polite mode for {name} to {'😇 Polite' if new_polite else '⚡ Priority'}")
                time.sleep(0.5)
                st.rerun()
            
            # Column 4: Pause/Resume toggle
            new_paused = row_cols[4].toggle(
                "Pause",
                value=is_paused,
                key=f"pause_toggle_{run_path}",
                label_visibility="collapsed"
            )
            if new_paused != is_paused:
                auto_data["goals"][run_path]["paused"] = new_paused
                if new_paused:
                    auto_data["goals"][run_path]["status"] = "Paused"
                else:
                    auto_data["goals"][run_path]["status"] = "Idle"
                save_auto_pilot(auto_data)
                st.toast(f"{'⏸️ Paused' if new_paused else '▶️ Resumed'} {name}")
                time.sleep(0.5)
                st.rerun()

            # Column 5: Last Submit time
            last_sub = goal.get("last_submitted")
            if last_sub:
                try:
                    dt_sub = datetime.datetime.fromisoformat(last_sub)
                    sub_str = dt_sub.strftime("%m/%d %H:%M")
                except:
                    sub_str = last_sub
            else:
                sub_str = "Never"
            row_cols[5].markdown(f"<div style='padding-top: 5px;'>{sub_str}</div>", unsafe_allow_html=True)
            
            # Column 6: Actions
            if row_cols[6].button("❌", key=f"rm_fleet_{run_path}", help="Remove from Auto-Pilot"):
                del auto_data["goals"][run_path]
                save_auto_pilot(auto_data)
                st.toast(f"Removed {name} from Auto-Pilot")
                time.sleep(0.5)
                st.rerun()
            
            # Inline snapshot and notes preview expander
            with st.expander(f"🖼️ Snapshot & Notes for {name}", expanded=False):
                col_note, col_snap = st.columns([1, 1])
                with col_note:
                    notes_db = load_notes()
                    current_note = notes_db.get(run_path, "")
                    new_note = st.text_area("🗒️ Run Notes", value=current_note, height=120, key=f"ap_note_{run_path}", help="Save observations or metadata for this run.")
                    if st.button("💾 Save Notes", key=f"ap_save_note_{run_path}"):
                        save_note(run_path, new_note)
                        st.toast("Notes saved!")
                with col_snap:
                    try:
                        from Pulse.snapshot_helper import generate_snapshot
                        force_refresh = st.button("🔄 Refresh Snapshot", key=f"ap_refresh_{run_path}")
                        snap_path = generate_snapshot(run_path, force=force_refresh)
                        if snap_path:
                            st.image(snap_path, caption="Latest Simulation State (Y-Z Plane)", use_container_width=True)
                        else:
                            st.info("No snapshot available (no dump files found).")
                    except Exception as e:
                        st.warning(f"Snapshot preview unavailable: {e}")

            st.markdown("<hr style='margin: 5px 0px 5px 0px; border-color: rgba(49, 51, 63, 0.1);'>", unsafe_allow_html=True)
        
        if st.button("Clear Completed Goals"):
            new_goals = {k: v for k, v in goals.items() if lineage.get(k, {}).get("steps", 0) < v["target_steps"]}
            auto_data["goals"] = new_goals
            save_auto_pilot(auto_data)
            st.rerun()
