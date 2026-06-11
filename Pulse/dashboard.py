import streamlit as st
import pandas as pd
import sys
import os
import json
import glob
from pathlib import Path

# Prevent "Errno 5" I/O errors by redirecting stdout/stderr to devnull
# This is required for libraries that try to flush standard streams in headless environments
try:
    if sys.stdout is not None:
        sys.stdout.flush()
except OSError:
    sys.stdout = open(os.devnull, 'w')
try:
    if sys.stderr is not None:
        sys.stderr.flush()
except OSError:
    sys.stderr = open(os.devnull, 'w')

import time
import subprocess
import psutil
import tempfile
import datetime

# Add project root to sys.path to allow importing from analysis module
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pulse_core import PBSManager, SimulationMonitor, SyncManager, AutoPilotManager
import streamlit.components.v1 as components
from analysis.controllers.sim_data import SimDataController
from analysis.controllers.renderer import SimulationRenderer
from analysis.controllers.highlighter import HighlightController
from stpyvista import stpyvista
import socket
import random

scratch_dir = os.path.join(ROOT_DIR, "scratch")
if scratch_dir not in sys.path:
    sys.path.append(scratch_dir)
import submit_flexible
from lineage_tracker import scan_dumping_yard
import importlib
import analysis.controllers.unified_renderer
importlib.reload(analysis.controllers.unified_renderer)
from analysis.controllers.unified_renderer import UnifiedRenderer
NOTES_FILE = os.path.join(ROOT_DIR, "Pulse/lineage_notes.json")
def load_notes():
    if os.path.exists(NOTES_FILE):
        try:
            with open(NOTES_FILE, "r") as f: return json.load(f)
        except: pass
    return {}

def save_note(run_id, note):
    notes = load_notes()
    notes[run_id] = note
    with open(NOTES_FILE, "w") as f: json.dump(notes, f)

VIZ_SETTINGS_FILE = os.path.join(ROOT_DIR, "Pulse/viz_settings.json")
def load_viz_settings():
    if os.path.exists(VIZ_SETTINGS_FILE):
        try:
            with open(VIZ_SETTINGS_FILE, "r") as f: return json.load(f)
        except: pass
    return {"offset": [0.0, 0.0, 0.0], "zoom": 1.0, "camera": None, "width": 1280, "height": 720}

def save_viz_settings(offset, zoom, camera=None, width=1280, height=720):
    try:
        with open(VIZ_SETTINGS_FILE, "w") as f:
            json.dump({"offset": offset, "zoom": zoom, "camera": camera, "width": width, "height": height}, f)
    except: pass

def wrap_sim_name(name, max_width=25):
    """Wraps simulation names by breaking on underscores to keep nodes compact."""
    if len(name) <= max_width: return name
    parts = name.split("_")
    lines = []
    current_line = ""
    for part in parts:
        if len(current_line) + len(part) + 1 > max_width and current_line:
            lines.append(current_line)
            current_line = part
        else:
            current_line = (current_line + "_" + part) if current_line else part
    if current_line: lines.append(current_line)
    return "<br/>".join(lines)

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

    # --- Navigation ---
    # We use a horizontal radio button to act as 'True Navigation'
    # Standard st.tabs execute all tabs' code on every rerun, which causes the 1-minute lag for the 3D visualiser.
    # This radio selector ensures ONLY the active page's code is executed.
    st.write('<style>div.row-widget.stRadio > div{flex-direction:row; justify-content: center; gap: 20px;} div.row-widget.stRadio label{background: #f0f2f6; padding: 10px 20px; border-radius: 5px; cursor: pointer;} div.row-widget.stRadio div[role="radiogroup"] > label[data-baseweb="radio"]{background: #f0f2f6; border: 1px solid #ddd;}</style>', unsafe_allow_html=True)
    
    # --- Navigation Persistence ---
    pages = ["📊 Active Queue", "🧬 Lineage", "🤖 Auto-Pilot", "🎬 Visualization", "📈 Analysis", "🧪 Results", "⏱️ ETA", "🕰️ History", "🔄 Sync"]
    
    # Initialize from URL or default
    query_nav = st.query_params.get("tab", pages[0])
    if query_nav not in pages:
        query_nav = pages[0]
    
    current_idx = pages.index(query_nav)
    
    nav = st.radio(
        "Select Page",
        pages,
        index=current_idx,
        horizontal=True,
        label_visibility="collapsed",
        key="nav_selector"
    )

    # Sync URL with selection
    if nav != query_nav:
        st.query_params["tab"] = nav

    st.divider()

    # --- AUTO-PILOT SETTINGS LOADER ---
    AUTO_PILOT_FILE = os.path.join(ROOT_DIR, "Pulse", "auto_pilot.json")
    def load_auto_pilot():
        if os.path.exists(AUTO_PILOT_FILE):
            try:
                with open(AUTO_PILOT_FILE, "r") as f: return json.load(f)
            except: pass
        return {"settings": {"enabled": False}, "goals": {}}

    def save_auto_pilot(data):
        with open(AUTO_PILOT_FILE, "w") as f:
            json.dump(data, f, indent=4)

    if nav == "📊 Active Queue":
        from Pulse.views.active_queue import render_active_queue
        render_active_queue(refresh_rate, user_filter)


    elif nav == "🕰️ History":
        from Pulse.views.history import render_history
        render_history(user_filter)


    elif nav == "📈 Analysis":
        from Pulse.views.analysis import render_analysis
        render_analysis()

    elif nav == "🧪 Results":
        from Pulse.views.results import render_results
        render_results(user_filter)

    elif nav == "⏱️ ETA":
        from Pulse.views.eta import render_eta
        render_eta()

    elif nav == "🎬 Visualization":
        from Pulse.views.visualization import render_visualization
        render_visualization()


    elif nav == "🧬 Lineage":
        from Pulse.views.lineage import render_lineage
        render_lineage(user_filter)


    elif nav == "🤖 Auto-Pilot":
        from Pulse.views.autopilot import render_autopilot
        render_autopilot(user_filter)


    elif nav == "🔄 Sync":
        from Pulse.views.sync import render_sync
        render_sync(user_filter)


    # 5. System Status (Master Node only)
    if "master" in socket.gethostname():

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
