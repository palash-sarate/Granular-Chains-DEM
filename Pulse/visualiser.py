import streamlit as st
import pyvista as pv
import os
import glob
import nest_asyncio2
import json
import pandas as pd
from stpyvista import stpyvista
from analysis.controllers.unified_renderer import UnifiedRenderer
nest_asyncio2.apply()

# The dashboard.py defines ROOT_DIR
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_FILE = os.path.join(ROOT_DIR, ".pulse_viz_config.json")

def load_config():
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, 'r') as f:
                return json.load(f)
        except: pass
    return {}

def save_config(config):
    try:
        with open(CONFIG_FILE, 'w') as f:
            json.dump(config, f)
    except: pass

@st.cache_data(show_spinner=False)
def parse_lammps_dump(dump_path):
    """Parses a LAMMPS dump file and returns a pandas DataFrame."""
    skip = 0
    cols = []
    try:
        with open(dump_path, 'r') as f:
            for i, line in enumerate(f):
                if line.startswith("ITEM: ATOMS"):
                    cols = line.strip().split()[2:]
                    skip = i + 1
                    break
        
        # Optimization: Only load the columns we actually need for visualization
        needed_cols = [c for c in cols if c in ['x', 'y', 'z', 'mol', 'diameter']]
        df = pd.read_csv(dump_path, skiprows=skip, sep=r'\s+', names=cols, usecols=needed_cols, engine='c')
        return df
    except Exception as e:
        st.error(f"Error parsing dump file: {e}")
        return None

@st.cache_resource(show_spinner=False)
def load_vtk_geometry(file_path):
    """Loads a VTK geometry file with caching."""
    return pv.read(file_path)

def st_directory_picker(label, key, base_path):
    """A simple directory picker for Streamlit."""
    if key not in st.session_state:
        # Try to load from persistent config first
        config = load_config()
        last_dir = config.get(key, base_path)
        # Ensure it exists, else fallback to base
        if os.path.exists(last_dir):
            st.session_state[key] = last_dir
        else:
            st.session_state[key] = base_path
        
    curr = st.session_state[key]
    
    # Ensure path exists (e.g. if deleted externally)
    if not os.path.exists(curr):
        curr = base_path
        st.session_state[key] = base_path

    st.markdown(f"**{label}**")
    st.code(curr, language="bash")
    
    c1, c2, c3 = st.columns([1, 1, 3])
    if c1.button("⬆️ Up", key=f"{key}_up"):
        new_path = os.path.dirname(curr)
        st.session_state[key] = new_path
        # Save to persistent config
        config = load_config()
        config[key] = new_path
        save_config(config)
        st.rerun()
        
    if c2.button("🏠 Home", key=f"{key}_home"):
        st.session_state[key] = base_path
        # Save to persistent config
        config = load_config()
        config[key] = base_path
        save_config(config)
        st.rerun()
        
    try:
        subdirs = sorted([d for d in os.listdir(curr) if os.path.isdir(os.path.join(curr, d)) and not d.startswith(".")])
        if subdirs:
            chosen = st.selectbox("Browse subdirectories:", ["-- Select to enter --"] + subdirs, key=f"{key}_browse")
            if chosen != "-- Select to enter --":
                new_path = os.path.join(curr, chosen)
                st.session_state[key] = new_path
                # Save to persistent config
                config = load_config()
                config[key] = new_path
                save_config(config)
                st.rerun()
    except Exception as e:
        st.error(f"Access error: {e}")
        
    return st.session_state[key]

def set_cam(pos):
    st.session_state.cam_pos = pos

def render_visualiser():
    # Force off-screen rendering to prevent Errno 5 in headless/remote environments
    pv.OFF_SCREEN = True
    
    st.subheader("🎥 Advanced PyVista Visualizer")

    st.markdown("Visualize particles and geometry from simulation run directories using PyVista.")
    
    base_yard = os.path.join(ROOT_DIR, "dumping_yard")
    run_dir = st_directory_picker("Select Simulation Run Directory", "viz_run_dir", base_yard)
    
    if not run_dir or run_dir == base_yard:
        st.info("Please navigate into a simulation run directory (e.g., Grid_Hopper_Filling_test/Grid_Fill_...) to begin visualization.")
        return
        
    # Check for Geometry_vtk and chain directories
    geo_dir = os.path.join(run_dir, "Geometry_vtk")
    chain_dir = os.path.join(run_dir, "chain")
    
    has_geo = os.path.exists(geo_dir)
    has_chain = os.path.exists(chain_dir)
    
    if not has_geo and not has_chain:
        st.warning(f"No 'Geometry_vtk' or 'chain' folder found in selected directory.")
        return
        
    # Find vtk files
    vtk_files = []
    if has_geo:
        vtk_files = glob.glob(os.path.join(geo_dir, "*.vtk"))
        
    # Find dump files
    dump_files = []
    if has_chain:
        try:
            dump_files = sorted(glob.glob(os.path.join(chain_dir, "chain_*.dump")), key=lambda x: int(os.path.basename(x).replace("chain_", "").replace(".dump", "")))
        except ValueError:
            # Fallback sort if naming is not perfectly strictly chain_X.dump
            dump_files = sorted(glob.glob(os.path.join(chain_dir, "*.dump")))
        
    col1, col2 = st.columns([1, 3])
    
    with col1:
        st.markdown("### 🕹️ View Options")
        
        # Geometry options
        show_geom = st.checkbox("Show Geometry", value=True, disabled=not vtk_files)
        geom_alpha = st.slider("Geometry Opacity", 0.0, 1.0, 0.4)
        
        # Particle options
        st.divider()
        show_particles = st.checkbox("Show Particles", value=True, disabled=not dump_files)
        
        selected_dump = None
        if dump_files and show_particles:
            dump_names = [os.path.basename(d) for d in dump_files]
            
            # Reset index to latest only if directory changed or not initialized
            viz_idx_key = f"scrub_idx_{st.session_state.viz_run_dir}"
            if viz_idx_key not in st.session_state:
                st.session_state[viz_idx_key] = len(dump_names) - 1
            
            # Ensure index is still valid (in case files were deleted)
            if st.session_state[viz_idx_key] >= len(dump_names):
                st.session_state[viz_idx_key] = len(dump_names) - 1
                
            selected_idx = st.select_slider(
                "Timestep Scrub", 
                options=range(len(dump_names)), 
                key=viz_idx_key,
                format_func=lambda i: dump_names[i]
            )
            selected_dump = os.path.join(chain_dir, dump_names[selected_idx])
            
        st.divider()
        c1, c2 = st.columns(2)
        show_axes = c1.checkbox("Show Corner Axes", value=True)
        show_grid = c2.checkbox("Show 3D Grid", value=False)
        
        st.divider()
        st.markdown("### 📷 Camera Presets")
        if "cam_pos" not in st.session_state:
            st.session_state.cam_pos = "iso"
            
        # Camera buttons
        c_p1, c_p2 = st.columns(2)
        c_p1.button("ISO", on_click=set_cam, args=("iso",), width="stretch")
        c_p2.button("Reset View", on_click=set_cam, args=("reset",), width="stretch")
        
        c_x1, c_x2 = st.columns(2)
        c_x1.button("X+", on_click=set_cam, args=("yz",), width="stretch")
        c_x2.button("X-", on_click=set_cam, args=("zy",), width="stretch")
        
        c_y1, c_y2 = st.columns(2)
        c_y1.button("Y+", on_click=set_cam, args=("zx",), width="stretch")
        c_y2.button("Y-", on_click=set_cam, args=("xz",), width="stretch")
        
        c_z1, c_z2 = st.columns(2)
        c_z1.button("Z+", on_click=set_cam, args=("xy",), width="stretch")
        c_z2.button("Z-", on_click=set_cam, args=("yx",), width="stretch")
        
        # Add a flag to apply preset in the next render
        if "cam_pos" in st.session_state and st.session_state.cam_pos != "last":
            st.session_state.apply_cam_preset = True
        
        st.caption("Scene renders automatically on change.")

    with col2:
        # Auto-render logic: Always render if a directory is selected
        if True: 
            with st.spinner("Rendering Scene with PyVista..."):
                try:
                    # Setup Plotter using UnifiedRenderer
                    plotter = UnifiedRenderer.setup_plotter(window_size=[800, 600])
                    
                    UnifiedRenderer.apply_scene(
                        plotter,
                        vtk_files=vtk_files if show_geom else None,
                        dump_path=selected_dump if show_particles else None,
                        show_geometry=show_geom,
                        geom_opacity=geom_alpha,
                        show_particles=show_particles,
                        show_axes=show_axes,
                        show_grid=show_grid,
                        camera_state=None # Handled below by custom logic
                    )

                    # Camera logic: Apply after all meshes/axes are added to prevent VTK auto-reset
                    # 1. Apply preset if requested
                    if st.session_state.get("apply_cam_preset"):
                        cam = st.session_state.cam_pos
                        if cam == "iso":
                            plotter.view_isometric()
                        elif cam == "reset":
                            plotter.reset_camera()
                        else:
                            plotter.camera_position = cam
                        st.session_state.apply_cam_preset = False
                        st.session_state.cam_pos = "last"
                        st.session_state.last_cam_pos = plotter.camera_position
                    
                    # 2. Otherwise apply the last known camera position if it exists
                    elif "last_cam_pos" in st.session_state and st.session_state.last_cam_pos:
                        try:
                            plotter.camera_position = st.session_state.last_cam_pos
                        except: pass
                        
                    # Render using patched stpyvista
                    stpv_state = stpyvista(plotter, key="pv_plot")
                    
                    if stpv_state and "camera_position" in stpv_state:
                        st.session_state.last_cam_pos = stpv_state["camera_position"]
                    
                except Exception as e:
                    import traceback
                    st.error(f"Failed to render 3D scene: {e}")
                    with st.expander("Show Detailed Error Log"):
                        st.code(traceback.format_exc())

