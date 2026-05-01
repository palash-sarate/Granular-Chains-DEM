import streamlit as st
import pyvista as pv
import os
import glob
import nest_asyncio2
import json
import pandas as pd
from stpyvista import stpyvista
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
        
        if not cols:
            return None
            
        df = pd.read_csv(dump_path, skiprows=skip, sep=r'\s+', names=cols)
        return df
    except Exception as e:
        st.error(f"Error parsing dump file: {e}")
        return None

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
            # Use select_slider for a scrollbar-like scrubbing experience
            selected_idx = st.select_slider(
                "Timestep Scrub", 
                options=range(len(dump_names)), 
                value=len(dump_names)-1,
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
                    # Setup Plotter
                    plotter = pv.Plotter(window_size=[800, 600])
                    plotter.background_color = "white"
                    
                    # Add Geometry
                    if show_geom and vtk_files:
                        # Pleasant color palette for cycling
                        COLORS = ["#ADD8E6", "#90EE90", "#FFB6C1", "#FFFFE0", "#E6E6FA", "#F08080", "#AFEEEE", "#FFE4E1"]
                        for i, vtk_file in enumerate(vtk_files):
                            mesh = pv.read(vtk_file)
                            color = COLORS[i % len(COLORS)]
                            plotter.add_mesh(mesh, color=color, opacity=geom_alpha, show_edges=False)
                            
                    # Add Particles
                    if show_particles and selected_dump:
                        df = parse_lammps_dump(selected_dump)
                        if df is not None and not df.empty:
                            # Use a separate palette for particles to distinguish from geometry
                            P_COLORS = ["#FF5733", "#33FF57", "#3357FF", "#F333FF", "#FF33A8", "#33FFF3", "#F3FF33", "#FF8C00"]
                            
                            if all(k in df.columns for k in ['x', 'y', 'z']):
                                # Check if we can group by 'mol' (molecule ID / chain)
                                if 'mol' in df.columns:
                                    mols = df['mol'].unique()
                                    for i, m_id in enumerate(mols):
                                        m_df = df[df['mol'] == m_id]
                                        points = m_df[['x', 'y', 'z']].values
                                        pdata = pv.PolyData(points)
                                        
                                        if 'diameter' in m_df.columns:
                                            pdata['radius'] = m_df['diameter'].values / 2.0
                                        else:
                                            pdata['radius'] = [0.001] * len(points)
                                            
                                        sphere = pv.Sphere(radius=1.0, theta_resolution=8, phi_resolution=8)
                                        particles = pdata.glyph(scale="radius", geom=sphere, orient=False)
                                        
                                        color = P_COLORS[i % len(P_COLORS)]
                                        plotter.add_mesh(particles, color=color, smooth_shading=True)
                                else:
                                    # Fallback to single group if no 'mol' column
                                    points = df[['x', 'y', 'z']].values
                                    pdata = pv.PolyData(points)
                                    if 'diameter' in df.columns:
                                        pdata['radius'] = df['diameter'].values / 2.0
                                    else:
                                        pdata['radius'] = [0.001] * len(points)
                                    sphere = pv.Sphere(radius=1.0, theta_resolution=8, phi_resolution=8)
                                    particles = pdata.glyph(scale="radius", geom=sphere, orient=False)
                                    plotter.add_mesh(particles, color="#FF5733", smooth_shading=True)
                    
                    # (Camera application moved below axes/grid)
                        
                    if show_axes:
                        plotter.add_axes()
                        
                    if show_grid:
                        plotter.show_grid(color='gray', font_size=10)

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
                    
                    # 2. Otherwise apply the last known camera position if it exists
                    elif "last_cam_pos" in st.session_state and st.session_state.last_cam_pos:
                        try:
                            plotter.camera_position = st.session_state.last_cam_pos
                        except: pass
                        
                    # Render using patched stpyvista
                    # We use a unique key to keep the component state
                    stpv_state = stpyvista(plotter, key="pv_plot")
                    
                    # Save the camera position returned by the component for the next rerun
                    if stpv_state and "camera_position" in stpv_state:
                        st.session_state.last_cam_pos = stpv_state["camera_position"]
                    
                except Exception as e:
                    st.error(f"Failed to render 3D scene: {e}")
