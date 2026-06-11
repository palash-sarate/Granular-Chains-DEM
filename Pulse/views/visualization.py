import streamlit as st
import os
import glob
import datetime
import time
import subprocess
import pandas as pd
from Pulse.pulse_core import PBSManager, SyncManager
from Pulse.viz_manager import VizManager
from analysis.controllers.unified_renderer import UnifiedRenderer
from analysis.controllers.sim_data import SimDataController
from analysis.controllers.highlighter import HighlightController
from stpyvista import stpyvista

# Settings File Paths and Helpers
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
VIZ_SETTINGS_FILE = os.path.join(ROOT_DIR, "Pulse/viz_settings.json")

def load_viz_settings():
    if os.path.exists(VIZ_SETTINGS_FILE):
        try:
            import json
            with open(VIZ_SETTINGS_FILE, "r") as f: return json.load(f)
        except: pass
    return {
        "offset": [0.0, 0.0, 0.0],
        "zoom": 1.0,
        "camera": None,
        "width": 1280,
        "height": 720,
        "framing_mode": "manual",
        "viewport_bounds": [-0.2, 0.2, -0.1, 0.5],
        "text_font_size": 14,
        "text_x": 20,
        "text_y": 40,
        "axes_size": 0.10,
        "axes_x": 0.00,
        "axes_y": 0.00
    }

def save_viz_settings(offset, zoom, camera=None, width=1280, height=720, framing_mode="manual", viewport_bounds=None, text_font_size=14, text_x=20, text_y=40, axes_size=0.10, axes_x=0.00, axes_y=0.00):
    try:
        import json
        with open(VIZ_SETTINGS_FILE, "w") as f:
            json.dump({
                "offset": offset,
                "zoom": zoom,
                "camera": camera,
                "width": width,
                "height": height,
                "framing_mode": framing_mode,
                "viewport_bounds": viewport_bounds,
                "text_font_size": text_font_size,
                "text_x": text_x,
                "text_y": text_y,
                "axes_size": axes_size,
                "axes_x": axes_x,
                "axes_y": axes_y
            }, f)
    except: pass

def render_visualization():
    st.subheader("🎬 Lineage Visualization & Movie Generation")
    st.markdown("Create high-quality movies across multiple simulation stages with automatic cloud-restoration.")
    
    lineage = PBSManager.load_lineage()
    
    if not lineage:
        st.info("No lineage data found. Please scan your simulations in the 'Lineage' tab first.")
    else:
        # 1. Chain Selection
        st.write("### 1. Select Lineage Chain")
        use_manual = st.checkbox("🧩 Use Manual Path (e.g. for Archived Runs)", key="viz_use_manual")
        
        start_node = end_node = None # Initialize for safety

        if use_manual:
            manual_path = st.text_input("Absolute Path to Run Directory", 
                                       placeholder="/Data/palash_data/dumping_yard/Archived_Runs/...", key="viz_manual_path")
            if manual_path:
                manual_path = manual_path.strip()
                if os.path.exists(manual_path):
                    chain = [manual_path]
                    start_node = manual_path
                    end_node = manual_path
                    names = {manual_path: os.path.basename(manual_path)}
                    st.success(f"Manual Path Resolved: `{os.path.basename(manual_path)}`")
                else:
                    st.error("❌ Path does not exist. Please provide a valid absolute path.")
                    chain = []
            else:
                chain = []
        else:
            all_paths = sorted(list(lineage.keys()))
            names = {p: lineage[p]['name'] for p in all_paths}
            
            c1, c2 = st.columns(2)
            start_node = c1.selectbox("Start Node", all_paths, format_func=lambda x: names.get(x, x), index=0, key="viz_start_node")
            
            # Filter end_node options based on descendants of start_node (including start_node itself)
            possible_ends = [start_node] + VizManager.get_descendants(start_node, lineage)
            end_node = c2.selectbox("End Node", possible_ends, format_func=lambda x: names.get(x, x), index=len(possible_ends)-1, key="viz_end_node")

            if end_node:
                chain = VizManager.get_lineage_chain(start_node, end_node)
            else:
                chain = []
        
        if not chain:
            st.error("No valid direct lineage chain found. Select a different Start/End node.")
        else:
            st.success(f"Resolved Chain: {' → '.join([names[p] for p in chain])}")
            
            # Check for missing synced runs in the chain
            missing_synced = []
            for p in chain:
                if not os.path.exists(p) and lineage.get(p, {}).get("sync_status") == "Synced":
                    missing_synced.append(p)
            
            if missing_synced:
                st.warning(f"⚠️ **{len(missing_synced)} runs** in this chain are synced but not local. They must be restored for scrubbing and preview.")
                if st.button("📥 Restore Missing Data in Chain", use_container_width=True, type="primary"):
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    for i, p in enumerate(missing_synced):
                        status_text.text(f"Restoring {os.path.basename(p)}... ({i+1}/{len(missing_synced)})")
                        SyncManager.restore_run(p)
                        progress_bar.progress((i + 1) / len(missing_synced))
                    st.success("Restoration complete!")
                    st.session_state.pop("viz_chain_frames", None)
                    st.session_state.pop("last_viz_chain", None) # Force re-scan
                    st.rerun()
            
            # Frame Scanning & Scrubbing
            if "last_viz_chain" not in st.session_state or st.session_state["last_viz_chain"] != chain:
                st.session_state["last_viz_chain"] = chain
                with st.spinner("Scanning chain frames..."):
                    st.session_state["viz_chain_frames"] = VizManager.get_chain_frames(chain)
            
            target_frame = None
            all_frames = st.session_state.get("viz_chain_frames", [])
            
            if not all_frames:
                st.warning("No dump frames found in this chain. Ensure simulations have data files.")
            else:
                st.divider()
                st.write("### ⏱️ Timeline Scrubber")
                if len(all_frames) == 1:
                    scrub_idx = 0
                    st.info("ℹ️ Only 1 frame found in this chain.")
                else:
                    scrub_idx = st.slider("Scrub through all frames in chain", 0, len(all_frames)-1, 0, 
                                         format="Frame %d", help="Drag to select a specific moment from any run in the chain.")
                target_frame = all_frames[scrub_idx]
                viz_dt = st.session_state.get("viz_dt", 1e-6)
                st.caption(f"📍 **Selected**: `{os.path.basename(target_frame['path'])}` | Timestep: `{target_frame['ts']}` | `{target_frame['ts']*viz_dt:.3f}s`")

            # 2. Configuration
            st.divider()
            st.write("### 2. Rendering Configuration")
            
            col1, col2 = st.columns(2)
            
            with col1:
                saved_settings = load_viz_settings()
                
                # Get default bounds from metadata if available
                meta_bounds = {"y": [-0.2, 0.2], "z": [-0.1, 0.5]}
                if chain:
                    run_dir = chain[-1]
                    meta_names = ["grid_metadata.json", "metadata.json", "sim_metadata.json"]
                    for name in meta_names:
                        p = os.path.join(run_dir, name)
                        if os.path.exists(p):
                            try:
                                import json
                                with open(p, 'r') as f:
                                    meta = json.load(f)
                                env = meta.get("envelope", {})
                                tb = env.get("total_bounds", {})
                                if "y" in tb and "z" in tb:
                                    meta_bounds["y"] = [tb["y"][0], tb["y"][1]]
                                    meta_bounds["z"] = [tb["z"][0], tb["z"][1]]
                                    break
                            except:
                                pass

                # Initialize states
                if "viz_framing_mode" not in st.session_state:
                    st.session_state["viz_framing_mode"] = saved_settings.get("framing_mode", "manual")
                if "viz_viewport_bounds" not in st.session_state:
                    default_bounds = saved_settings.get("viewport_bounds")
                    if not default_bounds:
                        default_bounds = [
                            meta_bounds["y"][0], meta_bounds["y"][1],
                            meta_bounds["z"][0], meta_bounds["z"][1]
                        ]
                    st.session_state["viz_viewport_bounds"] = default_bounds

                framing_mode = st.radio(
                    "Framing Mode",
                    ["Manual Translation & Zoom", "Coordinate-Based Bounds (Y-Z Plane)"],
                    index=0 if st.session_state["viz_framing_mode"] == "manual" else 1,
                    horizontal=True,
                    help="Choose whether to position the camera using manual translations/zoom or exact physical coordinates."
                )
                framing_mode_val = "manual" if framing_mode == "Manual Translation & Zoom" else "bounds"

                if framing_mode_val == "bounds":
                    st.write("**Physical Simulation Bounds (Y-Z Plane)**")
                    b_y1, b_y2 = st.columns(2)
                    b_z1, b_z2 = st.columns(2)
                    
                    curr_bounds = st.session_state["viz_viewport_bounds"]
                    if len(curr_bounds) < 4:
                        curr_bounds = [
                            meta_bounds["y"][0], meta_bounds["y"][1],
                            meta_bounds["z"][0], meta_bounds["z"][1]
                        ]
                    
                    ymin = b_y1.number_input("Min Y (m)", value=float(curr_bounds[0]), format="%.4f")
                    ymax = b_y2.number_input("Max Y (m)", value=float(curr_bounds[1]), format="%.4f")
                    zmin = b_z1.number_input("Min Z (m)", value=float(curr_bounds[2]), format="%.4f")
                    zmax = b_z2.number_input("Max Z (m)", value=float(curr_bounds[3]), format="%.4f")
                    
                    viewport_bounds = [ymin, ymax, zmin, zmax]
                    st.session_state["viz_viewport_bounds"] = viewport_bounds
                    
                    st.write("**Canvas & Resolution**")
                    c_w, c_h = st.columns(2)
                    canvas_width = c_w.number_input("Width (px)", 320, 3840, saved_settings.get("width", 1280))
                    
                    phys_w = ymax - ymin
                    phys_h = zmax - zmin
                    if phys_w > 0 and phys_h > 0:
                        aspect_ratio = phys_w / phys_h
                        canvas_height = int(round(canvas_width / aspect_ratio))
                        canvas_height = max(240, min(3840, canvas_height))
                    else:
                        canvas_height = saved_settings.get("height", 720)
                    
                    c_h.number_input("Height (px) [Auto]", value=canvas_height, disabled=True, key="viz_canvas_height_disabled")
                    
                    new_offset = [0.0, 0.0, 0.0]
                    st.session_state["viz_scene_offset"] = new_offset
                    zoom = 1.0
                    st.session_state["viz_zoom"] = zoom
                else:
                    st.write("**Canvas & Resolution**")
                    c_w, c_h = st.columns(2)
                    canvas_width = c_w.number_input("Width (px)", 320, 3840, saved_settings.get("width", 1280))
                    canvas_height = c_h.number_input("Height (px)", 240, 3840, saved_settings.get("height", 720))
                    
                    st.write("**Scene Translation (Offset)**")
                    off_x, off_y, off_z = st.columns(3)
                    
                    if "viz_scene_offset" not in st.session_state:
                        st.session_state["viz_scene_offset"] = saved_settings.get("offset", [0.0, 0.0, 0.0])
                        st.session_state["viz_zoom"] = saved_settings.get("zoom", 1.0)
                        if saved_settings.get("camera"):
                            st.session_state["viz_last_cam_pos"] = saved_settings["camera"]

                    trans_x = off_x.number_input("Move X", value=st.session_state["viz_scene_offset"][0], format="%.3f")
                    trans_y = off_y.number_input("Move Y", value=st.session_state["viz_scene_offset"][1], format="%.3f")
                    trans_z = off_z.number_input("Move Z", value=st.session_state["viz_scene_offset"][2], format="%.3f")
                    
                    new_offset = [trans_x, trans_y, trans_z]
                    st.session_state["viz_scene_offset"] = new_offset

                    zoom = st.number_input("Base Zoom", min_value=0.01, max_value=100.0, value=st.session_state["viz_zoom"], format="%.2f", help="Set the camera magnification level.")
                    st.session_state["viz_zoom"] = zoom
                    viewport_bounds = None

                st.write("**Text Overlay Settings**")
                tx_col1, tx_col2, tx_col3 = st.columns(3)
                
                if "viz_text_font_size" not in st.session_state:
                    st.session_state["viz_text_font_size"] = saved_settings.get("text_font_size", 14)
                if "viz_text_x" not in st.session_state:
                    st.session_state["viz_text_x"] = saved_settings.get("text_x", 20)
                if "viz_text_y" not in st.session_state:
                    st.session_state["viz_text_y"] = saved_settings.get("text_y", 40)
                    
                text_font_size = tx_col1.number_input("Font Size", 6, 72, st.session_state["viz_text_font_size"])
                text_x = tx_col2.number_input("X Offset (px)", 0, 1000, st.session_state["viz_text_x"])
                text_y = tx_col3.number_input("Y Offset (px)", 0, 1000, st.session_state["viz_text_y"])
                
                st.session_state["viz_text_font_size"] = text_font_size
                st.session_state["viz_text_x"] = text_x
                st.session_state["viz_text_y"] = text_y

                st.write("**Corner Axes Settings**")
                ax_col1, ax_col2, ax_col3 = st.columns(3)
                
                if "viz_axes_size" not in st.session_state:
                    st.session_state["viz_axes_size"] = saved_settings.get("axes_size", 0.10)
                if "viz_axes_x" not in st.session_state:
                    st.session_state["viz_axes_x"] = saved_settings.get("axes_x", 0.00)
                if "viz_axes_y" not in st.session_state:
                    st.session_state["viz_axes_y"] = saved_settings.get("axes_y", 0.00)
                    
                axes_size = ax_col1.number_input("Axes Size", 0.01, 1.00, st.session_state["viz_axes_size"], step=0.01, format="%.2f", help="Axes scale fraction (e.g. 0.10 = 10% of window size)")
                axes_x = ax_col2.number_input("Axes X (left)", 0.00, 1.00, st.session_state["viz_axes_x"], step=0.01, format="%.2f", help="X starting position (e.g. 0.00 = left edge)")
                axes_y = ax_col3.number_input("Axes Y (bottom)", 0.00, 1.00, st.session_state["viz_axes_y"], step=0.01, format="%.2f", help="Y starting position (e.g. 0.00 = bottom edge)")
                
                st.session_state["viz_axes_size"] = axes_size
                st.session_state["viz_axes_x"] = axes_x
                st.session_state["viz_axes_y"] = axes_y

                st.write("**Time & Rate**")
                fps = st.number_input("FPS", 1, 60, 30)
                dt = st.number_input("Time Step (dt)", 1e-8, 1e-3, 1e-6, format="%.2e", key="viz_dt")

                st.session_state["viz_canvas_res"] = [canvas_width, canvas_height]
                st.session_state["viz_framing_mode"] = framing_mode_val
                
                # If anything changed, save to disk
                if framing_mode_val != saved_settings.get("framing_mode") or \
                   viewport_bounds != saved_settings.get("viewport_bounds") or \
                   (framing_mode_val == "manual" and (new_offset != saved_settings.get("offset") or zoom != saved_settings.get("zoom"))) or \
                   canvas_width != saved_settings.get("width") or \
                   canvas_height != saved_settings.get("height") or \
                   text_font_size != saved_settings.get("text_font_size") or \
                   text_x != saved_settings.get("text_x") or \
                   text_y != saved_settings.get("text_y") or \
                   axes_size != saved_settings.get("axes_size") or \
                   axes_x != saved_settings.get("axes_x") or \
                   axes_y != saved_settings.get("axes_y"):
                    
                    save_viz_settings(
                        new_offset, 
                        zoom, 
                        camera=st.session_state.get("viz_last_cam_pos"), 
                        width=canvas_width, 
                        height=canvas_height,
                        framing_mode=framing_mode_val,
                        viewport_bounds=viewport_bounds,
                        text_font_size=text_font_size,
                        text_x=text_x,
                        text_y=text_y,
                        axes_size=axes_size,
                        axes_x=axes_x,
                        axes_y=axes_y
                    )
            
            with col2:
                st.write("**Visual Elements**")
                
                # Discover all available VTKs in the chain
                all_chain_vtks = []
                for p in chain:
                    geo_dir = os.path.join(p, "Geometry_vtk")
                    if os.path.exists(geo_dir):
                        all_chain_vtks.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))
                
                vtk_basenames = sorted(list(set([os.path.basename(v) for v in all_chain_vtks])))
                selected_vtk_names = st.multiselect("Visible Geometry Layers", vtk_basenames, default=vtk_basenames)
                
                show_geo = st.checkbox("Show Geometry (Global)", value=True)
                show_interactive = st.checkbox("Enable Interactive 3D View", value=False)
                show_axes = st.checkbox("Show Corner Axes", value=True)
                show_grid = st.checkbox("Show 3D Grid", value=False)

            # Interactive 3D Configuration Preview
            st.divider()
            st.write("#### 🕹️ Interactive Configuration Preview")
            st.caption("Rotate, zoom, and pan to set the perfect camera view for your movie.")
            
            # Use the target_frame from the scrubber above
            target_dump = None
            if target_frame:
                f_path, f_ts = target_frame['path'], target_frame['ts']
                s_dirs = [os.path.join(f_path, "chain"), f_path]
                for sd in s_dirs:
                    if os.path.isdir(sd):
                        matches = glob.glob(os.path.join(sd, f"*{f_ts}.dump"))
                        if matches:
                            target_dump = matches[0]
                            break
            
            # Collect geometry files
            vtk_files = []
            for p in chain:
                geo_dir = os.path.join(p, "Geometry_vtk")
                if os.path.exists(geo_dir):
                    vtk_files.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))

            # Check for missing data in the interactive preview
            geo_vtk_dir = os.path.join(target_frame['path'], "Geometry_vtk") if target_frame else None
            missing_vtks = geo_vtk_dir and not os.path.exists(geo_vtk_dir)

            if target_frame and not os.path.exists(target_frame['path']):
                st.error(f"Data not local for {os.path.basename(target_frame['path'])}. Restore it first for a preview.")
                if st.button(f"📥 Restore {os.path.basename(target_frame['path'])} Now", use_container_width=True, type="primary"):
                    with st.spinner("Restoring..."):
                        success, msg = SyncManager.restore_run(target_frame['path'])
                        if success:
                            st.success("Restore complete!")
                            st.rerun()
                        else:
                            st.error(msg)
            elif target_frame:
                # --- Advanced Geometry Tools Expander ---
                with st.expander("🛠️ Advanced Geometry Tools", expanded=missing_vtks):
                    st.write("#### 🏗️ VTK Mesh (Re)generation")
                    st.caption("Reconstruct the hopper geometry from simulation metadata. Useful if files are missing or you want higher resolution.")
                    
                    v_col1, v_col2 = st.columns(2)
                    lattice_spacing = v_col1.number_input("Lattice Spacing (m)", 0.0001, 0.02, 0.002, format="%.4f", help="Resolution of the sampling grid. Lower = Sharper but Slower.")
                    recon_radius = v_col2.number_input("Recon. Radius (m)", 0.0, 0.1, 0.0, format="%.4f", help="Radius for surface reconstruction. Leave 0.0 for auto-calculation (1.2 * spacing).")
                    
                    if st.button("🔄 (Re)generate Geometry Meshes", use_container_width=True, type="primary" if missing_vtks else "secondary"):
                        with st.spinner("Regenerating geometry meshes..."):
                            try:
                                import warnings
                                warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
                                
                                from simulation.grid_hopper_manager import GridHopperManager
                                from simulation.runner import SimulationRunner
                                from analysis.geometry_extractor import GeometryExtractor
                                
                                # 1. Robust Metadata Search
                                meta_path = None
                                curr_p = Path(target_frame['path'])
                                for _ in range(4):
                                    for m_name in ["grid_metadata.json", "metadata.json"]:
                                        test_p = curr_p / m_name
                                        if test_p.exists():
                                            meta_path = test_p
                                            break
                                    if meta_path: break
                                    curr_p = curr_p.parent
                                    if curr_p.name == "dumping_yard": break

                                if not meta_path:
                                    st.error("Metadata not found. Cannot regenerate geometry.")
                                else:
                                    import json
                                    with open(meta_path, 'r') as f:
                                        meta = json.load(f)
                                    
                                    runner = SimulationRunner(lammps_executable="lmp") 
                                    mgr = GridHopperManager(runner)
                                    
                                    inc_name = meta.get("geometry_inc", "replicated_geometry.inc")
                                    inc_path = os.path.join(target_frame['path'], os.path.basename(inc_name))
                                    
                                    if not os.path.exists(inc_path):
                                        mgr._generate_replicated_geometry(
                                            setup_path=Path(meta["hopper_template_data"]),
                                            n_hoppers=len(meta["metadata"]),
                                            spacing=meta["spacing"],
                                            job_dir=Path(target_frame['path']),
                                            normalized_geo_vars=meta["geometry_vars"],
                                            inc_name=os.path.basename(inc_name)
                                        )
                                    
                                    extractor = GeometryExtractor(lammps_cmd="lmp")
                                    env = meta.get("envelope", {})
                                    bounds = [env['total_bounds']['x'][0], env['total_bounds']['x'][1],
                                              env['total_bounds']['y'][0], env['total_bounds']['y'][1],
                                              env['total_bounds']['z'][0], env['total_bounds']['z'][1]] if 'total_bounds' in env else None
                                    
                                    extractor.extract(
                                        inc_file=Path(inc_path),
                                        outdir=Path(geo_vtk_dir),
                                        auto_vis=True,
                                        combined=False,
                                        bounds=bounds,
                                        spacing=lattice_spacing,
                                        radius=recon_radius if recon_radius > 0 else None
                                    )
                                    st.success("Geometry VTKs generated successfully!")
                                    st.rerun()
                            except Exception as e:
                                st.error(f"Failed to generate VTKs: {e}")
                                import traceback
                                st.code(traceback.format_exc())

                if missing_vtks:
                    st.warning("⚠️ **Geometry VTKs Missing**: Hopper boundaries will not be visible. Use the tools above to generate them.")

                # Filter VTKs based on selection
                active_vtk_files = []
                for v in vtk_files:
                    if os.path.basename(v) in selected_vtk_names:
                        active_vtk_files.append(v)

                # Setup Plotter with custom resolution
                plotter = UnifiedRenderer.setup_plotter(window_size=st.session_state["viz_canvas_res"])

                UnifiedRenderer.apply_scene(
                    plotter,
                    vtk_files=active_vtk_files,
                    dump_path=target_dump,
                    show_geometry=show_geo,
                    show_particles=True,
                    show_axes=show_axes,
                    show_grid=show_grid,
                    camera_state=st.session_state.get("viz_last_cam_pos") if framing_mode_val == "manual" else None,
                    offset=st.session_state["viz_scene_offset"] if framing_mode_val == "manual" else None,
                    zoom=zoom if framing_mode_val == "manual" else 1.0,
                    viewport_bounds=viewport_bounds if framing_mode_val == "bounds" else None,
                    axes_viewport=(axes_x, axes_y, axes_x + axes_size, axes_y + axes_size)
                )
                
                # Viewfinder Border CSS
                vw, vh = st.session_state["viz_canvas_res"]
                st.markdown(f"""
                    <style>
                    .viz-viewfinder {{
                        border: 2px solid #333;
                        border-radius: 5px;
                        padding: 10px;
                        background: #000;
                        width: 100%;
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        margin-bottom: 10px;
                    }}
                    .viz-info-overlay {{
                        color: #666;
                        font-family: monospace;
                        font-size: 0.8rem;
                        text-align: center;
                        margin-top: 5px;
                    }}
                    </style>
                    <div class="viz-info-overlay">Recording Frame: {vw}x{vh}px (Aspect Ratio: {round(vw/vh, 2)})</div>
                """, unsafe_allow_html=True)

                # Display interactive plotter in a bordered 'Viewfinder' container (if enabled)
                if show_interactive:
                    with st.container(border=True):
                        st_state = stpyvista(plotter, key="viz_preview_plot")
                    
                    # Persist camera state on interaction
                    if st_state is not None:
                        # Extract camera position safely
                        new_pos = st_state.get("camera_position")
                        if new_pos:
                            # Subtract current offset before saving to maintain 'subject-relative' view
                            import numpy as np
                            off = np.array(st.session_state["viz_scene_offset"])
                            clean_pos = list(np.array(new_pos[0]) - off)
                            clean_fp = list(np.array(new_pos[1]) - off)
                            st.session_state["viz_last_cam_pos"] = [clean_pos, clean_fp, new_pos[2]]
                else:
                    st.info("💡 **Interactive Preview Disabled**: Use 'High-Res Snapshot' below to see your current scene setup.")
                
                c_cam1, c_cam2, c_cam3 = st.columns(3)
                if c_cam1.button("💾 Save View", use_container_width=True, help="Save the current camera angle as the default for this session and future restarts."):
                    save_viz_settings(
                        st.session_state["viz_scene_offset"], 
                        st.session_state["viz_zoom"], 
                        camera=st.session_state.get("viz_last_cam_pos"),
                        width=canvas_width,
                        height=canvas_height,
                        framing_mode=framing_mode_val,
                        viewport_bounds=viewport_bounds,
                        text_font_size=text_font_size,
                        text_x=text_x,
                        text_y=text_y,
                        axes_size=axes_size,
                        axes_x=axes_x,
                        axes_y=axes_y
                    )
                    st.success("Camera orientation and resolution saved to disk!")
                
                if c_cam2.button("🔄 Reset View", use_container_width=True):
                    st.session_state.pop("viz_last_cam_pos", None)
                    st.rerun()
                
                if c_cam3.button("📸 High-Res Snapshot", use_container_width=True):
                    with st.spinner("Generating High-Res PNG..."):
                        params = {
                            "chain_paths": chain,
                            "fps": fps,
                            "dt": dt,
                            "show_geometry": show_geo,
                            "selected_vtks": selected_vtk_names,
                            "show_axes": show_axes,
                            "show_grid": show_grid,
                            "camera_position": st.session_state.get("viz_last_cam_pos") if framing_mode_val == "manual" else None,
                            "offset": st.session_state["viz_scene_offset"] if framing_mode_val == "manual" else None,
                            "resolution": st.session_state["viz_canvas_res"],
                            "zoom": zoom if framing_mode_val == "manual" else 1.0,
                            "viewport_bounds": viewport_bounds if framing_mode_val == "bounds" else None,
                            "text_font_size": text_font_size,
                            "text_x": text_x,
                            "text_y": text_y,
                            "axes_viewport": (axes_x, axes_y, axes_x + axes_size, axes_y + axes_size)
                        }
                        shot_path = VizManager.generate_high_res_snapshot(params, target_frame)
                        st.session_state["last_viz_snapshot"] = shot_path

                if "last_viz_snapshot" in st.session_state:
                    st.write("---")
                    st.write("🖼️ **Latest High-Res Snapshot**")
                    st.image(st.session_state["last_viz_snapshot"], use_container_width=True)
                    if st.button("🗑️ Clear Snapshot"):
                        st.session_state.pop("last_viz_snapshot")
                        st.rerun()

                st.info("📌 **Capture Active**: Your current orientation, zoom, resolution, and offset will be used for the movie.")
            else:
                st.info("No frames found in the selected lineage. Please check if simulation data exists.")

            # 3. Job Submission
            st.divider()
            st.write("### 3. Execution")
            
            # Check for existing job
            lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
            is_running = False
            job_id = ""
            if os.path.exists(lock_path):
                try:
                    with open(lock_path, "r") as f:
                        content = f.read().strip()
                        if content.startswith("PBS:"):
                            job_id = content.replace("PBS:", "")
                            is_running = True
                except: pass

            if is_running:
                st.warning(f"⚠️ **Visualization job is currently running** (PBS ID: `{job_id}`)")
                
                # Log Display
                log_path = os.path.join(ROOT_DIR, VizManager.VIZ_JOB_LOG)
                if os.path.exists(log_path):
                    with st.expander("📄 View Visualization Logs", expanded=True):
                        with open(log_path, "r", encoding="utf-8", errors="replace") as f:
                            st.code(f.read(), language="text")
                        if st.button("🔄 Refresh Logs"):
                            st.rerun()

                if st.button("🛑 Cancel Viz Job", type="secondary"):
                    subprocess.run(["qdel", job_id])
                    if os.path.exists(lock_path): os.remove(lock_path)
                    st.rerun()
            else:
                output_name = st.text_input("Movie Name", value=f"{names[start_node]}_to_{names[end_node]}.mp4")
                
                if st.button("🎬 Submit Visualization Job", type="primary", use_container_width=True):
                    params = {
                        "chain_paths": chain,
                        "output_name": output_name,
                        "fps": fps,
                        "dt": dt,
                        "show_geometry": show_geo,
                        "selected_vtks": selected_vtk_names,
                        "show_axes": show_axes,
                        "show_grid": show_grid,
                        "camera_position": st.session_state.get("viz_last_cam_pos") if framing_mode_val == "manual" else None,
                        "offset": st.session_state["viz_scene_offset"] if framing_mode_val == "manual" else None,
                        "resolution": st.session_state["viz_canvas_res"],
                        "zoom": zoom if framing_mode_val == "manual" else 1.0,
                        "viewport_bounds": viewport_bounds if framing_mode_val == "bounds" else None,
                        "text_font_size": text_font_size,
                        "text_x": text_x,
                        "text_y": text_y,
                        "axes_viewport": (axes_x, axes_y, axes_x + axes_size, axes_y + axes_size)
                    }
                    
                    success, msg = VizManager.submit_viz_job(params)
                    if success:
                        st.success(msg)
                        st.rerun()
                    else:
                        st.error(msg)
            
            # Show Restore Button if needed
            if st.session_state.get("viz_restore_error"):
                st.error(st.session_state["viz_restore_error"])
                t_path = st.session_state.get("viz_restore_path")
                if t_path and st.button(f"📥 Restore {os.path.basename(t_path)} Now", use_container_width=True, type="primary"):
                    with st.spinner("Restoring from Google Drive..."):
                        success, msg = SyncManager.restore_run(t_path)
                        if success:
                            st.success("Restore complete! You can now preview the snapshot.")
                            st.session_state.pop("viz_restore_error", None)
                            st.rerun()
                        else:
                            st.error(msg)

            st.info("💡 Note: Missing simulation data will be automatically restored from Google Drive and cleaned up after the job finishes.")

            # Produced Visualisation Movies section
            st.divider()
            st.write("### 📥 Produced Visualisation Movies")
            
            viz_dir = os.path.join(ROOT_DIR, "Visualisations")
            if os.path.exists(viz_dir):
                movie_files = sorted(
                    [f for f in os.listdir(viz_dir) if f.endswith(".mp4")],
                    key=lambda x: os.path.getmtime(os.path.join(viz_dir, x)),
                    reverse=True
                )
                if movie_files:
                    import pandas as pd
                    
                    # 1. Scrollable List showing files metadata
                    st.write("##### 📁 Available Files")
                    df_data = []
                    for movie in movie_files:
                        movie_path = os.path.join(viz_dir, movie)
                        file_size_mb = os.path.getsize(movie_path) / (1024 * 1024)
                        created_time = datetime.datetime.fromtimestamp(os.path.getmtime(movie_path)).strftime("%Y-%m-%d %H:%M")
                        df_data.append({
                            "Movie Name": movie,
                            "Size (MB)": round(file_size_mb, 1),
                            "Created": created_time
                        })
                    
                    st.dataframe(
                        pd.DataFrame(df_data), 
                        use_container_width=True, 
                        hide_index=True,
                        height=200
                    )
                    
                    # 2. Scrollable Checklist for actions
                    if "selected_movies_history" not in st.session_state:
                        st.session_state["selected_movies_history"] = []
                        
                    # Sanitize history in case movies were deleted from disk
                    st.session_state["selected_movies_history"] = [
                        m for m in st.session_state["selected_movies_history"] if m in movie_files
                    ]
                    
                    selected_movies = []
                    st.write("##### 🗳️ Select Movie(s) to Play or Download")
                    with st.container(height=200):
                        for movie in movie_files:
                            if st.checkbox(movie, key=f"ap_movie_chk_{movie}"):
                                selected_movies.append(movie)
                    
                    # Update selection history to determine "latest selected"
                    added = [m for m in selected_movies if m not in st.session_state["selected_movies_history"]]
                    removed = [m for m in st.session_state["selected_movies_history"] if m not in selected_movies]
                    new_history = [m for m in st.session_state["selected_movies_history"] if m not in removed]
                    for m in added:
                        new_history.append(m)
                    st.session_state["selected_movies_history"] = new_history
                    
                    if selected_movies:
                        latest_movie = new_history[-1] if new_history else selected_movies[-1]
                        
                        st.write("---")
                        st.markdown(f"#### ⚙️ Actions for Selected ({len(selected_movies)} item(s))")
                        
                        # Play/Download columns
                        col_dl, col_del = st.columns(2)
                        
                        # Download action
                        if len(selected_movies) == 1:
                            movie_path = os.path.join(viz_dir, latest_movie)
                            try:
                                with open(movie_path, "rb") as f:
                                    movie_bytes = f.read()
                                col_dl.download_button(
                                    label=f"📥 Download {latest_movie}",
                                    data=movie_bytes,
                                    file_name=latest_movie,
                                    mime="video/mp4",
                                    use_container_width=True,
                                    key="dl_single_movie"
                                )
                            except Exception as e:
                                col_dl.error(f"Error reading file: {e}")
                        else:
                            # Multiselect ZIP download
                            import zipfile
                            import io
                            try:
                                zip_buffer = io.BytesIO()
                                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                                    for movie in selected_movies:
                                        movie_path = os.path.join(viz_dir, movie)
                                        zip_file.write(movie_path, arcname=movie)
                                
                                col_dl.download_button(
                                    label=f"📥 Download All Selected ({len(selected_movies)} files as ZIP)",
                                    data=zip_buffer.getvalue(),
                                    file_name="selected_visualisations.zip",
                                    mime="application/zip",
                                    use_container_width=True,
                                    key="dl_zip_movies"
                                )
                            except Exception as e:
                                col_dl.error(f"Error zipping files: {e}")
                                
                        # Delete action
                        if col_del.button("🗑️ Delete Selected File(s)", type="secondary", use_container_width=True, key="del_selected_btn"):
                            st.session_state["confirm_del_selected"] = True
                            st.rerun()
                            
                        if st.session_state.get("confirm_del_selected"):
                            st.warning(f"Are you sure you want to permanently delete {len(selected_movies)} selected movie(s) from disk?")
                            cy, cn = st.columns([1, 4])
                            if cy.button("Yes, Delete", type="primary", key="confirm_del_sel_yes"):
                                for movie in selected_movies:
                                    movie_path = os.path.join(viz_dir, movie)
                                    if os.path.exists(movie_path):
                                        os.remove(movie_path)
                                    # Clear checkbox state
                                    st.session_state.pop(f"ap_movie_chk_{movie}", None)
                                st.toast(f"Deleted {len(selected_movies)} movie(s).")
                                st.session_state.pop("confirm_del_selected", None)
                                st.session_state["selected_movies_history"] = []
                                time.sleep(0.5)
                                st.rerun()
                            if cn.button("Cancel", key="confirm_del_sel_no"):
                                st.session_state.pop("confirm_del_selected", None)
                                st.rerun()
                                
                        # Video Player for Latest Movie
                        st.write(f"##### 🎥 Previewing Latest Selected: `{latest_movie}`")
                        st.video(os.path.join(viz_dir, latest_movie))
                    else:
                        st.info("💡 Select one or more movies above to play or download them.")
                else:
                    st.info("No produced movies found in the Visualisations folder.")
            else:
                st.info("Visualisations folder not found.")
