import os
import sys
import json
import subprocess
import threading
import glob
import re
from datetime import datetime
from typing import List, Dict, Optional

# Add project root to sys.path
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from Pulse.pulse_core import PBSManager, SyncManager

class VizManager:
    VIZ_DIR = "Visualisations"
    VIZ_JOB_LOG = "Pulse/viz_job.log"
    VIZ_LOCK = "Pulse/viz.lock"

    CAMERA_PRESETS = {
        "ISO": {"pos_factor": [1, 1, 1], "up": [0, 0, 1]},
        "X":   {"pos_factor": [1, 0, 0], "up": [0, 0, 1]},
        "-X":  {"pos_factor": [-1, 0, 0], "up": [0, 0, 1]},
        "Y":   {"pos_factor": [0, 1, 0], "up": [0, 0, 1]},
        "-Y":  {"pos_factor": [0, -1, 0], "up": [0, 0, 1]},
        "Z":   {"pos_factor": [0, 0, 1], "up": [0, 1, 0]},
        "-Z":  {"pos_factor": [0, 0, -1], "up": [0, 1, 0]}
    }

    @staticmethod
    def get_descendants(start_path: str, lineage: Dict) -> List[str]:
        """Finds all nodes that are descendants of the start_path."""
        descendants = []
        
        # Build child map
        child_map = {}
        for path, info in lineage.items():
            parent = info.get('parent')
            if parent:
                if parent not in child_map: child_map[parent] = []
                child_map[parent].append(path)
        
        # BFS to find all descendants
        stack = [start_path]
        while stack:
            curr = stack.pop()
            children = child_map.get(curr, [])
            for child in children:
                if child not in descendants:
                    descendants.append(child)
                    stack.append(child)
        
        return sorted(descendants)

    @staticmethod
    def get_lineage_chain(start_path: str, end_path: str) -> List[str]:
        """Resolves the full path of runs between two nodes in the lineage."""
        lineage = PBSManager.load_lineage()
        if start_path not in lineage or end_path not in lineage:
            return []

        # We go backwards from end_path to start_path
        chain = []
        curr = end_path
        while curr:
            chain.append(curr)
            if curr == start_path:
                return list(reversed(chain))
            curr = lineage[curr].get('parent')
        
        return [] # No path found

    @staticmethod
    def get_chain_frames(chain_paths: List[str]) -> List[dict]:
        """Scans the chain and returns a list of all available frames as (path, local_ts)."""
        frames = []
        for path in chain_paths:
            if not os.path.exists(path):
                continue
            
            # Look in 'chain' subfolder first (standard), fallback to root
            search_dirs = [os.path.join(path, "chain"), path]
            found_dumps = []
            
            for s_dir in search_dirs:
                if os.path.isdir(s_dir):
                    found_dumps = glob.glob(os.path.join(s_dir, "*.dump"))
                    if found_dumps:
                        break
            
            ts_list = []
            for f in found_dumps:
                # Extract number from filename (e.g. particles_100.dump or 100.dump)
                match = re.search(r'(\d+)\.dump', os.path.basename(f))
                if match:
                    ts_list.append(int(match.group(1)))
            
            for ts in sorted(ts_list):
                frames.append({'path': path, 'ts': ts})
        return frames

    @staticmethod
    def apply_camera_preset(plt, preset_name: str, zoom: float = 1.0):
        """Applies a standard camera view to the plotter."""
        import numpy as np
        if preset_name not in VizManager.CAMERA_PRESETS:
            plt.reset_camera()
            return

        # 1. Get current focal point and bounds
        plt.reset_camera() # Initial fit
        cam = plt.camera
        fp = np.array(cam.GetFocalPoint())
        curr_pos = np.array(cam.GetPosition())
        dist = np.linalg.norm(curr_pos - fp)
        
        # 2. Calculate new position based on preset
        preset = VizManager.CAMERA_PRESETS[preset_name]
        factor = np.array(preset["pos_factor"])
        # If factor is [1,1,1] (ISO), we need to normalize it
        if preset_name == "ISO":
            factor = factor / np.sqrt(3)
        
        # Apply zoom by modifying distance
        dist = dist / zoom
        
        new_pos = fp + factor * dist
        cam.SetPosition(new_pos)
        cam.SetViewUp(preset["up"])
        cam.SetFocalPoint(fp)
        plt.render()

    @staticmethod
    def submit_viz_job(params: Dict):
        """Submits the movie generation process as a PBS job."""
        os.makedirs(os.path.join(ROOT_DIR, VizManager.VIZ_DIR), exist_ok=True)
        
        job_script = f"""#!/bin/bash
#PBS -N Viz_Movie
#PBS -q workq
#PBS -l nodes=master:ppn=16
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o {os.path.join(ROOT_DIR, VizManager.VIZ_JOB_LOG)}

cd $PBS_O_WORKDIR
trap 'kill 0' EXIT

# Activate conda environment
if [ -f /home/guest/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/guest/miniconda3/etc/profile.d/conda.sh
    conda activate gchain
fi

# Run the visualization script
python Pulse/viz_manager.py --params '{json.dumps(params)}'
"""
        
        script_path = os.path.join(ROOT_DIR, "Pulse/temp_viz.pbs")
        with open(script_path, "w") as f:
            f.write(job_script)
        
        try:
            result = subprocess.run(["qsub", script_path], capture_output=True, text=True, check=True)
            job_id = result.stdout.strip()
            
            # Write to lock file
            lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
            with open(lock_path, "w") as f:
                f.write(f"PBS:{job_id}")
            
            return True, f"Visualization job submitted: {job_id}"
        except Exception as e:
            return False, f"Failed to submit Viz job: {e}"
        finally:
            if os.path.exists(script_path):
                os.remove(script_path)

    @staticmethod
    def generate_preview(params: Dict) -> Optional[str]:
        """Generates a single frame preview image based on the first run in the chain."""
        import os
        # Force headless mode for VTK/vedo
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
        
        import vedo
        vedo.settings.enable_default_mouse_callbacks = False
        vedo.settings.enable_default_keyboard_callbacks = False
        
        import importlib
        import analysis.controllers.renderer
        import analysis.controllers.sim_data
        import analysis.controllers.vtk_overlay
        
        # Force reload to pick up fixes in these modules
        importlib.reload(analysis.controllers.renderer)
        importlib.reload(analysis.controllers.sim_data)
        importlib.reload(analysis.controllers.vtk_overlay)
        
        from analysis.controllers.renderer import SimulationRenderer
        from analysis.controllers.sim_data import SimDataController
        from analysis.controllers.vtk_overlay import VtkOverlayController

        chain_paths = params.get('chain_paths', [])
        if not chain_paths: return None
        
        # Target frame resolution
        target_path = params.get('target_path')
        target_ts = params.get('target_ts')
        
        if not target_path:
            target_path = chain_paths[0]
            
        # For preview, we only work if the folder exists locally
        if not os.path.exists(target_path):
            return f"Error: Data not local for {os.path.basename(target_path)}. Restore it first for a preview."

        camera_params = params.get('camera', {})
        vtk_files = params.get('vtk_files', [])
        show_geometry = params.get('show_geometry', True)
        timestamp_cfg = params.get('timestamp', {})
        dt = params.get('dt', 1e-6)
        
        preview_path = os.path.join(ROOT_DIR, VizManager.VIZ_DIR, "preview.png")
        os.makedirs(os.path.join(ROOT_DIR, VizManager.VIZ_DIR), exist_ok=True)

        try:
            plt = vedo.Plotter(offscreen=True, size=(1280, 720), interactive=False)
            
            data_ctrl = SimDataController()
            data_ctrl.stop_caching()
            data_ctrl.load_folder(target_path, enable_preloading=False)
            
            vtk_ctrl = VtkOverlayController(plt)
            for vtk_p in vtk_files:
                full_vtk = vtk_p if os.path.isabs(vtk_p) else os.path.join(target_path, vtk_p)
                if os.path.exists(full_vtk): vtk_ctrl.add_mesh(full_vtk)
                
            renderer = SimulationRenderer(plt, data_ctrl, vtk_ctrl, None, lambda: {'show_geometry': show_geometry})
            
            # 1. Render target timestep
            ts = target_ts if target_ts is not None else (data_ctrl.timesteps[0] if data_ctrl.timesteps else 0)
            
            # Explicitly load the batch containing the target timestep
            if data_ctrl.sim_source:
                batch_idx = data_ctrl.sim_source.get_batch_index_for_timestep(ts)
                data_ctrl.sim_source.load_batch(batch_idx)
                
            renderer.show_timestep(ts)
            
            # 2. Scene Configuration
            bg_choice = params.get('bg_color', 'Black').lower()
            plt.background(bg_choice)
            contrast_color = 'white' if bg_choice in ['black', 'dark gray'] else 'black'
            
            if params.get('show_axes'):
                plt.axes = 1
            
            sb_cfg = params.get('scalebar', {})
            if sb_cfg.get('show'):
                length = sb_cfg.get('length', 0.1)
                thick = sb_cfg.get('thickness', 0.02)
                
                # Compute bounds to position the bar
                acts = plt.actors
                if acts:
                    try:
                        b = vedo.utils.get_bounds(acts)
                        p_name = sb_cfg.get('pos', 'bottom-left')
                        
                        # X: Left vs Right
                        x = b[0] + length/2 if 'left' in p_name else b[1] - length/2
                        # Z: Bottom vs Top
                        z = b[4] if 'bottom' in p_name else b[5]
                        # Y: Always Front (ymin)
                        y = b[2] - thick
                        
                        pos = [x, y, z]
                        sb_block = vedo.Cube(pos=pos, side_length=(length, thick, thick)).c(contrast_color)
                        sb_label = vedo.Text3D(f"{length}", pos=[pos[0], pos[1]-thick*2, pos[2]], s=thick*1.5, c=contrast_color, justify='center')
                        plt.add(sb_block, sb_label)
                    except:
                        plt.add_scale_indicator(s=length, c=contrast_color)
                else:
                    plt.add_scale_indicator(s=length, c=contrast_color)

            # 3. Apply Camera Preset NOW that actors exist
            preset = params.get('camera_preset', 'ISO')
            zoom = params.get('zoom', 1.0)
            VizManager.apply_camera_preset(plt, preset, zoom=zoom)
            
            # Overlay
            ts_pos = timestamp_cfg.get('pos', 'top-left')
            ts_size = timestamp_cfg.get('size', 0.8)
            ts_actor = vedo.Text2D(f"PREVIEW\nTime: {ts*dt:.4e}s\nRun: {os.path.basename(target_path)}", pos=ts_pos, s=ts_size, c='red')
            plt.add(ts_actor)
            
            # Final render to ensure overlays are captured
            plt.render()
            
            plt.screenshot(preview_path)
            plt.close()
            return preview_path
        except Exception as e:
            return f"Preview error: {str(e)}"

    @staticmethod
    def run_rendering(params: Dict):
        """The actual rendering logic executed on the compute node."""
        import os
        # Force headless mode
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
        
        import vedo
        vedo.settings.enable_default_mouse_callbacks = False
        vedo.settings.enable_default_keyboard_callbacks = False
        
        import pandas as pd
        import importlib
        import analysis.controllers.renderer
        import analysis.controllers.sim_data
        import analysis.controllers.vtk_overlay

        # Force reload for fixes
        importlib.reload(analysis.controllers.renderer)
        importlib.reload(analysis.controllers.sim_data)
        importlib.reload(analysis.controllers.vtk_overlay)

        from analysis.controllers.renderer import SimulationRenderer
        from analysis.controllers.sim_data import SimDataController
        from analysis.controllers.vtk_overlay import VtkOverlayController
        
        # 1. Setup Parameters
        chain_paths = params.get('chain_paths', [])
        camera_params = params.get('camera', {}) # pos, fp, up, scale
        vtk_files = params.get('vtk_files', [])
        show_geometry = params.get('show_geometry', True)
        timestamp_cfg = params.get('timestamp', {})
        fps = params.get('fps', 24)
        output_name = params.get('output_name', f"movie_{datetime.now().strftime('%Y%m%d_%H%M%S')}.mp4")
        output_path = os.path.join(ROOT_DIR, VizManager.VIZ_DIR, output_name)

        # 2. Ensure Data is Local (Restore if needed)
        restored_paths = []
        for path in chain_paths:
            if not os.path.exists(path):
                print(f"Restoring missing data for: {os.path.basename(path)}")
                ok, msg = SyncManager.restore_run(path)
                if ok: restored_paths.append(path)
                else: print(f"Error restoring {path}: {msg}")

        # 3. Initialize Plotter (Offscreen)
        # Vedo requires a backend for offscreen rendering. VTK usually handles this.
        plt = vedo.Plotter(offscreen=True, size=(1920, 1080), interactive=False)
        video = vedo.Video(output_path, fps=fps)

        # 4. Setup Overlays
        ts_pos = timestamp_cfg.get('pos', 'top-left')
        ts_size = timestamp_cfg.get('size', 0.8)
        ts_actor = vedo.Text2D("", pos=ts_pos, s=ts_size, c='black', bg='white', alpha=0.7)
        plt.add(ts_actor)

        # 5. Rendering Loop
        try:
            for path in chain_paths:
                print(f"Processing run: {os.path.basename(path)}")
                
                # Initialize Data Controller
                data_ctrl = SimDataController()
                # Disable background caching for simplicity in batch mode
                data_ctrl.stop_caching() 
                
                # Load folder
                ok, err, _ = data_ctrl.load_folder(path, enable_preloading=False)
                if not ok:
                    print(f"Skipping {path}: {err}")
                    continue
                
                # Setup VTK Controller
                vtk_ctrl = VtkOverlayController(plt)
                for vtk_p in vtk_files:
                    full_vtk = vtk_p if os.path.isabs(vtk_p) else os.path.join(path, vtk_p)
                    if os.path.exists(full_vtk):
                        vtk_ctrl.add_mesh(full_vtk)
                
                # Setup Renderer
                renderer = SimulationRenderer(plt, data_ctrl, vtk_ctrl, None, lambda: {'show_geometry': show_geometry})
                
                # 1. Render first frame to establish bounds
                first_ts = data_ctrl.timesteps[0] if data_ctrl.timesteps else 0
                renderer.show_timestep(first_ts)
                
                # 2. Scene Configuration
                bg_choice = params.get('bg_color', 'Black').lower()
                plt.background(bg_choice)
                contrast_color = 'white' if bg_choice in ['black', 'dark gray'] else 'black'

                if params.get('show_axes'):
                    plt.axes = 1
                
                sb_cfg = params.get('scalebar', {})
                if sb_cfg.get('show'):
                    length = sb_cfg.get('length', 0.1)
                    thick = sb_cfg.get('thickness', 0.02)
                    acts = plt.actors
                    if acts:
                        try:
                            b = vedo.utils.get_bounds(acts)
                            p_name = sb_cfg.get('pos', 'bottom-left')
                            x = b[0] + length/2 if 'left' in p_name else b[1] - length/2
                            z = b[4] if 'bottom' in p_name else b[5]
                            y = b[2] - thick
                            
                            pos = [x, y, z]
                            sb_block = vedo.Cube(pos=pos, side_length=(length, thick, thick)).c(contrast_color)
                            sb_label = vedo.Text3D(f"{length}", pos=[pos[0], pos[1]-thick*2, pos[2]], s=thick*1.5, c=contrast_color, justify='center')
                            plt.add(sb_block, sb_label)
                        except:
                            plt.add_scale_indicator(s=length, c=contrast_color)
                    else:
                        plt.add_scale_indicator(s=length, c=contrast_color)

                # 3. Apply Camera Preset
                preset = params.get('camera_preset', 'ISO')
                zoom = params.get('zoom', 1.0)
                VizManager.apply_camera_preset(plt, preset, zoom=zoom)

                # 4. Timestep Overlay setup
                ts_pos = timestamp_cfg.get('pos', 'top-left')
                ts_size = timestamp_cfg.get('size', 0.8)
                ts_actor = vedo.Text2D("", pos=ts_pos, s=ts_size, c='white')
                plt.add(ts_actor)

                # Render Timesteps
                timesteps = data_ctrl.timesteps
                for ts in timesteps:
                    renderer.show_timestep(ts)
                    # Update timestamp
                    dt = params.get('dt', 1e-6)
                    ts_actor.text(f"Time: {ts*dt:.4e}s\nStep: {ts}\nRun: {os.path.basename(path)}")
                    plt.render()
                    video.add_frame()
            
            video.close()
            print(f"Movie saved successfully to: {output_path}")

        except Exception as e:
            print(f"Rendering failed: {e}")
            import traceback
            traceback.print_exc()

        finally:
            # 6. Cleanup (Free space for restored folders)
            if restored_paths:
                print(f"Cleaning up {len(restored_paths)} restored folders...")
                SyncManager.free_restored_space(restored_paths)
            
            # Remove lock
            lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
            if os.path.exists(lock_path):
                os.remove(lock_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", help="JSON string of parameters")
    args = parser.parse_args()
    
    if args.params:
        params = json.loads(args.params)
        VizManager.run_rendering(params)
