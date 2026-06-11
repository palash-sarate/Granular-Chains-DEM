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

from analysis.controllers.unified_renderer import UnifiedRenderer

from Pulse.pulse_core import PBSManager, SyncManager

class VizManager:
    VIZ_DIR = "Visualisations"
    VIZ_JOB_LOG = "Pulse/viz_job.log"
    VIZ_LOCK = "Pulse/viz.lock"

    @staticmethod
    def get_descendants(start_path: str, lineage: Dict) -> List[str]:
        """Finds all nodes that are descendants of the start_path."""
        descendants = []
        child_map = {}
        for path, info in lineage.items():
            parent = info.get('parent')
            if parent:
                if parent not in child_map: child_map[parent] = []
                child_map[parent].append(path)
        
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

        chain = []
        curr = end_path
        while curr:
            chain.append(curr)
            if curr == start_path:
                return list(reversed(chain))
            curr = lineage[curr].get('parent')
        return []

    @staticmethod
    def get_chain_frames(chain_paths: List[str]) -> List[dict]:
        """Scans the chain and returns a list of all available frames."""
        frames = []
        for path in chain_paths:
            if not os.path.exists(path): continue
            
            search_dirs = [os.path.join(path, "chain"), path]
            found_dumps = []
            for s_dir in search_dirs:
                if os.path.isdir(s_dir):
                    found_dumps = glob.glob(os.path.join(s_dir, "*.dump"))
                    if found_dumps: break
            
            ts_list = []
            for f in found_dumps:
                match = re.search(r'(\d+)\.dump', os.path.basename(f))
                if match: ts_list.append(int(match.group(1)))
            
            for ts in sorted(ts_list):
                frames.append({'path': path, 'ts': ts})
        return frames

    @staticmethod
    def submit_viz_job(params: Dict):
        """Submits the movie generation process as a PBS job."""
        log_path = os.path.join(ROOT_DIR, VizManager.VIZ_JOB_LOG)
        if os.path.exists(log_path):
            try: os.remove(log_path)
            except: pass

        os.makedirs(os.path.join(ROOT_DIR, VizManager.VIZ_DIR), exist_ok=True)
        
        job_script = f"""#!/bin/bash
#PBS -N Viz_Movie
#PBS -q workq
#PBS -l nodes=master:ppn=16
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o {os.path.join(ROOT_DIR, VizManager.VIZ_JOB_LOG)}

cd $PBS_O_WORKDIR
export PYTHONPATH=$PYTHONPATH:$PBS_O_WORKDIR
trap 'kill 0' EXIT

if [ -f /home/guest/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/guest/miniconda3/etc/profile.d/conda.sh
    conda activate gchain
fi

python Pulse/viz_manager.py --params '{json.dumps(params)}'
"""
        script_path = os.path.join(ROOT_DIR, "Pulse/temp_viz.pbs")
        with open(script_path, "w") as f: f.write(job_script)
        
        try:
            result = subprocess.run(["qsub", script_path], capture_output=True, text=True, check=True)
            job_id = result.stdout.strip()
            lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
            with open(lock_path, "w") as f: f.write(f"PBS:{job_id}")
            return True, f"Visualization job submitted: {job_id}"
        except Exception as e:
            return False, f"Failed to submit Viz job: {e}"
        finally:
            if os.path.exists(script_path): os.remove(script_path)

    @staticmethod
    def generate_preview(params: Dict) -> Optional[str]:
        """Generates a single frame preview image using UnifiedRenderer."""
        import pyvista as pv
        chain_paths = params.get('chain_paths', [])
        if not chain_paths: return None
        
        target_path = params.get('target_path') or chain_paths[0]
        target_ts = params.get('target_ts')
        
        if not os.path.exists(target_path):
            return f"Error: Data not local for {os.path.basename(target_path)}. Restore it first for a preview."

        plotter = UnifiedRenderer.setup_plotter(off_screen=True)
        
        # Find target dump
        target_dump = None
        if target_ts is not None:
             dumps = glob.glob(os.path.join(target_path, "chain", f"*{target_ts}*.dump")) or \
                     glob.glob(os.path.join(target_path, f"*{target_ts}*.dump"))
             if dumps: target_dump = dumps[0]
        
        if not target_dump:
            dumps = sorted(glob.glob(os.path.join(target_path, "chain", "*.dump"))) or \
                    sorted(glob.glob(os.path.join(target_path, "*.dump")))
            if dumps: target_dump = dumps[-1]

        # Find vtk files
        vtk_files = []
        for p in chain_paths:
            geo_dir = os.path.join(p, "Geometry_vtk")
            if os.path.exists(geo_dir):
                vtk_files.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))

        UnifiedRenderer.apply_scene(
            plotter,
            vtk_files=vtk_files,
            dump_path=target_dump,
            show_geometry=params.get('show_geometry', True),
            show_particles=True,
            show_axes=params.get('show_axes', True),
            show_grid=params.get('show_grid', False),
            camera_state=params.get('camera_position'),
            offset=params.get('offset'),
            zoom=params.get('zoom', 1.0),
            viewport_bounds=params.get('viewport_bounds'),
            axes_viewport=params.get('axes_viewport')
        )

        preview_path = os.path.join(ROOT_DIR, VizManager.VIZ_DIR, "preview.png")
        os.makedirs(os.path.dirname(preview_path), exist_ok=True)
        plotter.screenshot(preview_path)
        plotter.close()
        return preview_path

    @staticmethod
    def generate_high_res_snapshot(params: Dict, target_frame: Dict):
        """Generates a high-resolution PNG snapshot using the movie rendering pipeline."""
        res = params.get('resolution', [1920, 1080])
        offset = params.get('offset', [0.0, 0.0, 0.0])
        zoom = params.get('zoom', 1.0)
        dt = params.get('dt', 1e-6)
        
        plotter = UnifiedRenderer.setup_plotter(off_screen=True, window_size=res)
        
        # Collect VTKs
        vtk_files = []
        selected_names = params.get('selected_vtks', [])
        chain_paths = params.get('chain_paths', [])
        for p in chain_paths:
            geo_dir = os.path.join(p, "Geometry_vtk")
            if os.path.exists(geo_dir):
                vtks_in_dir = glob.glob(os.path.join(geo_dir, "*.vtk"))
                for v in vtks_in_dir:
                    if not selected_names or os.path.basename(v) in selected_names:
                        vtk_files.append(v)
        
        # print(f"DEBUG Snapshot: Found {len(vtk_files)} VTKs for rendering")
        
        # Find dump file
        dump_file = None
        f_path, f_ts = target_frame['path'], target_frame['ts']
        s_dirs = [os.path.join(f_path, "chain"), f_path]
        for sd in s_dirs:
            if os.path.isdir(sd):
                matches = glob.glob(os.path.join(sd, f"*{f_ts}.dump"))
                if matches:
                    dump_file = matches[0]
                    break
        
        UnifiedRenderer.apply_scene(
            plotter,
            vtk_files=vtk_files,
            dump_path=dump_file,
            show_geometry=params.get('show_geometry', True),
            show_particles=True,
            show_axes=params.get('show_axes', True),
            show_grid=params.get('show_grid', False),
            camera_state=params.get('camera_position'),
            offset=offset,
            zoom=zoom,
            viewport_bounds=params.get('viewport_bounds'),
            axes_viewport=params.get('axes_viewport')
        )
        
        # Add Overlay (using exact pixel positioning to prevent scaling with resolution)
        t_font = params.get('text_font_size', 14)
        t_x = params.get('text_x', 20)
        t_y = params.get('text_y', 40)
        plotter.add_text(f"{f_ts*dt:.3f}s", position=(t_x, res[1] - t_y), font_size=t_font, color='red')
        
        snapshot_path = os.path.join(ROOT_DIR, VizManager.VIZ_DIR, f"snapshot_{f_ts}.png")
        os.makedirs(os.path.dirname(snapshot_path), exist_ok=True)
        plotter.screenshot(snapshot_path)
        plotter.close()
        return snapshot_path

    @staticmethod
    def run_rendering(params: Dict):
        """The actual rendering logic executed on the compute node."""
        import pyvista as pv
        chain_paths = params.get('chain_paths', [])
        output_name = params.get('output_name', "movie.mp4")
        fps = params.get('fps', 30)
        dt = params.get('dt', 1e-6)
        
        # Ensure data is local
        restored_paths = []
        for p in chain_paths:
            if not os.path.exists(p):
                ok, _ = SyncManager.restore_run(p)
                if ok: restored_paths.append(p)

        movie_dir = os.path.join(ROOT_DIR, VizManager.VIZ_DIR)
        os.makedirs(movie_dir, exist_ok=True)
        frames_dir = os.path.join(movie_dir, "temp_frames")
        if os.path.exists(frames_dir):
            import shutil
            shutil.rmtree(frames_dir)
        os.makedirs(frames_dir, exist_ok=True)
        
        frames = VizManager.get_chain_frames(chain_paths)
        if not frames: return

        vtk_files = []
        for p in chain_paths:
            geo_dir = os.path.join(p, "Geometry_vtk")
            if os.path.exists(geo_dir):
                vtk_files.extend(glob.glob(os.path.join(geo_dir, "*.vtk")))

        res = params.get('resolution', [1920, 1080])
        plotter = UnifiedRenderer.setup_plotter(off_screen=True, window_size=res)
        offset = params.get('offset', [0.0, 0.0, 0.0])
        
        print(f"Starting rendering of {len(frames)} frames...")
        selected_names = params.get('selected_vtks', [])
        
        for i, frame in enumerate(frames):
            plotter.clear()
            
            # Filter VTKs for this frame
            active_vtks = []
            if selected_names:
                for v in vtk_files:
                    if os.path.basename(v) in selected_names:
                        active_vtks.append(v)
            else:
                active_vtks = vtk_files

            f_path, f_ts = frame['path'], frame['ts']
            dump_file = None
            s_dirs = [os.path.join(f_path, "chain"), f_path]
            for sd in s_dirs:
                if os.path.isdir(sd):
                    matches = glob.glob(os.path.join(sd, f"*{f_ts}.dump"))
                    if matches:
                        dump_file = matches[0]
                        break
            if not dump_file: continue
            
            UnifiedRenderer.apply_scene(
                plotter,
                vtk_files=active_vtks,
                dump_path=dump_file,
                show_geometry=params.get('show_geometry', True),
                show_particles=True,
                show_axes=params.get('show_axes', True),
                show_grid=params.get('show_grid', False),
                camera_state=params.get('camera_position'),
                offset=offset,
                zoom=params.get('zoom', 1.0),
                viewport_bounds=params.get('viewport_bounds'),
                axes_viewport=params.get('axes_viewport')
            )
            # Use exact pixel positioning to prevent scaling with resolution
            t_font = params.get('text_font_size', 14)
            t_x = params.get('text_x', 20)
            t_y = params.get('text_y', 40)
            plotter.add_text(f"{f_ts*dt:.3f}s\nFrame: {i}", position=(t_x, res[1] - t_y - t_font), font_size=t_font, color='red')
            plotter.screenshot(os.path.join(frames_dir, f"frame_{i:06d}.png"))
            if i % 10 == 0: print(f"Rendered {i}/{len(frames)} frames...")

        plotter.close()
        
        movie_path = os.path.join(movie_dir, output_name)
        ffmpeg_cmd = ["ffmpeg", "-y", "-framerate", str(fps), "-i", os.path.join(frames_dir, "frame_%06d.png"), 
                      "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", movie_path]
        try:
            subprocess.run(ffmpeg_cmd, check=True)
            print(f"Movie saved: {movie_path}")
        except Exception as e:
            print(f"FFmpeg failed: {e}")
        
        # Cleanup
        if restored_paths:
            SyncManager.free_restored_space(restored_paths)
        lock_path = os.path.join(ROOT_DIR, VizManager.VIZ_LOCK)
        if os.path.exists(lock_path): os.remove(lock_path)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", type=str, required=True)
    args = parser.parse_args()
    VizManager.run_rendering(json.loads(args.params))
