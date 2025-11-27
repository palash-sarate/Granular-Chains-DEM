import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import pandas as pd
import time
import shutil

class Animator:
    def __init__(self, dataframe, output_file="simulation.mp4"):
        """
        Initialize the Animator with a pandas DataFrame.
        
        :param dataframe: Pandas DataFrame containing simulation data.
        :param output_file: Name of the output video file.
        """
        # Create a copy to avoid modifying the original data
        self.df = dataframe.copy()
        
        # Check if we have a MultiIndex (likely timestep, id) and reset it to make them columns
        if isinstance(self.df.index, pd.MultiIndex):
            self.df.reset_index(inplace=True)
        # Check if 'timestep' or 'id' are in index but not columns (single index case)
        elif self.df.index.name in ['timestep', 'id']:
             self.df.reset_index(inplace=True)

        self.output_file = output_file
        
        # Validate required columns
        required_cols = ['timestep', 'id', 'x', 'y', 'z']
        missing = [col for col in required_cols if col not in self.df.columns]
        if missing:
            raise ValueError(f"DataFrame must contain columns: {required_cols}. Missing: {missing}")

        # Pre-calculate global bounds for consistent axis scaling
        self.x_limits = self._safe_limits(self.df['x'].min(), self.df['x'].max())
        self.y_limits = self._safe_limits(self.df['y'].min(), self.df['y'].max())
        self.z_limits = self._safe_limits(self.df['z'].min(), self.df['z'].max())

    def _safe_limits(self, vmin, vmax, padding=0.5):
        """
        Ensures limits are not identical. 
        If min == max, expands the range by +/- padding.
        """
        if vmin == vmax:
            return (vmin - padding, vmax + padding)
        return (vmin, vmax)
    
    @staticmethod
    def _render_single_frame(frame_idx, timestep, group_data, limits, color_by, cmap, resolution_dpi, point_size, temp_dir, view_angles):
        """
        Static method to render a single frame. 
        Must be static to be picklable for ProcessPoolExecutor.
        """
        if len(view_angles) == 3:
            elev, azim, roll = view_angles
        else:
            elev, azim = view_angles
            roll = 0
        
        # Extract data
        x = group_data['x'].values
        y = group_data['y'].values
        z = group_data['z'].values
        
        # Determine colors
        if color_by in group_data.columns:
            c_values = group_data[color_by].values
        elif color_by == 'velocity_magnitude':
            if all(c in group_data.columns for c in ['vx', 'vy', 'vz']):
                vx, vy, vz = group_data['vx'].values, group_data['vy'].values, group_data['vz'].values
                c_values = np.sqrt(vx**2 + vy**2 + vz**2)
            else:
                c_values = np.zeros_like(x)
        else:
            # Fallback if color column not found
            c_values = np.zeros_like(x)

        # Determine marker sizes (render as spheres with diameter)
        if 'diameter' in group_data.columns:
            diameters = group_data['diameter'].values
            
            # Calculate scaling factor to map data units to point sizes
            # Figure size is 10x10 inches. 3D axis typically takes ~70% of space.
            # 1 inch = 72 points.
            # Scale = (Figure Size in Points * Fraction) / Axis Data Range
            
            x_range = limits[0][1] - limits[0][0]
            y_range = limits[1][1] - limits[1][0]
            z_range = limits[2][1] - limits[2][0]
            max_range = max(x_range, y_range, z_range)
            if max_range == 0: max_range = 1.0
            
            # Approximate scaling: 10 inches * 72 points/inch * 0.4 (axis fraction)
            # Reduced from 0.7 to 0.4 to prevent visual intersection of non-intersecting particles
            axis_length_points = 10 * 72 * 0.4
            
            # Size in points = (Diameter / Data Range) * Axis Length in Points
            # Scatter 's' argument is area in points^2
            s_values = ((diameters / max_range) * axis_length_points) ** 2
            # print(f"Frame {frame_idx}: Using diameter-based sizes.")
            # print(f"s_values {s_values}")
        else:
            s_values = point_size

        # Setup plot
        fig = plt.figure(figsize=(10, 10), dpi=resolution_dpi)
        ax = fig.add_subplot(111, projection='3d')
        
        # Set view angle
        try:
            ax.view_init(elev=elev, azim=azim, roll=roll)
        except TypeError:
            # Fallback for older matplotlib versions that don't support roll
            if roll != 0:
                print(f"Warning: 'roll' parameter ignored. Matplotlib version may be too old.")
            ax.view_init(elev=elev, azim=azim)
            
        # Ensure the aspect ratio is cubic so the fixed limits look correct
        try:
            ax.set_box_aspect((1, 1, 1))
        except AttributeError:
            # Older matplotlib versions might use pbaspect or not support this
            pass
        
        # Scatter plot
        sc = ax.scatter(x, y, z, c=c_values, cmap=cmap, s=s_values)
        
        # Set consistent limits
        ax.set_xlim(limits[0])
        ax.set_ylim(limits[1])
        ax.set_zlim(limits[2])
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f"Timestep: {timestep}")
        
        # Save temp file
        temp_filename = os.path.join(temp_dir, f"temp_frame_{frame_idx:05d}.png")
        plt.savefig(temp_filename)
        plt.close(fig)
        return temp_filename

    def create_animation(self, fps=30, resolution_dpi=100, color_by='id', cmap='viridis', point_size=50, start_frame=0, end_frame=None, view='isometric', axis_limits=None):
        """
        Main function to create the animation using the loaded DataFrame.
        
        :param start_frame: Index of the first frame to include (0-based).
        :param end_frame: Index of the last frame to include (exclusive). If None, includes all.
        :param view: Camera view. Can be a tuple (elev, azim) or string: 'isometric', 'top'/'xy', 'front'/'xz', 'side'/'yz'.
        :param axis_limits: Optional. Can be a float (range size centered on data) or a list of tuples [(xmin, xmax), (ymin, ymax), (zmin, zmax)].
        """
        # Determine view angles
        roll = 0
        if isinstance(view, (tuple, list)):
            if len(view) == 2:
                elev, azim = view
            elif len(view) == 3:
                elev, azim, roll = view
        elif view in ['top', 'xy']:
            elev, azim = 90, -90
        elif view in ['front', 'xz']:
            elev, azim = 0, -90
        elif view in ['front_rev', '-xz']:
            elev, azim = 0, 90
        elif view in ['side', 'yz']:
            elev, azim = 0, 0
        elif view == 'z_left_x_down':
            elev, azim, roll = 0, 90, 90
        else: # isometric or default
            elev, azim = 30, -60
        
        view_angles = (elev, azim, roll)

        # Get all unique timesteps sorted
        all_timesteps = sorted(self.df['timestep'].unique())
        
        # Handle end_frame default
        if end_frame is None:
            end_frame = len(all_timesteps)
            
        # Slice the timesteps
        selected_timesteps = all_timesteps[start_frame:end_frame]
        
        if not selected_timesteps:
            print("No frames selected for animation.")
            return

        print(f"Preparing animation for {len(selected_timesteps)} frames (indices {start_frame} to {end_frame})...")
        print(f"View: {view} (Elev: {elev}, Azim: {azim})")
        
        # Create temp directory
        temp_dir = "temp_frames"
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        os.makedirs(temp_dir)

        # Filter DataFrame for selected timesteps
        subset_df = self.df[self.df['timestep'].isin(selected_timesteps)]
        
        # Group data by timestep
        grouped = subset_df.groupby('timestep')
        
        # Prepare arguments for parallel processing
        tasks = []
        
        # Determine limits
        if axis_limits is not None:
            if isinstance(axis_limits, (int, float)):
                # Calculate center of the data
                x_center = (self.x_limits[0] + self.x_limits[1]) / 2
                y_center = (self.y_limits[0] + self.y_limits[1]) / 2
                z_center = (self.z_limits[0] + self.z_limits[1]) / 2
                
                half_range = axis_limits / 2
                limits = (
                    (x_center - half_range, x_center + half_range),
                    (y_center - half_range, y_center + half_range),
                    (z_center - half_range, z_center + half_range)
                )
                print(f"Using fixed axis range of {axis_limits} centered at ({x_center:.4f}, {y_center:.4f}, {z_center:.4f})")
            elif isinstance(axis_limits, (list, tuple)) and len(axis_limits) == 3:
                limits = axis_limits
                print(f"Using explicit axis limits: {limits}")
            else:
                print("Warning: Invalid axis_limits format. Using auto-calculated limits.")
                limits = (self.x_limits, self.y_limits, self.z_limits)
        else:
            limits = (self.x_limits, self.y_limits, self.z_limits)
        
        for i, (timestep, group) in enumerate(grouped):
            tasks.append((i, timestep, group, limits, color_by, cmap, resolution_dpi, point_size, temp_dir, view_angles))
            
        temp_files = []
        
        print(f"Rendering {len(tasks)} frames in parallel...")
        start_time = time.time()

        # Use ProcessPoolExecutor for CPU-bound plotting tasks
        with ProcessPoolExecutor() as executor:
            # Submit all tasks
            futures = [executor.submit(self._render_wrapper, task) for task in tasks]
            
            # Track progress
            total_frames = len(futures)
            for i, future in enumerate(as_completed(futures)):
                temp_files.append(future.result())
                
                # Progress update
                completed = i + 1
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                remaining = (total_frames - completed) / rate if rate > 0 else 0
                
                print(f"\rProgress: {completed}/{total_frames} ({completed/total_frames*100:.1f}%) - "
                      f"Elapsed: {elapsed:.1f}s - ETA: {remaining:.1f}s", end="", flush=True)
        
        print() # Newline after progress bar
            
        print("Encoding video with ffmpeg...")
        
        # Try NVIDIA GPU encoding first
        ffmpeg_cmd_gpu = [
            'ffmpeg',
            '-y',
            '-framerate', str(fps),
            '-i', os.path.join(temp_dir, 'temp_frame_%05d.png'),
            '-c:v', 'h264_nvenc',  # NVIDIA GPU encoder
            '-pix_fmt', 'yuv420p',
            self.output_file
        ]

        # Fallback to CPU encoding
        ffmpeg_cmd_cpu = [
            'ffmpeg',
            '-y',
            '-framerate', str(fps),
            '-i', os.path.join(temp_dir, 'temp_frame_%05d.png'),
            '-c:v', 'libx264',     # CPU encoder
            '-pix_fmt', 'yuv420p',
            self.output_file
        ]
        
        try:
            print("Attempting GPU encoding (h264_nvenc)...")
            subprocess.run(ffmpeg_cmd_gpu, check=True, stderr=subprocess.PIPE)
            print(f"Animation saved to {self.output_file} (GPU Encoded)")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("GPU encoding failed or not available. Falling back to CPU encoding (libx264)...")
            try:
                subprocess.run(ffmpeg_cmd_cpu, check=True)
                print(f"Animation saved to {self.output_file} (CPU Encoded)")
            except FileNotFoundError:
                print("Error: ffmpeg not found. Please install ffmpeg and add it to PATH.")
            except subprocess.CalledProcessError as e:
                print(f"Error running ffmpeg: {e}")
        finally:
            print("Cleaning up temporary files...")
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    @staticmethod
    def _render_wrapper(args):
        """Helper to unpack arguments for the static method"""
        return Animator._render_single_frame(*args)