import os
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
from vedo import Plotter, Video, Text2D
from typing import Optional, List, Callable

class MovieExporterController:
    """Manages simulation movie export with custom ranges, scalebars, and timestamps."""
    
    def __init__(self, parent: tk.Frame, data_ctrl, renderer, plotter: Plotter):
        self.parent = parent
        self.data_ctrl = data_ctrl
        self.renderer = renderer
        self.plotter = plotter
        
        self.preview_active = False
        self._exporting = False
        self._movie_actors = []
        
        # UI Variables
        self.start_ts_var = tk.StringVar(value="0")
        self.end_ts_var = tk.StringVar(value="0")
        self.fps_var = tk.StringVar(value="24")
        self.filename_var = tk.StringVar(value="simulation_movie.mp4")
        self.show_scalebar_var = tk.BooleanVar(value=False)
        self.show_timestamp_var = tk.BooleanVar(value=True)
        self.dt_var = tk.StringVar(value="0.001") # Default dt
        
        self._setup_ui()

    def _setup_ui(self):
        tk.Label(self.parent, text='──── Movie Exporter ────', fg='gray').pack(fill=tk.X, pady=(4, 2))
        
        # Timestep Range
        range_frame = tk.LabelFrame(self.parent, text="Timestep Range", padx=5, pady=5)
        range_frame.pack(fill=tk.X, pady=2)
        
        tk.Label(range_frame, text="Start:").grid(row=0, column=0, sticky='w')
        self.start_entry = tk.Entry(range_frame, textvariable=self.start_ts_var, width=10)
        self.start_entry.grid(row=0, column=1, padx=2, pady=2)
        
        tk.Label(range_frame, text="End:").grid(row=1, column=0, sticky='w')
        self.end_entry = tk.Entry(range_frame, textvariable=self.end_ts_var, width=10)
        self.end_entry.grid(row=1, column=1, padx=2, pady=2)
        
        tk.Button(range_frame, text="Set Current", command=self._set_range_from_loader, font=('Arial', 8)).grid(row=2, column=0, columnspan=2, sticky='ew', pady=2)

        # Settings
        settings_frame = tk.LabelFrame(self.parent, text="Settings", padx=5, pady=5)
        settings_frame.pack(fill=tk.X, pady=2)
        
        tk.Label(settings_frame, text="FPS:").grid(row=0, column=0, sticky='w')
        tk.Entry(settings_frame, textvariable=self.fps_var, width=5).grid(row=0, column=1, padx=2, pady=2, sticky='w')
        
        tk.Label(settings_frame, text="dt:").grid(row=1, column=0, sticky='w')
        tk.Entry(settings_frame, textvariable=self.dt_var, width=8).grid(row=1, column=1, padx=2, pady=2, sticky='w')
        
        tk.Checkbutton(settings_frame, text="Scalebar", variable=self.show_scalebar_var).grid(row=2, column=0, columnspan=2, sticky='w')
        tk.Checkbutton(settings_frame, text="Timestamp", variable=self.show_timestamp_var).grid(row=3, column=0, columnspan=2, sticky='w')

        # File
        file_frame = tk.Frame(self.parent)
        file_frame.pack(fill=tk.X, pady=4)
        tk.Label(file_frame, text="File:").pack(side=tk.LEFT)
        tk.Entry(file_frame, textvariable=self.filename_var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)
        tk.Button(file_frame, text="...", command=self._browse_file, width=2).pack(side=tk.LEFT)

        # Actions
        btn_frame = tk.Frame(self.parent)
        btn_frame.pack(fill=tk.X, pady=4)
        
        self.preview_btn = tk.Button(btn_frame, text="Preview", command=self.toggle_preview, bg='#e1f5fe')
        self.preview_btn.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 2))
        
        self.export_btn = tk.Button(btn_frame, text="Export Movie", command=self.start_export, bg='#e8f5e9')
        self.export_btn.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(2, 0))

        self.status_label = tk.Label(self.parent, text="Ready", fg='gray')
        self.status_label.pack(fill=tk.X)

    def _set_range_from_loader(self):
        timesteps = self.data_ctrl.timesteps
        if timesteps:
            self.start_ts_var.set(str(timesteps[0]))
            self.end_ts_var.set(str(timesteps[-1]))
            # Also try to update dt from geometry data if available
            if self.data_ctrl.geometry_data and 'dt' in self.data_ctrl.geometry_data:
                self.dt_var.set(str(self.data_ctrl.geometry_data['dt']))

    def _browse_file(self):
        f = filedialog.asksaveasfilename(defaultextension=".mp4", filetypes=[("MP4 Video", "*.mp4"), ("GIF Image", "*.gif"), ("All Files", "*.*")])
        if f:
            self.filename_var.set(f)

    def _update_overlays(self, ts):
        """Adds or updates temporary overlays for recording/preview."""
        # Clear previous movie-specific actors
        for act in self._movie_actors:
            try: self.plotter.remove(act)
            except: pass
        self._movie_actors = []

        if self.show_scalebar_var.get():
            # In vedo, add_scalebar is a method of Plotter that returns the actor
            sb = self.plotter.add_scalebar(pos=(0.8, 0.05), c='black')
            if sb:
                self._movie_actors.append(sb)

        if self.show_timestamp_var.get():
            try:
                dt = float(self.dt_var.get())
                time = ts * dt
                txt = f"Time: {time:.4f} s\nStep: {ts}"
            except:
                txt = f"Step: {ts}"
            
            t2d = Text2D(txt, pos='top-right', s=0.8, c='black', bg='white', alpha=0.7)
            self.plotter.add(t2d)
            self._movie_actors.append(t2d)

    def toggle_preview(self):
        if self.preview_active:
            self.preview_active = False
            self.preview_btn.config(text="Preview", bg='#e1f5fe')
            self.status_label.config(text="Preview stopped")
        else:
            self.preview_active = True
            self.preview_btn.config(text="Stop Preview", bg='#ffecb3')
            self._run_preview()

    def _run_preview(self):
        if not self.preview_active: return
        
        try:
            start = int(self.start_ts_var.get())
            end = int(self.end_ts_var.get())
            timesteps = [t for t in self.data_ctrl.timesteps if start <= t <= end]
            
            if not timesteps:
                messagebox.showwarning("Warning", "No timesteps in selected range.")
                self.toggle_preview()
                return
            
            self.status_label.config(text=f"Previewing {len(timesteps)} frames...")
            
            def step_loop(idx):
                if not self.preview_active or idx >= len(timesteps):
                    if self.preview_active: # Finished normally
                        self.toggle_preview()
                    return
                
                ts = timesteps[idx]
                self.renderer.show_timestep(ts)
                self._update_overlays(ts)
                self.plotter.render()
                
                # Schedule next frame (approx 30fps for preview if possible)
                self.parent.after(10, lambda: step_loop(idx + 1))
            
            step_loop(0)
            
        except Exception as e:
            messagebox.showerror("Error", f"Preview failed: {e}")
            self.toggle_preview()

    def start_export(self):
        if self._exporting: return
        
        if self.preview_active:
            self.toggle_preview() # Stop preview first
            
        filename = self.filename_var.get()
        if not filename:
            messagebox.showerror("Error", "No filename specified.")
            return

        try:
            start = int(self.start_ts_var.get())
            end = int(self.end_ts_var.get())
            fps = int(self.fps_var.get())
            timesteps = [t for t in self.data_ctrl.timesteps if start <= t <= end]
            
            if not timesteps:
                messagebox.showerror("Error", "No timesteps in selected range.")
                return

            if tk.messagebox.askokcancel("Confirm Export", f"Export {len(timesteps)} frames to {os.path.basename(filename)}?\nUI will be locked during export."):
                self._run_export_sync(timesteps, filename, fps)
                
        except Exception as e:
            messagebox.showerror("Error", f"Invalid input: {e}")

    def _run_export_sync(self, timesteps, filename, fps):
        self._exporting = True
        self.status_label.config(text="Exporting... (UI locked)", fg='red')
        self.export_btn.config(state=tk.DISABLED)
        self.preview_btn.config(state=tk.DISABLED)
        
        # vedo Video setup
        # Note: on Windows, ffmpeg must be installed.
        video = Video(filename, fps=fps)
        
        try:
            total = len(timesteps)
            for i, ts in enumerate(timesteps):
                self.status_label.config(text=f"Exporting frame {i+1}/{total}...")
                self.parent.update() # Keep UI responsive enough to show status
                
                # Render frame
                # We use the renderer's normal flow, but we can't use 'after' here
                # because we want to stick to the export speed.
                self.renderer.show_timestep(ts)
                self._update_overlays(ts)
                self.plotter.render()
                
                # Capture
                video.add_frame()
                
            video.close()
            messagebox.showinfo("Success", f"Movie exported successfully to:\n{filename}")
            
        except Exception as e:
            messagebox.showerror("Export Failed", f"An error occurred during export:\n{e}")
            try: video.close()
            except: pass
            
        finally:
            self._exporting = False
            self.status_label.config(text="Ready", fg='gray')
            self.export_btn.config(state=tk.NORMAL)
            self.preview_btn.config(state=tk.NORMAL)
            # Cleanup overlays
            for act in self._movie_actors:
                try: self.plotter.remove(act)
                except: pass
            self._movie_actors = []
            self.plotter.render()
