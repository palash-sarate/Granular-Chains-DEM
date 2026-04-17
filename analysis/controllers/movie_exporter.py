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
        self.dt_var = tk.StringVar(value="0.000001") # Default dt
        # Scalebar and timestamp UI options
        self.scalebar_width_var = tk.StringVar(value="0.2")   # fraction of view width
        self.scalebar_height_var = tk.StringVar(value="0.02") # fraction of view height / font scale for Text2D fallback
        self.timestamp_pos_var = tk.StringVar(value="top-right")
        self.timestamp_font_var = tk.StringVar(value="0.8")   # Text2D font size scalar
        
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
        # Scalebar size controls
        tk.Label(settings_frame, text="Scalebar W:").grid(row=4, column=0, sticky='w')
        tk.Entry(settings_frame, textvariable=self.scalebar_width_var, width=8).grid(row=4, column=1, padx=2, pady=2, sticky='w')
        tk.Label(settings_frame, text="Scalebar H:").grid(row=5, column=0, sticky='w')
        tk.Entry(settings_frame, textvariable=self.scalebar_height_var, width=8).grid(row=5, column=1, padx=2, pady=2, sticky='w')
        # Timestamp controls
        tk.Label(settings_frame, text="Timestamp Pos:").grid(row=6, column=0, sticky='w')
        tk.OptionMenu(settings_frame, self.timestamp_pos_var, "top-right", "top-left", "bottom-right", "bottom-left").grid(row=6, column=1, padx=2, pady=2, sticky='w')
        tk.Label(settings_frame, text="Timestamp Size:").grid(row=7, column=0, sticky='w')
        tk.Entry(settings_frame, textvariable=self.timestamp_font_var, width=8).grid(row=7, column=1, padx=2, pady=2, sticky='w')

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
            # if self.data_ctrl.geometry_data and 'dt' in self.data_ctrl.geometry_data:
            #     self.dt_var.set(str(self.data_ctrl.geometry_data['dt']))

    def _browse_file(self):
        # Default initial directory to the current simulation folder (if available)
        initial_dir = None
        if getattr(self.data_ctrl, 'sim_source', None):
            initial_dir = getattr(self.data_ctrl.sim_source, 'data_dir', None)
        if not initial_dir and getattr(self.data_ctrl, 'current_sim_folder', None):
            initial_dir = self.data_ctrl.current_sim_folder
        if not initial_dir:
            initial_dir = os.getcwd()

        initial_file = os.path.basename(self.filename_var.get()) if self.filename_var.get() else "simulation_movie.mp4"
        f = filedialog.asksaveasfilename(initialdir=initial_dir, initialfile=initial_file, defaultextension=".mp4", filetypes=[("MP4 Video", "*.mp4"), ("GIF Image", "*.gif"), ("All Files", "*.*")])
        if f:
            # Normalize and store the chosen path
            self.filename_var.set(os.path.normpath(f))

    def _set_plotter_ui_visible(self, visible: bool):
        """Show or hide interactive UI elements that live inside the vedo Plotter (buttons, sliders)."""
        try:
            # Buttons (returned by add_button)
            for b in getattr(self.plotter, 'buttons', []) or []:
                try:
                    if visible:
                        if hasattr(b, 'on') and callable(b.on):
                            b.on()
                        elif hasattr(b, 'enable') and callable(b.enable):
                            b.enable()
                        else:
                            actor = getattr(b, 'actor', None)
                            if actor is not None:
                                actor.SetVisibility(1)
                    else:
                        if hasattr(b, 'off') and callable(b.off):
                            b.off()
                        elif hasattr(b, 'disable') and callable(b.disable):
                            b.disable()
                        else:
                            actor = getattr(b, 'actor', None)
                            if actor is not None:
                                actor.SetVisibility(0)
                except Exception:
                    try:
                        actor = getattr(b, 'actor', None)
                        if actor is not None:
                            actor.SetVisibility(1 if visible else 0)
                    except Exception:
                        pass

            # Sliders (returned by add_slider)
            for s in getattr(self.plotter, 'sliders', []) or []:
                try:
                    if visible and hasattr(s, 'on') and callable(s.on):
                        s.on()
                    elif (not visible) and hasattr(s, 'off') and callable(s.off):
                        s.off()
                except Exception:
                    try:
                        actor = getattr(s, 'actor', None)
                        if actor is not None:
                            actor.SetVisibility(1 if visible else 0)
                    except Exception:
                        pass
        except Exception:
            # Best-effort: silently ignore if plotter doesn't expose these attributes
            pass

    def _init_overlays(self):
        """Initializes constant overlays and creates the text actor for timestamps once."""
        # Clear previous movie-specific actors
        for act in self._movie_actors:
            try: self.plotter.remove(act)
            except: pass
        self._movie_actors = []
        self._timestamp_actor = None
        self._scalebar_actor = None

        if self.show_scalebar_var.get():
            sb = None
            try:
                sb_w = float(self.scalebar_width_var.get())
            except Exception:
                sb_w = 0.2
            try:
                sb_h = float(self.scalebar_height_var.get())
            except Exception:
                sb_h = 0.02

            # Try a few possible vedo Plotter method names for adding a scalebar
            for method in ("add_scalebar", "addScalarBar", "add_scalar_bar", "addScalarbar"):
                fn = getattr(self.plotter, method, None)
                if callable(fn):
                    try:
                        try:
                            sb = fn(pos=(0.8, 0.05), c='black', width=sb_w, height=sb_h)
                        except TypeError:
                            try:
                                sb = fn(pos=(0.8, 0.05), c='black', size=(sb_w, sb_h))
                            except TypeError:
                                try:
                                    sb = fn(pos=(0.8, 0.05), c='black')
                                except TypeError:
                                    sb = fn()
                    except Exception:
                        sb = None
                    break

            if sb:
                self._movie_actors.append(sb)
                self._scalebar_actor = sb
            else:
                # Fallback: create a Text2D based bar using block characters
                chars = max(3, int(sb_w * 40))
                bar_str = '█' * chars
                s_font = max(0.4, sb_h * 20)
                s_text = Text2D(bar_str, pos='bottom-right', s=s_font, c='black', bg='white', alpha=0.9)
                try:
                    self.plotter.add(s_text)
                except Exception:
                    pass
                self._movie_actors.append(s_text)
                self._scalebar_actor = s_text

        if self.show_timestamp_var.get():
            pos = self.timestamp_pos_var.get() if hasattr(self, 'timestamp_pos_var') else 'top-right'
            try: s_font = float(self.timestamp_font_var.get())
            except Exception: s_font = 0.8
            t2d = Text2D("Preparing...", pos=pos, s=s_font, c='black', bg='white', alpha=0.7)
            self.plotter.add(t2d)
            self._movie_actors.append(t2d)
            self._timestamp_actor = t2d

    def _update_timestamp(self, ts):
        """Fast update of timestamp text without recreating VTK actors."""
        if getattr(self, '_timestamp_actor', None):
            try:
                dt = float(self.dt_var.get())
                time = ts * dt
                txt = f"Time: {time:.4f} s\nStep: {ts}"
            except:
                txt = f"Step: {ts}"
            self._timestamp_actor.text(txt)

    def toggle_preview(self):
        if self.preview_active:
            # Stop preview and restore UI
            self.preview_active = False
            self.preview_btn.config(text="Preview", bg='#e1f5fe')
            self.status_label.config(text="Preview stopped")
            try: self._set_plotter_ui_visible(True)
            except: pass
        else:
            # Hide vedo on-screen buttons while previewing
            try: self._set_plotter_ui_visible(False)
            except: pass
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
            try:
                fps_preview = int(self.fps_var.get())
                if fps_preview <= 0:
                    fps_preview = 30
            except Exception:
                fps_preview = 30
            delay_ms = max(1, int(1000 / fps_preview))
            self._init_overlays()

            def step_loop(idx):
                if not self.preview_active or idx >= len(timesteps):
                    if self.preview_active: # Finished normally
                        self.toggle_preview()
                    return
                
                ts = timesteps[idx]
                self.renderer.show_timestep(ts)
                self._update_timestamp(ts)
                self.plotter.render()
                
                # Schedule next frame according to FPS setting
                self.parent.after(delay_ms, lambda: step_loop(idx + 1))
            
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

        # If the user provided a relative filename (or just a basename),
        # save it inside the simulation output directory when available.
        if not os.path.isabs(filename):
            sim_dir = None
            if getattr(self.data_ctrl, 'sim_source', None):
                sim_dir = getattr(self.data_ctrl.sim_source, 'data_dir', None)
            if not sim_dir and getattr(self.data_ctrl, 'current_sim_folder', None):
                sim_dir = self.data_ctrl.current_sim_folder
            if sim_dir:
                filename = os.path.join(sim_dir, filename)
            else:
                filename = os.path.join(os.getcwd(), filename)
            filename = os.path.normpath(filename)
            self.filename_var.set(filename)

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
        # Hide vedo UI (buttons/sliders) during export so they are not captured
        try: self._set_plotter_ui_visible(False)
        except: pass
        
        # vedo Video setup
        # Note: on Windows, ffmpeg must be installed.
        video = Video(filename, fps=fps)
        
        try:
            self._init_overlays()
            total = len(timesteps)
            for i, ts in enumerate(timesteps):
                self.status_label.config(text=f"Exporting frame {i+1}/{total}...")
                
                if i % 5 == 0:
                    # Update UI less frequently to avoid main thread event pump bottleneck
                    self.parent.update() 
                
                # Render frame
                # We use the renderer's normal flow, but we can't use 'after' here
                # because we want to stick to the export speed.
                self.renderer.show_timestep(ts)
                self._update_timestamp(ts)
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
            # Restore vedo UI now that export is done
            try: self._set_plotter_ui_visible(True)
            except: pass
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
