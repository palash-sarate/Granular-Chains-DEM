import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional, Dict, Callable
import pandas as pd
from vedo import Plotter
from analysis.controllers import (
    SimDataController, SimulationRenderer, HighlightController,
    PlaybackController, SimulationLoader, VtkOverlayController,
    AnalysisToolManager
)
from analysis.utils.progress_bar import ProgressBar

# Optional drag-and-drop support via tkinterdnd2. If not available,
# the UI will show an instruction and Open buttons remain functional.
try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    BaseTk = TkinterDnD.Tk
    DND_CONST = DND_FILES
    DND_AVAILABLE = True
except Exception:
    BaseTk = tk.Tk
    DND_AVAILABLE = False

class ViewerApp(BaseTk):
    def __init__(self):
        super().__init__()
        self.title('Chains Simulation Viewer')
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        # 1. Base Controller & Data
        self.data_ctrl = SimDataController()
        self.data_ctrl.on_data_loaded_cb = self._on_data_loaded
        self.data_ctrl.on_batch_ready_cb = self._on_batch_ready
        
        # 2. Progress State
        self.progress_bar = None
        self.batches_to_load = 0
        self.batches_loaded = 0
        self.enable_preloading_var = tk.BooleanVar(value=True)
        
        # 2. Vedo Scene & Render Control
        self.plotter = Plotter(bg='white', interactive=True)
        self.vtk_ctrl = VtkOverlayController(self.plotter)
        self.hl_ctrl = HighlightController(self.plotter, self._get_current_rendering_data, self._get_chain_size)
        self.renderer = SimulationRenderer(self.plotter, self.data_ctrl, self.vtk_ctrl, self.hl_ctrl, self._get_ui_settings)
        
        # 3. Logic Controllers
        self.playback_ctrl = PlaybackController(self, lambda: self.data_ctrl.timesteps, self.renderer.show_timestep)
        self.analysis_ctrl = AnalysisToolManager(self, self.data_ctrl, self._get_ui_vars)
        
        # UI Callbacks for Loader
        loader_cbs = {
            'on_load_success': self._on_data_loaded,
            'get_current_ts': lambda: self.current_timestep,
            'sync_ui_to_frame': self._sync_ui_to_frame,
            'clear_vtk_list': lambda: self.vtk_listbox.delete(0, tk.END),
            'add_vtk': self._add_vtk_mesh_ui,
            'enable_preloading': self.enable_preloading_var
        }
        self.loader_ctrl = SimulationLoader(self, self.data_ctrl, self.playback_ctrl, self.vtk_ctrl, self.renderer, loader_cbs)

        # ── UI LAYOUT ───────────────────────────────────────────
        ctrl = tk.Frame(self)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        col1 = tk.Frame(ctrl); col1.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
        try:
            import tkinter.ttk as ttk
            ttk.Separator(ctrl, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=4)
        except Exception:
            tk.Frame(ctrl, width=1, bg='gray').pack(side=tk.LEFT, fill=tk.Y, padx=4)
        col2 = tk.Frame(ctrl); col2.pack(side=tk.LEFT, fill=tk.Y)

        # ── Column 1 : File & Playback ──────────────────────────
        open_row1 = tk.Frame(col1)
        open_row1.pack(fill=tk.X, pady=(4, 0))
        tk.Button(open_row1, text='Open Folder...', command=self.loader_ctrl.open_dump_folder).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(open_row1, text='Open Dumps...', command=self.loader_ctrl.open_dump_files).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

        open_row2 = tk.Frame(col1)
        open_row2.pack(fill=tk.X, pady=(4, 0))
        tk.Button(open_row2, text='Open Data...', command=self.loader_ctrl.open_data_file).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(open_row2, text='Open VTKs...', command=self._open_vtk_files_ui).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        self.refresh_button = tk.Button(col1, text='Refresh', command=lambda: self.loader_ctrl.load_simulation_folder(self.data_ctrl.current_sim_folder, False), state=tk.DISABLED)
        self.refresh_button.pack(fill=tk.X, pady=(8, 0))
        tk.Checkbutton(col1, text='Enable Preloading', variable=self.enable_preloading_var).pack(anchor='w', pady=(2, 0))

        drop_text = 'Drop files here' if DND_AVAILABLE else 'DND disabled'
        self.drop_label = tk.Label(col1, text=drop_text, relief='ridge', width=22, height=4)
        self.drop_label.pack(fill=tk.X, pady=(8, 4))
        if DND_AVAILABLE:
            try:
                self.drop_label.drop_target_register(DND_CONST)
                self.drop_label.dnd_bind('<<Drop>>', self._on_drop_ui)
            except Exception: pass

        tk.Label(col1, text='Timesteps:').pack(anchor='w', pady=(8, 0))
        self.ts_listbox = tk.Listbox(col1, width=22, height=6)
        self.ts_listbox.pack(fill=tk.Y)
        self.ts_listbox.bind('<<ListboxSelect>>', self._on_ts_select_ui)

        tk.Label(col1, text='Frame Navigation:').pack(anchor='w', pady=(8, 0))
        playback_frame = tk.Frame(col1); playback_frame.pack(fill=tk.X, pady=(2, 0))
        btn_row = tk.Frame(playback_frame); btn_row.pack(fill=tk.X)
        tk.Button(btn_row, text='|◀', width=3, command=self.playback_ctrl.jump_start).pack(side=tk.LEFT)
        tk.Button(btn_row, text='◀', width=3, command=self.playback_ctrl.step_prev).pack(side=tk.LEFT, padx=2)
        self.play_btn = tk.Button(btn_row, text='▶ Play', width=8, command=self.playback_ctrl.toggle_play); self.play_btn.pack(side=tk.LEFT, padx=2)
        tk.Button(btn_row, text='▶', width=3, command=self.playback_ctrl.step_next).pack(side=tk.LEFT, padx=2)
        tk.Button(btn_row, text='▶|', width=3, command=self.playback_ctrl.jump_end).pack(side=tk.LEFT)

        self.frame_slider = tk.Scale(playback_frame, from_=0, to=0, orient=tk.HORIZONTAL, command=self._on_slider_ui)
        self.frame_slider.pack(fill=tk.X, pady=(4, 0))
        self.frame_label = tk.Label(playback_frame, text='Frame: 0 / 0', fg='gray'); self.frame_label.pack(fill=tk.X)
        
        playback_frame = tk.Frame(col1)
        playback_frame.pack(fill=tk.X, pady=(2, 0))
        tk.Label(playback_frame, text='FPS:').pack(side=tk.LEFT)
        self.fps_entry = tk.Entry(playback_frame, width=4)
        self.fps_entry.insert(0, '10')
        self.fps_entry.pack(side=tk.LEFT, padx=2)

        self.loop_var = tk.BooleanVar(value=True)
        tk.Checkbutton(playback_frame, text='Loop', variable=self.loop_var).pack(side=tk.LEFT, padx=(4, 0))

        self.playback_ctrl.link_widgets(self.frame_slider, self.fps_entry, self.play_btn, self.frame_label, self.loop_var)

        view_btns = tk.Frame(col1)
        view_btns.pack(fill=tk.X, pady=(12, 0))
        tk.Button(view_btns, text='Fit View', command=self.renderer.fit_view).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(view_btns, text='Reset View', command=self.renderer.reset_view).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        # (Geometry checkbox moved to col2 next to Delete VTKs)

        # ── Column 2 : VTK, Highlight & Analysis ───────────────
        tk.Label(col2, text='Loaded VTKs:', anchor='w').pack(fill=tk.X, pady=(0, 2))
        vtk_frame = tk.Frame(col2); vtk_frame.pack(fill=tk.X)
        self.vtk_listbox = tk.Listbox(vtk_frame, width=22, height=5, selectmode=tk.MULTIPLE, exportselection=False)
        self.vtk_listbox.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.vtk_listbox.bind('<<ListboxSelect>>', self._on_vtk_select_ui)
        vtk_scroll = tk.Scrollbar(vtk_frame, orient=tk.VERTICAL); vtk_scroll.config(command=self.vtk_listbox.yview); vtk_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.vtk_listbox.config(yscrollcommand=vtk_scroll.set)
        vtk_btns = tk.Frame(col2)
        vtk_btns.pack(fill=tk.X, pady=(2, 0))
        tk.Button(vtk_btns, text='Delete Selected VTKs', command=self._delete_selected_vtks_ui).pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.show_geometry_var = tk.BooleanVar(value=True)
        tk.Checkbutton(vtk_btns, text='Show/Hide', variable=self.show_geometry_var, command=self._on_draw_toggle_ui).pack(side=tk.LEFT, padx=(4, 0))

        tk.Label(col2, text='──── Highlight ────', fg='gray').pack(fill=tk.X, pady=(12, 2))
        hl_mode_frame = tk.Frame(col2); hl_mode_frame.pack(fill=tk.X)
        self.highlight_mode_var = tk.StringVar(value='atom')
        self.highlight_id_var = tk.StringVar(value='1')
        self.chain_size_var = tk.StringVar(value='4')
        for m, l in [('atom', 'Atom'), ('bond', 'Bond'), ('angle', 'Angle'), ('chain', 'Chain')]:
            tk.Radiobutton(hl_mode_frame, text=l, variable=self.highlight_mode_var, value=m).pack(side=tk.LEFT)

        hl_input_row = tk.Frame(col2)
        hl_input_row.pack(fill=tk.X, pady=(4, 0))
        tk.Label(hl_input_row, text='ID:').pack(side=tk.LEFT)
        
        def adjust_hl_id(delta):
            try:
                val = int(self.highlight_id_var.get())
                self.highlight_id_var.set(str(max(1, val + delta)))
                self._apply_highlight_ui()
            except ValueError:
                pass

        tk.Button(hl_input_row, text='-', width=2, command=lambda: adjust_hl_id(-1)).pack(side=tk.LEFT)
        self.hl_id_entry = tk.Entry(hl_input_row, textvariable=self.highlight_id_var, width=5)
        self.hl_id_entry.pack(side=tk.LEFT, padx=2)
        tk.Button(hl_input_row, text='+', width=2, command=lambda: adjust_hl_id(1)).pack(side=tk.LEFT)

        tk.Label(hl_input_row, text='N:').pack(side=tk.LEFT, padx=(10, 0))
        self.chain_size_entry = tk.Entry(hl_input_row, textvariable=self.chain_size_var, width=6)
        self.chain_size_entry.pack(side=tk.LEFT, padx=5)

        hl_btn_row = tk.Frame(col2)
        hl_btn_row.pack(fill=tk.X, pady=(8, 0))
        tk.Button(hl_btn_row, text='Apply', command=self._apply_highlight_ui).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(hl_btn_row, text='Clear', command=self._clear_highlight_ui).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        
        self.highlight_status = tk.Label(col2, text='—', fg='gray', wraplength=160, justify='left'); self.highlight_status.pack(anchor='w', pady=(4, 0))

        tk.Label(col2, text='──── Analysis ────', fg='gray').pack(fill=tk.X, pady=(12, 2))
        analysis_row1 = tk.Frame(col2)
        analysis_row1.pack(fill=tk.X, pady=(2, 0))
        tk.Button(analysis_row1, text='Bonds…', command=self.analysis_ctrl.open_bond_win).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(analysis_row1, text='Angles…', command=self.analysis_ctrl.open_angle_win).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

        analysis_row2 = tk.Frame(col2)
        analysis_row2.pack(fill=tk.X, pady=(4, 0))
        tk.Button(analysis_row2, text='Atoms…', command=self.analysis_ctrl.open_atom_win).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(analysis_row2, text='Lepton…', command=self.analysis_ctrl.open_lepton_win).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

        self.current_timestep = None
        self._autosize_and_center()
        self.bind('<space>', lambda e: self.playback_ctrl.toggle_play())
        self.bind('<Left>', lambda e: self.playback_ctrl.step_prev())
        self.bind('<Right>', lambda e: self.playback_ctrl.step_next())

    # ── Internal Helpers ─────────────────────────────────────
    def _get_ui_settings(self):
        return {'show_geometry': self.show_geometry_var.get()}

    def _get_ui_vars(self):
        return {'chain_size': self.chain_size_var, 'enable_preloading': self.enable_preloading_var}

    def _get_current_rendering_data(self):
        return self.data_ctrl.df, self.current_timestep

    def _get_chain_size(self):
        try: return int(self.chain_size_var.get())
        except: return 4

    def _on_data_loaded(self, num_queued=0):
        # Handle progress bar if we have background batches to wait for
        if num_queued > 0:
            self.batches_to_load = num_queued
            self.batches_loaded = 0
            self.progress_bar = ProgressBar(self, title="Preloading Cache", message=f"Loading {num_queued} cached batches...", modal=True)
            self.progress_bar.set_progress(0, num_queued)
        else:
            self.batches_to_load = 0
            
        old_ts = self.current_timestep
        self.analysis_ctrl.refresh_windows()
        self.ts_listbox.delete(0, tk.END)
        for t in self.data_ctrl.timesteps: 
            self.ts_listbox.insert(tk.END, str(t))
            
        timesteps = self.data_ctrl.timesteps
        total_frames = len(timesteps)
        self.frame_slider.config(from_=0, to=max(0, total_frames - 1))
        
        target_idx = 0
        if old_ts is not None and old_ts in timesteps:
            target_idx = list(timesteps).index(old_ts)
        
        self.frame_slider.set(target_idx)
        self.playback_ctrl.update_status(target_idx)
        self.ts_listbox.selection_clear(0, tk.END)
        self.ts_listbox.selection_set(target_idx)
        self.ts_listbox.activate(target_idx)
        
        self.renderer.update_persistent_bounds()
        if timesteps:
            ts = timesteps[target_idx]
            self.current_timestep = ts
            self.renderer.show_timestep(ts)
        self.refresh_button.config(state=tk.NORMAL)

    def _on_batch_ready(self, batch_idx):
        """Callback from data controller when a background batch is loaded."""
        # Ensure UI updates happen on the main thread
        self.after(0, lambda: self._on_batch_ready_safe(batch_idx))

    def _on_batch_ready_safe(self, batch_idx):
        pb = getattr(self, 'progress_bar', None)
        if pb is not None:
            self.batches_loaded += 1
            pb.set_progress(self.batches_loaded, self.batches_to_load)
            pb.set_text(f"Loaded batch {batch_idx+1} of {self.batches_to_load + 1}...")
            
            if self.batches_loaded >= self.batches_to_load:
                self.progress_bar = None
                pb.finish()
        
        print(f"Batch {batch_idx} ready. Refreshing UI...")
        self._refresh_current_view()

    def _refresh_current_view(self):
        # Determine current timestep
        idx = int(self.frame_slider.get())
        if 0 <= idx < len(self.data_ctrl.timesteps):
            ts = self.data_ctrl.timesteps[idx]
            self.current_timestep = ts
            self.renderer.show_timestep(ts)
            self.analysis_ctrl.refresh_windows()

    def _sync_ui_to_frame(self, idx, ts):
        self.frame_slider.set(idx)
        self.ts_listbox.selection_clear(0, tk.END)
        self.ts_listbox.selection_set(idx)
        self.ts_listbox.activate(idx)
        self.current_timestep = ts
        self.renderer.show_timestep(ts)

    def _add_vtk_mesh_ui(self, path):
        ok, res = self.vtk_ctrl.add_mesh(path)
        if ok:
            self.vtk_listbox.insert(tk.END, res)
            self.vtk_listbox.selection_set(self.vtk_listbox.size() - 1)
            self.renderer.update_persistent_bounds()
            if self.current_timestep: self.renderer.show_timestep(self.current_timestep)
        else: print(f"VTK load failed: {res}")

    def _open_vtk_files_ui(self):
        files = filedialog.askopenfilenames(title='Select VTK files', filetypes=[('VTK files', '*.vtk'), ('All', '*.*')])
        for f in files: self._add_vtk_mesh_ui(f)
        self.plotter.render()

    def _on_drop_ui(self, event):
        files = list(self.tk.splitlist(event.data))
        self.loader_ctrl.handle_dropped_files(files)

    def _on_ts_select_ui(self, event):
        sel = self.ts_listbox.curselection()
        if sel:
            idx = sel[0]; ts = self.data_ctrl.timesteps[idx]
            self.current_timestep = ts
            self.frame_slider.set(idx); self.renderer.show_timestep(ts)

    def _on_slider_ui(self, val):
        idx = int(float(val))
        if 0 <= idx < len(self.data_ctrl.timesteps):
            ts = self.data_ctrl.timesteps[idx]
            self.current_timestep = ts
            self.ts_listbox.selection_clear(0, tk.END); self.ts_listbox.selection_set(idx); self.ts_listbox.activate(idx)
            self.playback_ctrl.update_status(idx); self.renderer.show_timestep(ts)

    def _on_vtk_select_ui(self, event):
        sel = self.vtk_listbox.curselection()
        for i in range(self.vtk_listbox.size()):
            self.vtk_ctrl.set_visibility(self.vtk_listbox.get(i), i in sel)
        if self.current_timestep: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def _delete_selected_vtks_ui(self):
        sel = self.vtk_listbox.curselection()
        if not sel: return
        names = [self.vtk_listbox.get(i) for i in sel]
        self.vtk_ctrl.remove_meshes(names)
        for i in reversed(sel): self.vtk_listbox.delete(i)
        self.renderer.update_persistent_bounds()
        if self.current_timestep: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def _apply_highlight_ui(self):
        try:
            ok, res = self.hl_ctrl.apply(self.highlight_mode_var.get(), int(self.highlight_id_var.get().strip()))
            self.highlight_status.config(text=res, fg='black' if ok else 'red')
            if ok and self.current_timestep:
                self.renderer.show_timestep(self.current_timestep)
            else:
                self.plotter.render()
        except: self.highlight_status.config(text='Invalid ID', fg='red')

    def _clear_highlight_ui(self):
        self.hl_ctrl.clear(); self.highlight_status.config(text='—', fg='gray'); self.plotter.render()

    def _on_draw_toggle_ui(self):
        self.vtk_ctrl.set_master_visibility(self.show_geometry_var.get())
        if self.current_timestep: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()


    def _autosize_and_center(self) -> None:
        try:
            self.update_idletasks()
            w = max(self.winfo_reqwidth(), self.winfo_width())
            h = max(self.winfo_reqheight(), self.winfo_height())
            sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
            x, y = max(0, int((sw - w) / 2)), max(0, int((sh - h) / 2))
            self.geometry(f"{w}x{h}+{x}+{y}")
            self.minsize(w, h)
        except: pass

    def _cleanup(self):
        self.playback_ctrl.pause()
        try:
            import matplotlib.pyplot as plt
            plt.close('all')
        except: pass

    def _on_close(self):
        self._cleanup()
        try: self.destroy()
        except: pass
        sys.exit(0)


def run():
    app = ViewerApp()
    try:
        app.mainloop()
    finally:
        # In case mainloop returns, ensure cleanup and exit
        try:
            app._cleanup()
        except Exception:
            pass
        sys.exit(0)


if __name__ == '__main__':
    run()
