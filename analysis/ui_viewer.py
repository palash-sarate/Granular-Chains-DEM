import os
import sys
import argparse
import math

# Add project root to sys.path if running as script to support absolute imports
if __name__ == '__main__':
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _root_dir = os.path.dirname(_script_dir)
    if _root_dir not in sys.path:
        sys.path.insert(0, _root_dir)

import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional, Dict, Callable
import pandas as pd
from vedo import Plotter
from analysis.controllers import (
    SimDataController, SimulationRenderer, HighlightController,
    PlaybackController, SimulationLoader, VtkOverlayController,
    AnalysisToolManager, RestartEditorController, MovieExporterController
)
from analysis.utils.progress_bar import ProgressBar
import tempfile
import shutil
from simulation.orchestrator import SimulationOrchestrator
from analysis.controllers.simulation_launcher import SimulationLauncherController

class ColumnMapDialog(tk.Toplevel):
    def __init__(self, parent, filepath):
        super().__init__(parent)
        self.title("Manual Column Mapping")
        self.filepath = filepath
        self.result = None
        self.transient(parent)
        self.grab_set()

        # Read preview
        preview_text = ""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            # Find atoms section for better preview
            atom_start = 0
            for i, line in enumerate(lines):
                if "Atoms" in line:
                    atom_start = i + 1
                    break
            preview_lines = lines[atom_start:atom_start+10]
            preview_text = "".join(preview_lines)
        except Exception as e:
            preview_text = f"Error reading preview: {e}"

        tk.Label(self, text=f"File: {os.path.basename(filepath)}", font=('Arial', 10, 'bold')).pack(pady=5)
        tk.Label(self, text="Preview (first 10 atom lines):").pack(anchor='w', padx=10)
        
        preview_box = tk.Text(self, height=10, width=80, font=('Courier', 9))
        preview_box.insert('1.0', preview_text)
        preview_box.config(state='disabled')
        preview_box.pack(padx=10, pady=5)

        tk.Label(self, text="Enter 0-based column indices (LAMMPS default: id=0, type=1, x=2, y=3, z=4):").pack(pady=5)

        form = tk.Frame(self)
        form.pack(pady=10)

        self.entries = {}
        fields = [('id', 0), ('type', 1), ('x', 2), ('y', 3), ('z', 4), ('mol', 5), ('diameter', 6)]
        
        for i, (field, default) in enumerate(fields):
            tk.Label(form, text=f"{field}:").grid(row=i//2, column=(i%2)*2, padx=5, pady=2, sticky='e')
            e = tk.Entry(form, width=5)
            e.insert(0, str(default))
            e.grid(row=i//2, column=(i%2)*2 + 1, padx=5, pady=2, sticky='w')
            self.entries[field] = e

        btn_frame = tk.Frame(self)
        btn_frame.pack(pady=10)
        tk.Button(btn_frame, text="Apply Mapping", command=self._on_apply, width=15, bg='#4CAF50', fg='white').pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text="Cancel", command=self.destroy, width=10).pack(side=tk.LEFT, padx=5)

        self.geometry("650x500")
        self._center_window()

    def _center_window(self):
        self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (width // 2)
        y = (self.winfo_screenheight() // 2) - (height // 2)
        self.geometry(f'{width}x{height}+{x}+{y}')

    def _on_apply(self):
        try:
            self.result = {k: int(e.get()) for k, e in self.entries.items() if e.get().strip()}
            self.destroy()
        except ValueError:
            messagebox.showerror("Error", "Column indices must be integers.")

class HopperCalculatorController:
    def __init__(self, parent):
        self.parent = parent
        frame = tk.LabelFrame(parent, text="Hopper Fill Calculator", padx=5, pady=5)
        frame.pack(fill=tk.X, pady=5)

        # Inputs
        self.vars = {
            'orifice': tk.DoubleVar(value=10.0),
            'width': tk.DoubleVar(value=31.0),
            'angle': tk.DoubleVar(value=60.0),
            'height': tk.DoubleVar(value=40.0),
            'fill': tk.DoubleVar(value=80.0),
            'bead_d': tk.DoubleVar(value=2.0),
            'n_hoppers': tk.IntVar(value=5)
        }

        fields = [
            ('Orifice (cm)', 'orifice'),
            ('Width (cm)', 'width'),
            ('Angle (°)', 'angle'),
            ('Total H (cm)', 'height'),
            ('Fill %', 'fill'),
            ('Bead D (mm)', 'bead_d'),
            ('N Hoppers', 'n_hoppers')
        ]

        for i, (label, var_name) in enumerate(fields):
            row = tk.Frame(frame)
            row.pack(fill=tk.X, pady=1)
            tk.Label(row, text=label, width=12, anchor='w', font=('Arial', 8)).pack(side=tk.LEFT)
            tk.Entry(row, textvariable=self.vars[var_name], width=8, font=('Arial', 8)).pack(side=tk.RIGHT)

        tk.Button(frame, text="Calculate Counts", command=self.calculate, bg='#2196F3', fg='white', font=('Arial', 9, 'bold')).pack(fill=tk.X, pady=10)

        self.result_text = tk.Text(frame, height=10, width=28, font=('Courier', 9), bg='#f8f8f8')
        self.result_text.pack(fill=tk.X)

    def calculate(self):
        try:
            o_w = self.vars['orifice'].get()
            h_w = self.vars['width'].get()
            ang = self.vars['angle'].get()
            h_tot = self.vars['height'].get()
            f_pct = self.vars['fill'].get()
            b_d = self.vars['bead_d'].get()
            n_h = self.vars['n_hoppers'].get()

            # Area fractions from graph
            phi_values = {4: 0.755, 12: 0.58, 24: 0.56, 48: 0.50}
            b_r_cm = (b_d / 10.0) / 2.0
            b_area = math.pi * (b_r_cm**2)
            ang_rad = math.radians(ang)
            h_conv = math.tan(ang_rad) * (h_w - o_w) / 2.0
            t_h = (f_pct / 100.0) * h_tot

            if t_h <= h_conv:
                w_at_h = o_w + (2.0 * t_h / math.tan(ang_rad))
                area = (o_w + w_at_h) / 2.0 * t_h
            else:
                area_conv = (o_w + h_w) / 2.0 * h_conv
                area_str = h_w * (t_h - h_conv)
                area = area_conv + area_str

            self.result_text.delete('1.0', tk.END)
            self.result_text.insert(tk.END, f"Target Area: {area:.1f} cm2\n")
            self.result_text.insert(tk.END, f"{'N':<4} | {'Chains (per)':<12}\n")
            self.result_text.insert(tk.END, "-" * 20 + "\n")
            
            counts = []
            for N in sorted(phi_values.keys()):
                phi = phi_values[N]
                total_beads = (area * phi) / b_area
                chains = round(total_beads / N)
                counts.append(chains * n_h)
                self.result_text.insert(tk.END, f"{N:<4} | {chains:<12}\n")
            
            self.result_text.insert(tk.END, f"\nTotal for {n_h} hoppers:\n")
            self.result_text.insert(tk.END, ",".join(map(str, counts)) + "\n")
            self.result_text.insert(tk.END, "\n(Copy to --n_fill)")

        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {e}")

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
        
        # UI State Persistence
        self.config_path = os.path.join(os.path.dirname(__file__), "ui_config.json")
        self.ui_state = self._load_ui_state()

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
        
        # 3. Status Bar for Picking
        self.status_bar = tk.Label(self, text='Click an atom to identify it', bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # 3.5 Restart Editor logic (initialized early for layout)
        def on_re_close():
            # Refresh view after closing editor
            if self.data_ctrl.current_sim_folder:
                self.renderer.show_timestep(self.current_timestep)
            else:
                self.renderer.clear()
            self.update_idletasks()
            self._autosize_and_center()
        
        # 3.6 Simulation Launcher
        self.orchestrator = SimulationOrchestrator(lammps_executable="lmp")

        # 4. Loader Interface and specialized controllers
        loader_cbs = {
            'on_load_success': self._on_data_loaded,
            'get_current_ts': lambda: self.current_timestep,
            'sync_ui_to_frame': self._sync_ui_to_frame,
            'clear_vtk_list': lambda: self.vtk_listbox.delete(0, tk.END),
            'add_vtk': self._add_vtk_mesh_ui,
            'enable_preloading': self.enable_preloading_var,
            'open_restart_editor': self._on_restart_editor_open,
            'ask_column_mapping': self._ask_column_mapping
        }
        self.loader_ctrl = SimulationLoader(self, self.data_ctrl, self.playback_ctrl, self.vtk_ctrl, self.renderer, loader_cbs)
        
        # Link Picking Callback
        self.renderer.on_pick_cb = self._on_atom_picked

        # ── UI LAYOUT ───────────────────────────────────────────
        self.ctrl_panel = tk.Frame(self)
        self.ctrl_panel.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        # Helper to create collapsible columns
        def create_collapsible(parent, name, title_text):
            is_visible = self.ui_state.get(f"{name}_visible", True)
            
            # Container for handle + content
            container = tk.Frame(parent)
            container.pack(side=tk.LEFT, fill=tk.Y)
            
            # Content frame
            content = tk.Frame(container)
            
            # Toggle handle
            handle_frame = tk.Frame(container, width=15, bg='#f0f0f0', relief=tk.FLAT)
            handle_frame.pack(side=tk.LEFT, fill=tk.Y, padx=2)
            
            toggle_btn = tk.Button(handle_frame, text="«" if is_visible else "»", 
                                  font=('Arial', 8), bg='#e0e0e0', relief=tk.FLAT, bd=0,
                                  command=lambda: toggle())
            toggle_btn.pack(side=tk.TOP, fill=tk.X)
            
            # Vertical title on handle (optional but nice)
            tk.Label(handle_frame, text=title_text, font=('Arial', 7), bg='#f0f0f0', fg='#888').pack(side=tk.TOP, pady=10)

            def toggle():
                nonlocal is_visible
                if is_visible:
                    content.pack_forget()
                    toggle_btn.config(text="»")
                    is_visible = False
                else:
                    content.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
                    toggle_btn.config(text="«")
                    is_visible = True
                self.ui_state[f"{name}_visible"] = is_visible
                self._save_ui_state()
                self._autosize_and_center()

            # Initial state
            if is_visible:
                content.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
            
            return content

        # Column 1 is always visible
        col1 = tk.Frame(self.ctrl_panel)
        col1.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))
        
        # Columns 2, 3, 4 are collapsible
        col2 = create_collapsible(self.ctrl_panel, "col2", "EDIT")
        col3 = create_collapsible(self.ctrl_panel, "col3", "LAUNCH")
        col4 = create_collapsible(self.ctrl_panel, "col4", "EXPORT")
        col5 = create_collapsible(self.ctrl_panel, "col5", "CALC")
        
        # 4. Now we can fully init specialized controllers that need UI parents
        self.re_ctrl = RestartEditorController(col2, self.data_ctrl, self.loader_ctrl, self.renderer, on_re_close, self._sync_restart_selection)
        self.sim_launcher = SimulationLauncherController(col3, self.orchestrator, on_simulation_started=lambda: self.refresh_button.invoke())
        self.movie_exporter = MovieExporterController(col4, self.data_ctrl, self.renderer, self.plotter)
        self.hopper_calc = HopperCalculatorController(col5)

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
        self.ts_listbox = tk.Listbox(col1, width=22, height=6, exportselection=False)
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
        tk.Button(playback_frame, text='Cache All', command=self._on_cache_all_ui).pack(side=tk.LEFT, padx=(4, 0))

        self.playback_ctrl.link_widgets(self.frame_slider, self.fps_entry, self.play_btn, self.frame_label, self.loop_var)

        view_btns = tk.Frame(col1)
        view_btns.pack(fill=tk.X, pady=(12, 0))
        tk.Button(view_btns, text='Fit View', command=self.renderer.fit_view).pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(view_btns, text='Reset View', command=self.renderer.reset_view).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

        # ── Column 1 : Lifecycle History ──────────────────────────
        tk.Label(col1, text='Simulation Lifecycle:', fg='#555').pack(anchor='w', pady=(12, 0))
        self.history_listbox = tk.Listbox(col1, width=22, height=4, selectmode=tk.MULTIPLE, exportselection=False, font=('Arial', 8), bg='#fafafa')
        self.history_listbox.pack(fill=tk.Y)
        self.history_listbox.bind('<Double-Button-1>', self._on_history_double_click_ui)
        tk.Button(col1, text='Load Selected Lifecycle', command=self._on_load_lifecycle_selection_ui, font=('Arial', 8)).pack(fill=tk.X, pady=(2, 0))

        # ── Column 2: Mesh & Utilities ──────────────────────────
        tk.Label(col2, text='VTK Meshes / Geometry:').pack(anchor='w', pady=(4, 0))
        vtk_frame = tk.Frame(col2); vtk_frame.pack(fill=tk.BOTH, expand=True)
        self.vtk_listbox = tk.Listbox(vtk_frame, selectmode=tk.MULTIPLE, height=5, exportselection=False)
        self.vtk_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
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
                # Only work if it's a single integer
                val = int(self.highlight_id_var.get().strip())
                self.highlight_id_var.set(str(max(1, val + delta)))
                self._apply_highlight_ui()
            except ValueError:
                # If it's a range like "1-3,5", incrementing doesn't make sense
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

        tk.Label(col2, text='──── Restart Editor ────', fg='gray').pack(fill=tk.X, pady=(12, 2))
        tk.Button(col2, text='Edit Restart File', command=lambda: self.loader_ctrl.open_restart_editor(filedialog.askopenfilename(filetypes=[('Restart Files', '*.bin')]))).pack(fill=tk.X)
        # Note: self.re_ctrl.panel is ready and will be packed by re_ctrl.open()

        # 4. Finalize UI & Show Plotter
        self.plotter.show(interactive=False)
        self.update_idletasks()
        self._autosize_and_center()
        
        # Keybinds
        self.bind('<space>', lambda e: self.playback_ctrl.toggle_play())
        self.bind('<Left>', lambda e: self.playback_ctrl.step_prev())
        self.bind('<Right>', lambda e: self.playback_ctrl.step_next())

    # ── Internal Helpers ─────────────────────────────────────
    def _get_ui_settings(self):
        return {'show_geometry': self.show_geometry_var.get()}

    def _get_ui_vars(self):
        return {'chain_size': self.chain_size_var, 'enable_preloading': self.enable_preloading_var}

    def _get_current_rendering_data(self):
        """Standardizes data retrieval for highlighter/export during both normal and preview modes."""
        ts = self.current_timestep
        # For simulation folders, we need to get the actual frame data for the current TS
        if self.data_ctrl.sim_source:
            df = self.data_ctrl.get_atom_data_at_timestep(ts)
            return df, ts
        return self.data_ctrl.df_mi, ts

    def _get_chain_size(self):
        try: return int(self.chain_size_var.get())
        except: return 4

    @property
    def current_timestep(self):
        # Safety for call during initialization
        if not hasattr(self, 'playback_ctrl') or self.playback_ctrl is None:
            return 0
            
        # In preview mode, we treat everything as step 0
        if self.data_ctrl.is_preview_mode:
            return 0
        try:
            return int(self.data_ctrl.timesteps[self.playback_ctrl.current_frame])
        except (AttributeError, IndexError, TypeError):
            return 0

    def _on_data_loaded(self, num_queued=0, clear_vtk=True):
        # Handle progress bar if we have background batches to wait for
        if num_queued > 0:
            self.batches_to_load = num_queued
            self.batches_loaded = 0
            self.progress_bar = ProgressBar(self, title="Preloading Cache", message=f"Loading {num_queued} cached batches...", modal=True)
            self.progress_bar.set_progress(0, num_queued)
        else:
            self.batches_to_load = 0
            
        # Ensure fresh VTK state for new folder metadata (optional)
        if clear_vtk:
            self.vtk_ctrl.metadata = None
            self.vtk_ctrl.clear()
            self.vtk_listbox.delete(0, tk.END)

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
            self.renderer.show_timestep(ts)
            
        # NEW: Restore metadata-driven state from simulation folder
        if self.data_ctrl.metadata:
            # NEW: Restore from full history, but save only to current run
            self.vtk_ctrl.restore_from_metadata(self.data_ctrl.metadata)
            self.vtk_ctrl.metadata = self.data_ctrl.writable_metadata
            
        # 2. Update Lifecycle History
        self.history_listbox.delete(0, tk.END)
        if hasattr(self.data_ctrl.sim_source, 'lineage'):
            lineage = self.data_ctrl.sim_source.lineage
            # Determine which folders are currently part of the active view
            active_paths = []
            if hasattr(self.data_ctrl.sim_source, 'data_dirs'):
                active_paths = [os.path.abspath(p) for p in self.data_ctrl.sim_source.data_dirs]
            else:
                active_paths = [os.path.abspath(self.data_ctrl.current_sim_folder)]

            for i, d in enumerate(lineage):
                name = os.path.basename(d)
                is_primary = os.path.abspath(d) == os.path.abspath(self.data_ctrl.current_sim_folder)
                
                if is_primary:
                    name = f"➤ {name}"
                
                self.history_listbox.insert(tk.END, name)
                
                # Restore selection for all folders in the active composite set
                if os.path.abspath(d) in active_paths:
                    self.history_listbox.selection_set(i)
                    if is_primary:
                        self.history_listbox.see(i)
            
            # 2. Sync the VTK listbox in UI with the restored meshes and visibility
            self.vtk_listbox.delete(0, tk.END)
            for i, (name, data) in enumerate(self.vtk_ctrl.vtk_meshes.items()):
                self.vtk_listbox.insert(tk.END, name)
                if data['visible']:
                    self.vtk_listbox.selection_set(i)
                else:
                    self.vtk_listbox.selection_clear(i)
            
            # 3. Restore Camera state
            self.renderer.restore_camera_state()

        # 4. Sync Movie Exporter
        self.movie_exporter._set_range_from_loader()

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
            self.renderer.show_timestep(ts)
            self.analysis_ctrl.refresh_windows()

    def _sync_ui_to_frame(self, idx, ts):
        self.frame_slider.set(idx)
        self.ts_listbox.selection_clear(0, tk.END)
        self.ts_listbox.selection_set(idx)
        self.ts_listbox.activate(idx)
        self.renderer.show_timestep(ts)

    def _add_vtk_mesh_ui(self, path):
        ok, res = self.vtk_ctrl.add_mesh(path)
        if ok:
            self.vtk_listbox.insert(tk.END, res)
            self.vtk_listbox.selection_set(self.vtk_listbox.size() - 1)
            self.renderer.update_persistent_bounds()
            if self.current_timestep is not None: self.renderer.show_timestep(self.current_timestep)
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
            self.frame_slider.set(idx); self.renderer.show_timestep(ts)

    def _on_slider_ui(self, val):
        idx = int(float(val))
        if 0 <= idx < len(self.data_ctrl.timesteps):
            ts = self.data_ctrl.timesteps[idx]
            self.ts_listbox.selection_clear(0, tk.END); self.ts_listbox.selection_set(idx); self.ts_listbox.activate(idx)
            self.playback_ctrl.update_status(idx); self.renderer.show_timestep(ts)

    def _on_vtk_select_ui(self, event):
        sel = self.vtk_listbox.curselection()
        for i in range(self.vtk_listbox.size()):
            self.vtk_ctrl.set_visibility(self.vtk_listbox.get(i), i in sel)
        if self.current_timestep is not None: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def _delete_selected_vtks_ui(self):
        sel = self.vtk_listbox.curselection()
        if not sel: return
        names = [self.vtk_listbox.get(i) for i in sel]
        self.vtk_ctrl.remove_meshes(names)
        for i in reversed(sel): self.vtk_listbox.delete(i)
        self.renderer.update_persistent_bounds()
        if self.current_timestep is not None: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def _apply_highlight_ui(self):
        try:
            # Pass the raw string to support ranges like "1-4,6"
            spec = self.highlight_id_var.get().strip()
            if not spec: return
            
            ok, res = self.hl_ctrl.apply(self.highlight_mode_var.get(), spec)
            self.highlight_status.config(text=res, fg='black' if ok else 'red')
            if ok:
                self.renderer.show_timestep(self.current_timestep)
            else:
                self.plotter.render()
        except Exception as e:
            self.highlight_status.config(text=f'Highlight Error: {str(e)}', fg='red')
            import traceback
            traceback.print_exc()

    def _clear_highlight_ui(self):
        self.hl_ctrl.clear(); self.highlight_status.config(text='—', fg='gray'); self.plotter.render()

    def _on_draw_toggle_ui(self):
        self.vtk_ctrl.set_master_visibility(self.show_geometry_var.get())
        if self.current_timestep is not None: self.renderer.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def _on_cache_all_ui(self):
        if not self.data_ctrl.sim_source:
            return
            
        new_queued = self.data_ctrl.cache_all_remaining()
        if new_queued > 0:
            if self.progress_bar is None:
                self.batches_to_load = new_queued
                self.batches_loaded = 0
                self.progress_bar = ProgressBar(self, title="Caching", message=f"Loading {new_queued} cached batches...", modal=True)
                self.progress_bar.set_progress(0, new_queued)
            else:
                self.batches_to_load += new_queued
                self.progress_bar.set_progress(self.batches_loaded, self.batches_to_load)

    # ── Restart Editor Logic ──────────────────────────────
    def _ask_column_mapping(self, filepath):
        dialog = ColumnMapDialog(self, filepath)
        self.wait_window(dialog)
        return dialog.result

    def _on_atom_picked(self, info: dict):
        """Callback from renderer when an atom is clicked."""
        msg = f"Atom ID: {info['id']} | Mol/Chain: {info['mol']} | Type: {info['type']} | Pos: {info['pos']}"
        self.status_bar.config(text=msg, fg='blue')
        
        # In Restart Editor mode, toggle selection
        if self.data_ctrl.is_preview_mode:
            self.re_ctrl.toggle_molecule(info['mol'])
        else:
            # Traditional behavior: just log and copy
            self.clipboard_clear()
            self.clipboard_append(str(info['mol']))

    def _on_history_double_click_ui(self, event):
        idx = self.history_listbox.curselection()
        if not idx: return
        
        # Resolve the directory from the index
        if hasattr(self.data_ctrl.sim_source, 'lineage'):
            target_dir = self.data_ctrl.sim_source.lineage[idx[0]]
            if os.path.isdir(target_dir):
                self.loader_ctrl.load_simulation_folder(target_dir)
            else:
                messagebox.showerror("Error", f"Folder no longer exists: {target_dir}")

    def _on_load_lifecycle_selection_ui(self):
        """Loads multiple selected folders as a single composite simulation."""
        indices = self.history_listbox.curselection()
        if not indices: return
        
        if hasattr(self.data_ctrl.sim_source, 'lineage'):
            lineage = self.data_ctrl.sim_source.lineage
            # Get all selected paths in chronological order (by lineage index)
            target_dirs = [lineage[i] for i in sorted(indices)]
            
            # Filter valid directories
            valid_dirs = [d for d in target_dirs if os.path.isdir(d)]
            if not valid_dirs:
                messagebox.showerror("Error", "None of the selected folders exist.")
                return
                
            self.loader_ctrl.load_simulation_folder(valid_dirs)

    def _on_restart_editor_open(self, path: str):
        """Callback when a restart file is loaded for editing."""
        # Clear highlights from any previous work
        self._clear_highlight_ui()
        self.re_ctrl.open(path)
        self.update_idletasks()
        self._autosize_and_center()

    def _sync_restart_selection(self, ids_str: str):
        """Syncs the selection from Restart Editor to the main Highlighter."""
        self.highlight_mode_var.set('chain')
        self.highlight_id_var.set(ids_str)
        self._apply_highlight_ui()


    def _autosize_and_center(self) -> None:
        try:
            self.update_idletasks()
            # Use requested dimensions to allow shrinking when columns collapse
            w = self.winfo_reqwidth()
            h = self.winfo_reqheight()
            
            # Get current position to stay in the same place while resizing
            geom = self.geometry().split('+')
            if len(geom) == 3:
                x, y = geom[1], geom[2]
                self.geometry(f"{w}x{h}+{x}+{y}")
            else:
                sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
                x, y = max(0, int((sw - w) / 2)), max(0, int((sh - h) / 2))
                self.geometry(f"{w}x{h}+{x}+{y}")
            
            # Update minsize to allow the new smaller size
            self.minsize(200, 200) 
        except: pass

    def _load_ui_state(self) -> dict:
        import json
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as f:
                    return json.load(f)
            except: pass
        return {}

    def _save_ui_state(self):
        import json
        try:
            with open(self.config_path, 'w') as f:
                json.dump(self.ui_state, f)
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
    parser = argparse.ArgumentParser(description="Chains Simulation Viewer")
    parser.add_argument("sim_path", nargs="?", default=None, help="Path to simulation folder or file to load")
    args, unknown = parser.parse_known_args()

    app = ViewerApp()
    if args.sim_path:
        app.after(100, lambda: app.loader_ctrl.handle_dropped_files([args.sim_path]))
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
