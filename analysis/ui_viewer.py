import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional
import pandas as pd
import numpy as np
from vedo import Plotter, Spheres, Lines, Axes, Box, Cylinder, Cone, Plane, Mesh
from analysis.utilities import validate_chain_spacing
from analysis.ts_windows import (BondPlotWindow, AnglePlotWindow, AtomPlotWindow,
                                   parse_range_spec)
from analysis.data_manager import SimulationData, parse_simple_data_file, load_lammps_geometry

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

class SimDataController:
    """Manages loading, parsing, and caching of simulation data."""
    def __init__(self):
        self.df: Optional[pd.DataFrame] = None
        self.df_mi: Optional[pd.DataFrame] = None
        self.df_bonds: Optional[pd.DataFrame] = None
        self.df_angles: Optional[pd.DataFrame] = None
        self.current_sim_folder: Optional[str] = None
        self.geometry_data: Optional[dict] = None
        self.geometry_script_path: Optional[str] = None
        self.timesteps: list = []
        self._init_limits: Optional[tuple] = None
        self.on_data_loaded_cb = None

    def load_folder(self, folder: str, force_reload: bool = False):
        try:
            sim = SimulationData(folder)
            data_dict = sim.load_data(force_reload=force_reload)
            df_atoms = data_dict['atoms']
            if df_atoms is not None and not df_atoms.empty:
                df_atoms = df_atoms.reset_index()
            
            if df_atoms is None or df_atoms.empty:
                return False, "No dump files found or parsing failed."
            
            # Reset index for common processing if needed, 
            # though load_dataframe will handle the atoms DF.
            self.current_sim_folder = folder
            self.geometry_data, self.geometry_script_path = load_lammps_geometry(folder)
            
            # Multi-dataset state
            self.df_bonds = data_dict['bonds']
            self.df_angles = data_dict['angles']
            
            self.load_dataframe(df_atoms)
            return True, None
        except Exception as e:
            import traceback
            traceback.print_exc()
            return False, str(e)

    def load_dump_files(self, files: list):
        frames = []
        for f in files:
            try:
                frame = self._parse_single_dump(f)
                if frame is not None:
                    frames.append(frame)
            except Exception as e:
                print(f"Failed to parse {f}: {e}")
        
        if not frames:
            return False, "No valid dump frames parsed."
        
        full = pd.concat(frames).reset_index()
        self.current_sim_folder = None
        self.geometry_data = None
        self.geometry_script_path = None
        self.load_dataframe(full)
        return True, None

    def load_data_file(self, path: str):
        try:
            df = parse_simple_data_file(path)
            if df.empty:
                return False, "Could not parse the selected data file"
            df = df.reset_index()
            self.current_sim_folder = None
            self.geometry_data = None
            self.geometry_script_path = None
            self.load_dataframe(df)
            return True, None
        except Exception as e:
            return False, str(e)

    def load_dataframe(self, df: pd.DataFrame):
        self.df = df.copy()
        # Ensure critical columns are numeric for calculation
        for c in ['x', 'y', 'z', 'id', 'timestep']:
            if c in self.df.columns:
                self.df[c] = pd.to_numeric(self.df[c], errors='coerce')
            else:
                raise ValueError(f"Missing column: {c}")
        
        # Provide fallback diameter if missing
        if 'diameter' not in self.df.columns:
            self.df['diameter'] = 0.01
        else:
            self.df['diameter'] = pd.to_numeric(self.df['diameter'], errors='coerce').fillna(0.01)

        self.df.dropna(subset=['x', 'y', 'z', 'id', 'timestep'], inplace=True)

        x_min, x_max = self.df['x'].min(), self.df['x'].max()
        y_min, y_max = self.df['y'].min(), self.df['y'].max()
        z_min, z_max = self.df['z'].min(), self.df['z'].max()
        self._init_limits = ((x_min, x_max), (y_min, y_max), (z_min, z_max))

        try:
            self.df_mi = self.df.set_index(['timestep', 'id']).sort_index()
        except Exception:
            self.df_mi = None

        self.timesteps = sorted(self.df['timestep'].unique())
        if self.on_data_loaded_cb:
            self.on_data_loaded_cb()

    def _parse_single_dump(self, f: str):
        try:
            with open(f, 'r') as fh:
                lines = fh.readlines()
            timestep = 0
            for i, line in enumerate(lines):
                if 'ITEM: TIMESTEP' in line:
                    timestep = int(lines[i+1])
                if 'ITEM: ATOMS' in line:
                    cols = line.split()[2:]
                    df = pd.read_csv(f, skiprows=i+1, names=cols, sep=r'\s+', engine='python')
                    df['timestep'] = timestep
                    return df.set_index(['timestep', 'id'])
        except Exception:
            return None
        return None

    def _normalize_dropped_path(self, p: str) -> str:
        p = (p or "").strip()
        if p.startswith("{") and p.endswith("}"):
            p = p[1:-1]
        try:
            return os.path.normpath(p)
        except Exception:
            return p

    def _resolve_sim_root(self, dropped_dir: str) -> Optional[str]:
        d = self._normalize_dropped_path(dropped_dir)
        if not d or not os.path.isdir(d):
            return None
        if os.path.isdir(os.path.join(d, "chain")):
            return d
        if os.path.basename(d).lower() == "chain":
            parent = os.path.dirname(d)
            if parent and os.path.isdir(os.path.join(parent, "chain")):
                return parent
        return None


class VtkOverlayController:
    """Manages VTK mesh loading, visibility, and scene bounds."""
    def __init__(self, plotter):
        self.plotter = plotter
        self.vtk_meshes: dict = {}
        self.vtk_color_idx: int = 0
        self._scene_bounds: Optional[list] = None

    def add_mesh(self, path: str):
        name = os.path.basename(path)
        if name in self.vtk_meshes:
            return False, "Already loaded"
        try:
            import vedo
            obj = vedo.load(path)
            mesh = obj.tomesh() if hasattr(obj, "tomesh") else obj
            colors = ['red', 'green', 'blue', 'gold', 'cyan', 'magenta', 'orange', 'purple', 'lime', 'pink']
            color = colors[self.vtk_color_idx % len(colors)]
            self.vtk_color_idx += 1
            mesh.c(color).alpha(0.8)
            self.vtk_meshes[name] = {
                'path': path, 'actor': mesh, 'visible': True, 'color': color
            }
            self.plotter.add(mesh)
            mesh.on()
            return True, name
        except Exception as e:
            return False, str(e)

    def set_visibility(self, name: str, visible: bool):
        if name in self.vtk_meshes:
            self.vtk_meshes[name]['visible'] = visible
            self.sync_visibility()

    def sync_visibility(self):
        for mesh_data in self.vtk_meshes.values():
            act = mesh_data['actor']
            if mesh_data['visible']:
                if act not in self.plotter.actors:
                    self.plotter.add(act)
                act.on()
            else:
                if act in self.plotter.actors:
                    self.plotter.remove(act)
                act.off()
        self.plotter.render()

    def remove_meshes(self, names: list):
        for name in names:
            if name in self.vtk_meshes:
                act = self.vtk_meshes[name]['actor']
                try:
                    self.plotter.remove(act)
                except Exception: pass
                del self.vtk_meshes[name]

    def clear(self):
        for mesh_data in self.vtk_meshes.values():
            if mesh_data['actor'] in self.plotter.actors:
                self.plotter.remove(mesh_data['actor'])
        self.vtk_meshes.clear()
        self.vtk_color_idx = 0
        self._scene_bounds = None

    def recompute_bounds(self, init_limits: Optional[tuple], only_visible=True):
        if init_limits is not None:
            xmin, xmax = init_limits[0]
            ymin, ymax = init_limits[1]
            zmin, zmax = init_limits[2]
        else:
            xmin, xmax = float('inf'), float('-inf')
            ymin, ymax = float('inf'), float('-inf')
            zmin, zmax = float('inf'), float('-inf')

        for mesh_data in self.vtk_meshes.values():
            if not only_visible or mesh_data['visible']:
                try:
                    bnds = mesh_data['actor'].bounds()
                    if len(bnds) == 6:
                        xmin = min(xmin, bnds[0]); xmax = max(xmax, bnds[1])
                        ymin = min(ymin, bnds[2]); ymax = max(ymax, bnds[3])
                        zmin = min(zmin, bnds[4]); zmax = max(zmax, bnds[5])
                except Exception: pass

        if xmin == float('inf'):
            xmin, xmax, ymin, ymax, zmin, zmax = -1, 1, -1, 1, -1, 1

        x_pad = (xmax - xmin) * 0.05 if xmax > xmin else 0.1
        y_pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.1
        z_pad = (zmax - zmin) * 0.05 if zmax > zmin else 0.1
        
        self._scene_bounds = [
            xmin - x_pad, xmax + x_pad,
            ymin - y_pad, ymax + y_pad,
            zmin - z_pad, zmax + z_pad
        ]
        return self._scene_bounds


class HighlightController:
    """Manages highlight state and rendering of highlighted atoms."""
    _HL_COLORS = {'atom': 'yellow', 'bond': 'orange', 'angle': 'violet', 'chain': 'lime'}
    
    def __init__(self, plotter, get_data_cb, get_cs_cb):
        self.plotter = plotter
        self.get_data_cb = get_data_cb
        self.get_cs_cb = get_cs_cb
        self._highlighted_ids: set = set()
        self._highlight_color: str = 'yellow'
        self._highlight_actors: list = []

    def clear(self):
        self._highlighted_ids = set()
        for act in self._highlight_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._highlight_actors = []

    def apply(self, mode: str, n: int):
        df, ts = self.get_data_cb()
        if df is None or ts is None:
            return False, "Data not available"
        
        subset = df[df['timestep'] == ts]
        cs = self.get_cs_cb()
        
        if mode == 'atom':
            ids, label = {n}, f'Atom {n}'
        elif mode == 'bond':
            nb = cs - 1
            if nb <= 0: return False, "Chain size must be ≥ 2"
            ci, binc = (n-1)//nb, (n-1)%nb
            cf = ci * cs + 1
            a1, a2 = cf + binc, cf + binc + 1
            ids, label = {a1, a2}, f'Bond {n} (atoms {a1}, {a2})'
        elif mode == 'angle':
            na = cs - 2
            if na <= 0: return False, "Chain size must be ≥ 3"
            ci, ainc = (n-1)//na, (n-1)%na
            cf = ci * cs + 1
            a1, a2, a3 = cf + ainc, cf + ainc + 1, cf + ainc + 2
            ids, label = {a1, a2, a3}, f'Angle {n} (atoms {a1}, {a2}, {a3})'
        elif mode == 'chain':
            if 'mol' in subset.columns:
                ids = set(subset[subset['mol'] == n]['id'].values.tolist())
            else:
                start = (n-1)*cs + 1
                ids = set(range(start, start + cs))
            label = f'Chain {n} ({len(ids)} atoms)'
        else: return False, "Unknown mode"

        self._highlighted_ids = ids
        self._highlight_color = self._HL_COLORS.get(mode, 'yellow')
        found = ids & set(subset['id'].values)
        if not found:
            return False, f'{label} - no atoms in frame'
        
        self.render(df, ts)
        return True, label

    def render(self, df, ts):
        for act in self._highlight_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._highlight_actors = []

        if not self._highlighted_ids or df is None or ts is None:
            return
        
        subset = df[df['timestep'] == ts]
        hi = subset[subset['id'].isin(self._highlighted_ids)]
        if hi.empty: return

        pts = hi[['x', 'y', 'z']].values
        r = (hi['diameter'].values / 2.0) * 1.30
        actor = Spheres(pts, r=r, c=self._highlight_color, alpha=0.95)
        self._highlight_actors = [actor]
        self.plotter.add(actor)

class PlaybackController:
    """Manages playback state, timer, and transport controls."""
    def __init__(self, root, get_timesteps_cb, on_frame_change_cb):
        self.root = root
        self.get_timesteps_cb = get_timesteps_cb
        self.on_frame_change_cb = on_frame_change_cb
        
        self.playing = False
        self.job = None
        
        # UI References
        self.frame_slider = None
        self.fps_slider = None
        self.play_btn = None
        self.frame_label = None
        self.loop_var = None

    def link_widgets(self, frame_slider, fps_slider, play_btn, frame_label, loop_var):
        self.frame_slider = frame_slider
        self.fps_slider = fps_slider
        self.play_btn = play_btn
        self.frame_label = frame_label
        self.loop_var = loop_var

    def update_status(self, idx: int):
        if not self.frame_label: return
        timesteps = self.get_timesteps_cb()
        total = len(timesteps)
        if total == 0:
            self.frame_label.config(text="No frames loaded")
            return
        ts = timesteps[idx]
        self.frame_label.config(text=f"Frame: {idx+1} / {total} (TS: {ts})")

    def jump_start(self):
        self.pause()
        if self.frame_slider and self.get_timesteps_cb():
            self.frame_slider.set(0)

    def jump_end(self):
        self.pause()
        ts = self.get_timesteps_cb()
        if self.frame_slider and ts:
            self.frame_slider.set(len(ts) - 1)

    def step_next(self):
        self.pause()
        if not self.frame_slider: return
        curr = int(self.frame_slider.get())
        ts = self.get_timesteps_cb()
        if curr < len(ts) - 1:
            self.frame_slider.set(curr + 1)
        elif self.loop_var and self.loop_var.get():
            self.frame_slider.set(0)

    def step_prev(self):
        self.pause()
        if not self.frame_slider: return
        curr = int(self.frame_slider.get())
        ts = self.get_timesteps_cb()
        if curr > 0:
            self.frame_slider.set(curr - 1)
        elif self.loop_var and self.loop_var.get():
            self.frame_slider.set(len(ts)-1)

    def pause(self):
        self.playing = False
        if self.job:
            self.root.after_cancel(self.job)
            self.job = None
        if self.play_btn:
            self.play_btn.config(text="▶ Play")

    def toggle_play(self):
        if self.playing:
            self.pause()
        else:
            if not self.get_timesteps_cb(): return
            self.playing = True
            if self.play_btn:
                self.play_btn.config(text="⏸ Pause")
            self._tick()

    def _tick(self):
        if not self.playing or not self.frame_slider: return
        
        curr = int(self.frame_slider.get())
        ts = self.get_timesteps_cb()
        total = len(ts)
        
        next_idx = curr + 1
        if next_idx >= total:
            if self.loop_var and self.loop_var.get():
                next_idx = 0
            else:
                self.pause()
                return

        self.frame_slider.set(next_idx)
        
        fps = self.fps_slider.get() if self.fps_slider else 10
        delay = int(1000 / max(fps, 1))
        self.job = self.root.after(delay, self._tick)

class ViewerApp(BaseTk):
    def __init__(self):
        super().__init__()
        self.title('Chains Simulation Viewer')
        # Let the window size to its content, then center.
        # (Avoid a fixed width that leaves empty space to the right.)
        # Ensure clean shutdown when window is closed
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        # Controllers (Pre-initialize logic controllers for UI binding)
        self.data_ctrl = SimDataController()
        self.data_ctrl.on_data_loaded_cb = self._on_data_loaded
        self.playback_ctrl = PlaybackController(self, lambda: self.data_ctrl.timesteps, self.show_timestep)

        # Left controls
        ctrl = tk.Frame(self)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        # Two side-by-side columns inside ctrl
        col1 = tk.Frame(ctrl)
        col1.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 4))

        try:
            import tkinter.ttk as ttk
            ttk.Separator(ctrl, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=4)
        except Exception:
            tk.Frame(ctrl, width=1, bg='gray').pack(side=tk.LEFT, fill=tk.Y, padx=4)

        col2 = tk.Frame(ctrl)
        col2.pack(side=tk.LEFT, fill=tk.Y)

        # ── Column 1 : File loading & navigation ─────────────────
        tk.Button(col1, text='Open Dump Folder', command=self.open_dump_folder).pack(fill=tk.X)
        tk.Button(col1, text='Open Dump Files...', command=self.open_dump_files).pack(fill=tk.X, pady=(4, 0))
        tk.Button(col1, text='Open Data File...', command=self.open_data_file).pack(fill=tk.X, pady=(4, 0))
        tk.Button(col1, text='Open VTK Files...', command=self.open_vtk_files).pack(fill=tk.X, pady=(4, 0))
        self.refresh_button = tk.Button(col1, text='Refresh', command=self.refresh_current_source, state=tk.DISABLED)
        self.refresh_button.pack(fill=tk.X, pady=(8, 0))

        # Drag-and-drop area
        drop_text = ('Drop dump/.data files here' if DND_AVAILABLE
                     else 'Drag-and-drop disabled\n(install tkinterdnd2)')
        self.drop_label = tk.Label(col1, text=drop_text, relief='ridge', width=22, height=4)
        self.drop_label.pack(fill=tk.X, pady=(8, 4))
        if DND_AVAILABLE:
            try:
                self.drop_label.drop_target_register(DND_CONST)
                self.drop_label.dnd_bind('<<Drop>>', self.on_drop)
            except Exception:
                self.drop_label.config(text='Drag-and-drop not available on this platform')

        tk.Label(col1, text='Timesteps:').pack(anchor='w', pady=(8, 0))
        self.ts_listbox = tk.Listbox(col1, width=22, height=6)
        self.ts_listbox.pack(fill=tk.Y)
        self.ts_listbox.bind('<<ListboxSelect>>', self.on_ts_select)

        tk.Label(col1, text='Frame Navigation:').pack(anchor='w', pady=(8, 0))
        
        # New Playback Control Frame
        playback_frame = tk.Frame(col1)
        playback_frame.pack(fill=tk.X, pady=(2, 0))

        btn_row = tk.Frame(playback_frame)
        btn_row.pack(fill=tk.X)
        
        tk.Button(btn_row, text='|◀', width=3, command=self.playback_ctrl.jump_start).pack(side=tk.LEFT)
        tk.Button(btn_row, text='◀', width=3, command=self.playback_ctrl.step_prev).pack(side=tk.LEFT, padx=2)
        self.play_btn = tk.Button(btn_row, text='▶ Play', width=8, command=self.playback_ctrl.toggle_play)
        self.play_btn.pack(side=tk.LEFT, padx=2)
        tk.Button(btn_row, text='▶', width=3, command=self.playback_ctrl.step_next).pack(side=tk.LEFT, padx=2)
        tk.Button(btn_row, text='▶|', width=3, command=self.playback_ctrl.jump_end).pack(side=tk.LEFT)

        self.frame_slider = tk.Scale(playback_frame, from_=0, to=0, orient=tk.HORIZONTAL, command=self.on_slider)
        self.frame_slider.pack(fill=tk.X, pady=(4, 0))

        self.frame_label = tk.Label(playback_frame, text='Frame: 0 / 0', fg='gray')
        self.frame_label.pack(fill=tk.X)

        speed_frame = tk.Frame(playback_frame)
        speed_frame.pack(fill=tk.X, pady=(4, 0))
        tk.Label(speed_frame, text='FPS:').pack(side=tk.LEFT)
        self.fps_slider = tk.Scale(speed_frame, from_=1, to=60, orient=tk.HORIZONTAL)
        self.fps_slider.set(10)
        self.fps_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))
        
        self.loop_var = tk.BooleanVar(value=True)
        tk.Checkbutton(playback_frame, text='Loop Playback', variable=self.loop_var).pack(anchor='w')

        # Link Playback Controller to its widgets
        self.playback_ctrl.link_widgets(self.frame_slider, self.fps_slider, self.play_btn, self.frame_label, self.loop_var)

        tk.Button(col1, text='Fit View', command=self.fit_view).pack(fill=tk.X, pady=(12, 0))
        tk.Button(col1, text='Reset View', command=self.reset_view).pack(fill=tk.X, pady=(4, 0))
        self.draw_chains_var = tk.BooleanVar(value=False)
        tk.Checkbutton(col1, text='Draw Chains', variable=self.draw_chains_var,
                       command=self._on_draw_toggle).pack(fill=tk.X, pady=(4, 0))
        self.show_geometry_var = tk.BooleanVar(value=False)
        tk.Checkbutton(col1, text='Show Geometry', variable=self.show_geometry_var,
                       command=self._on_draw_toggle).pack(fill=tk.X, pady=(4, 0))

        # ── Column 2 : VTK management & Highlight ────────────────
        tk.Label(col2, text='Loaded VTKs:', anchor='w').pack(fill=tk.X, pady=(0, 2))
        vtk_frame = tk.Frame(col2)
        vtk_frame.pack(fill=tk.X)
        self.vtk_listbox = tk.Listbox(vtk_frame, width=22, height=5, selectmode=tk.MULTIPLE,
                                      exportselection=False)
        self.vtk_listbox.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.vtk_listbox.bind('<<ListboxSelect>>', self.on_vtk_select)
        vtk_scroll = tk.Scrollbar(vtk_frame, orient=tk.VERTICAL)
        vtk_scroll.config(command=self.vtk_listbox.yview)
        vtk_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.vtk_listbox.config(yscrollcommand=vtk_scroll.set)
        tk.Button(col2, text='Delete Selected VTKs',
                  command=self.delete_selected_vtks).pack(fill=tk.X, pady=(2, 0))

        # Highlight sub-section
        tk.Label(col2, text='──── Highlight ────', fg='gray').pack(fill=tk.X, pady=(12, 2))

        hl_mode_frame = tk.Frame(col2)
        hl_mode_frame.pack(fill=tk.X)
        self.highlight_mode_var = tk.StringVar(value='atom')
        for _mode, _lbl in [('atom', 'Atom'), ('bond', 'Bond'),
                             ('angle', 'Angle'), ('chain', 'Chain')]:
            tk.Radiobutton(hl_mode_frame, text=_lbl,
                           variable=self.highlight_mode_var,
                           value=_mode).pack(side=tk.LEFT)

        hl_id_frame = tk.Frame(col2)
        hl_id_frame.pack(fill=tk.X, pady=(4, 0))
        tk.Label(hl_id_frame, text='ID:').pack(side=tk.LEFT)
        self.highlight_id_var = tk.StringVar()
        _hl_entry = tk.Entry(hl_id_frame, textvariable=self.highlight_id_var, width=7)
        _hl_entry.pack(side=tk.LEFT, padx=(4, 0))
        _hl_entry.bind('<Return>', lambda e: self._apply_highlight())
        tk.Button(hl_id_frame, text='Apply',
                  command=self._apply_highlight).pack(side=tk.LEFT, padx=(4, 0))

        hl_cs_frame = tk.Frame(col2)
        hl_cs_frame.pack(fill=tk.X, pady=(2, 0))
        tk.Label(hl_cs_frame, text='Chain size:').pack(side=tk.LEFT)
        self.chain_size_var = tk.StringVar(value='4')
        tk.Entry(hl_cs_frame, textvariable=self.chain_size_var,
                 width=4).pack(side=tk.LEFT, padx=(4, 0))

        tk.Button(col2, text='Clear Highlight',
                  command=self._clear_highlight).pack(fill=tk.X, pady=(6, 0))
        self.highlight_status = tk.Label(col2, text='—', fg='gray',
                                         wraplength=160, justify='left')
        self.highlight_status.pack(anchor='w', pady=(4, 0))

        # ── Time Series section ───────────────────────────────────
        tk.Label(col2, text='──── Time Series ────', fg='gray').pack(fill=tk.X, pady=(12, 2))
        tk.Button(col2, text='Bond Distances…',
                  command=self._open_bond_win).pack(fill=tk.X)
        tk.Button(col2, text='Angle Values…',
                  command=self._open_angle_win).pack(fill=tk.X, pady=(4, 0))
        tk.Button(col2, text='Atom Properties…',
                  command=self._open_atom_win).pack(fill=tk.X, pady=(4, 0))

        self.plotter = Plotter(
            bg='white',
            interactive=True
        )

        # Controllers
        self.vtk_ctrl = VtkOverlayController(self.plotter)
        self.hl_ctrl = HighlightController(self.plotter, self._get_current_rendering_data, self._get_chain_size)

        # Persistent View Bounds: Stores the max extents of all loaded objects to keep axes fixed
        self._persistent_view_bounds = None

        # Size-to-content and center on screen after widgets are laid out
        self._autosize_and_center()

        # Data holders
        self._dynamic_actors = []
        self.current_timestep = None
        # Store initial axis limits and default view angles
        self._default_view = (30, -60)
        try:
            self.ax.view_init(elev=self._default_view[0], azim=self._default_view[1])
        except Exception:
            pass

        # Time-series window references (None when closed)
        self._bond_win  = None
        self._angle_win = None
        self._atom_win  = None
        # Saved range strings — restored when windows are reopened
        self._bond_range_spec  = ''
        self._angle_range_spec = ''
        self._atom_range_spec  = ''

        # Global Keyboard Bindings
        self.bind('<space>', lambda e: self.playback_ctrl.toggle_play())
        self.bind('<Left>', lambda e: self.playback_ctrl.step_prev())
        self.bind('<Right>', lambda e: self.playback_ctrl.step_next())

    def _get_current_rendering_data(self):
        return self.data_ctrl.df, self.current_timestep

    def _get_chain_size(self):
        try: return int(self.chain_size_var.get())
        except: return 4

    def _on_data_loaded(self):
        # Notify any open time-series windows of new data
        if self._bond_win and self._bond_win.winfo_exists():
            self._bond_win.refresh_df(self.data_ctrl.df_bonds)
        if self._angle_win and self._angle_win.winfo_exists():
            self._angle_win.refresh_df(self.data_ctrl.df_angles)
        if self._atom_win and self._atom_win.winfo_exists():
            self._atom_win.refresh_df(self.data_ctrl.df_mi)

        self.ts_listbox.delete(0, tk.END)
        for t in self.data_ctrl.timesteps:
            self.ts_listbox.insert(tk.END, str(t))

        self.frame_slider.config(from_=0, to=max(0, len(self.data_ctrl.timesteps)-1))
        self.frame_slider.set(0)
        self.playback_ctrl.update_status(0)
        self._update_persistent_bounds()
        if self.data_ctrl.timesteps:
            self.show_timestep(self.data_ctrl.timesteps[0])
            self.plotter.reset_camera()
            self.plotter.render()
    def _autosize_and_center(self) -> None:
        # Ensure geometry requests are computed
        try:
            self.update_idletasks()

            w = max(self.winfo_reqwidth(), self.winfo_width())
            h = max(self.winfo_reqheight(), self.winfo_height())

            sw = self.winfo_screenwidth()
            sh = self.winfo_screenheight()

            x = max(0, int((sw - w) / 2))
            y = max(0, int((sh - h) / 2))

            self.geometry(f"{w}x{h}+{x}+{y}")
            self.minsize(w, h)
        except Exception:
            # If anything goes wrong (platform quirks), keep default behavior.
            pass

    def _on_draw_toggle(self):
        # refresh current timestep view to show/hide chains
        if self.current_timestep is not None:
            try:
                self.show_timestep(self.current_timestep)
            except Exception:
                pass

    def open_dump_folder(self):
        folder = filedialog.askdirectory(title='Select data directory (contains chain/ subfolder)')
        if not folder:
            return
        self._load_simulation_folder(folder)

    def _load_simulation_folder(self, folder: str, force_reload: bool = False) -> None:
        """Load a simulation folder."""
        self.playback_ctrl.pause()
        if self.data_ctrl.current_sim_folder is not None:
            if os.path.normpath(self.data_ctrl.current_sim_folder) != os.path.normpath(folder):
                self.clear_vtk_meshes()
        
        prev_ts = self.current_timestep
        ok, err = self.data_ctrl.load_folder(folder, force_reload=force_reload)
        if not ok:
            messagebox.showerror('Error', f'Failed to load dumps: {err}')
            return

        try: self.refresh_button.config(state=tk.NORMAL)
        except: pass

        if prev_ts is not None and prev_ts in self.data_ctrl.timesteps:
            try:
                idx = self.data_ctrl.timesteps.index(prev_ts)
                self.frame_slider.set(idx)
                self.ts_listbox.selection_clear(0, tk.END)
                self.ts_listbox.selection_set(idx)
                self.ts_listbox.activate(idx)
                self.show_timestep(prev_ts)
            except Exception: pass

    def refresh_current_source(self) -> None:
        """Reload the last opened/dropped simulation folder from disk."""
        if not self.data_ctrl.current_sim_folder:
            messagebox.showinfo('Refresh', 'No simulation folder loaded yet.')
            return
        self._load_simulation_folder(self.data_ctrl.current_sim_folder, force_reload=True)


    def open_dump_files(self):
        files = filedialog.askopenfilenames(title='Select dump files', filetypes=[('Dump files','*.dump'),('All','*.*')])
        if not files: return
        self.playback_ctrl.pause()
        self.clear_vtk_meshes()
        ok, err = self.data_ctrl.load_dump_files(files)
        if not ok: messagebox.showerror('Error', err)

    def on_drop(self, event):
        # event.data may be a list-like string; use tk splitlist for safety
        try:
            files = list(self.tk.splitlist(event.data))
        except Exception:
            # Fallback parsing
            raw = event.data.strip()
            files = [p.strip('{}') for p in raw.split()]

        # Normalize paths
        files = [self.data_ctrl._normalize_dropped_path(f) for f in files]
        self.handle_dropped_files(files)

    def handle_dropped_files(self, files):
        if not files:
            return

        # Enforce single dropped item (folder or file)
        if len(files) != 1:
            messagebox.showerror('Drop one item', 'Please drop a single simulation folder (or chain/ folder) or a single file.')
            return

        p = files[0]

        if os.path.isdir(p):
            sim_root = self.data_ctrl._resolve_sim_root(p)
            if sim_root is None:
                vtk_files = []
                for root, _, fs in os.walk(p):
                    for f in fs:
                        if f.lower().endswith('.vtk'):
                            vtk_files.append(os.path.join(root, f))
                if vtk_files:
                    for f in vtk_files: self._add_vtk_mesh(f)
                    if self.current_timestep is not None: self.show_timestep(self.current_timestep)
                    return
                messagebox.showerror('Invalid folder','Dropped folder must contain a chain/ subfolder or .vtk files.')
                return
            self._load_simulation_folder(sim_root)
            return

        if p.lower().endswith('.data'):
            ok, err = self.data_ctrl.load_data_file(p)
            if not ok: messagebox.showerror('Error', err)
            return
            return

        vtk_files = [f for f in files if f.lower().endswith('.vtk')]
        if vtk_files:
            for f in vtk_files: self._add_vtk_mesh(f)
            files = [f for f in files if not f.lower().endswith('.vtk')]
            if len(files) == 0:
                if self.current_timestep is not None: self.show_timestep(self.current_timestep)
                return
        else:
            self.clear_vtk_meshes()

        if files:
            ok, err = self.data_ctrl.load_dump_files(files)
            if not ok: messagebox.showerror('Error', err)

    def open_data_file(self):
        path = filedialog.askopenfilename(title='Select LAMMPS data file', filetypes=[('Data files','*.data'),('All','*.*')])
        if not path: return
        self.playback_ctrl.pause()
        self.clear_vtk_meshes()
        ok, err = self.data_ctrl.load_data_file(path)
        if not ok: messagebox.showerror('Error', err)


    def clear_vtk_meshes(self):
        self.vtk_ctrl.clear()
        if hasattr(self, 'vtk_listbox'):
            self.vtk_listbox.delete(0, tk.END)

    def open_vtk_files(self):
        files = filedialog.askopenfilenames(title='Select VTK files', filetypes=[('VTK files', '*.vtk'), ('All', '*.*')])
        if not files: return
        for f in files: self._add_vtk_mesh(f)
        if self.current_timestep is not None:
            self.show_timestep(self.current_timestep)
        self.plotter.render()

    def _update_persistent_bounds(self):
        """Update the max extents of the axes to include particles and all loaded VTKs, then keep them fixed."""
        self._persistent_view_bounds = self.vtk_ctrl.recompute_bounds(self.data_ctrl._init_limits, only_visible=False)

    def _add_vtk_mesh(self, path):
        ok, res = self.vtk_ctrl.add_mesh(path)
        if ok and hasattr(self, 'vtk_listbox'):
            self.vtk_listbox.insert(tk.END, res)
            self.vtk_listbox.selection_set(self.vtk_listbox.size() - 1)
            self._update_persistent_bounds()
        elif not ok: print(f"Failed to load VTK {path}: {res}")

    def on_vtk_select(self, event):
        if not hasattr(self, 'vtk_listbox'): return
        sel = self.vtk_listbox.curselection()
        for i in range(self.vtk_listbox.size()):
            self.vtk_ctrl.set_visibility(self.vtk_listbox.get(i), i in sel)
        if self.current_timestep is not None: self.show_timestep(self.current_timestep)
        else: self.plotter.render()

    def delete_selected_vtks(self):
        if not hasattr(self, 'vtk_listbox'): return
        sel = self.vtk_listbox.curselection()
        if not sel: return
        names = [self.vtk_listbox.get(i) for i in sel]
        self.vtk_ctrl.remove_meshes(names)
        for i in reversed(sel): self.vtk_listbox.delete(i)
        if self.current_timestep is not None: self.show_timestep(self.current_timestep)
        else: self.plotter.render()


    def on_ts_select(self, event):
        sel = self.ts_listbox.curselection()
        if not sel: return
        idx = sel[0]
        t = self.data_ctrl.timesteps[idx]
        self.frame_slider.set(idx)
        self.show_timestep(t)

    def on_slider(self, val):
        idx = int(float(val))
        if idx < 0 or idx >= len(self.data_ctrl.timesteps): return
        t = self.data_ctrl.timesteps[idx]
        self.ts_listbox.selection_clear(0, tk.END)
        self.ts_listbox.selection_set(idx)
        self.ts_listbox.activate(idx)
        self.playback_ctrl.update_status(idx)
        self.show_timestep(t)

    def fit_view(self):
        if self.plotter:
            self.plotter.reset_camera()
            self.plotter.render()

    def reset_view(self):
        if self.plotter:
            self.plotter.reset_camera()
            self.plotter.render()
        
    def _draw_chains_on_subset(self, subset: pd.DataFrame):
        """Draw lines connecting consecutive particles.

        Grouping strategy:
        - If a column named 'mol' or 'molecule' or 'chain' exists, group by that.
        - Otherwise, treat the whole subset as a single chain sorted by 'id'.
        """
        if subset.empty:
            return

        # Determine possible grouping column
        group_col = None
        for c in ('mol', 'molecule', 'chain', 'chain_id'):
            if c in subset.columns:
                group_col = c
                break

        if group_col is None:
            # single chain: sort by id
            pts = subset.sort_values('id')
            xs = pts['x'].values
            ys = pts['y'].values
            zs = pts['z'].values
            if xs.size > 1:
                self.ax.plot(xs, ys, zs, color='k', linewidth=0.8, alpha=0.8)
            return

        # multiple chains/groups
        for key, grp in subset.groupby(group_col):
            pts = grp.sort_values('id')
    def show_timestep(self, timestep):
        if self.data_ctrl.df is None: return
        subset = self.data_ctrl.df[self.data_ctrl.df['timestep'] == timestep]
        self.current_timestep = timestep

        for act in self._dynamic_actors: self.plotter.remove(act)
        self._dynamic_actors = []

        actors = self.plot_chain_data(subset)
        
        # Use persistent bounds so the axes stay fixed regardless of particle motion or visibility toggles
        if self._persistent_view_bounds is None:
            self._update_persistent_bounds()
        
        b = self._persistent_view_bounds
        axes = Axes(xrange=(b[0], b[1]), yrange=(b[2], b[3]), zrange=(b[4], b[5]),
                    xtitle='X', ytitle='Y', ztitle='Z', c='black')

        self._dynamic_actors.extend(actors)
        self._dynamic_actors.append(axes)

        if self.show_geometry_var.get():
            self._dynamic_actors.extend(self._build_geometry_actors(b))

        self.plotter.add(*self._dynamic_actors)
        self.vtk_ctrl.sync_visibility()
        self.hl_ctrl.render(self.data_ctrl.df, self.current_timestep)
        self.plotter.render()


    def _build_geometry_actors(self, bounds):
        if not self.data_ctrl.geometry_data: return []
        regions = self.data_ctrl.geometry_data.get('regions', {})
        box_id = self.data_ctrl.geometry_data.get('box_region')
        actors = []
        colors = ['red', 'green', 'blue', 'gold', 'cyan', 'magenta']
        color_idx = 0
        max_span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        plane_size = max_span if max_span > 0 else 1.0

        for r_id, r_data in regions.items():
            style = r_data.get('style')
            params = r_data.get('params', {})
            if style == 'union':
                continue

            is_box = (r_id == box_id)
            color = 'black' if is_box else colors[color_idx % len(colors)]
            if not is_box:
                color_idx += 1

            try:
                if style == 'block':
                    xlo, xhi = params['xlo'], params['xhi']
                    ylo, yhi = params['ylo'], params['yhi']
                    zlo, zhi = params['zlo'], params['zhi']
                    center = ((xlo + xhi) / 2.0, (ylo + yhi) / 2.0, (zlo + zhi) / 2.0)
                    dims = (xhi - xlo, yhi - ylo, zhi - zlo)
                    actor = Box(pos=center, length=dims[0], width=dims[1], height=dims[2]).c(color).alpha(0.10)
                    if is_box:
                        actor = actor.wireframe().lw(1).alpha(1.0)
                    actors.append(actor)
                elif style == 'cylinder':
                    dim = params['dim']
                    c1, c2 = params['c1'], params['c2']
                    radius = params['radius']
                    lo, hi = params['lo'], params['hi']
                    if dim == 'x':
                        pos = ((lo + hi) / 2.0, c1, c2)
                        axis = (1, 0, 0)
                    elif dim == 'y':
                        pos = (c1, (lo + hi) / 2.0, c2)
                        axis = (0, 1, 0)
                    else:
                        pos = (c1, c2, (lo + hi) / 2.0)
                        axis = (0, 0, 1)
                    actors.append(Cylinder(pos=pos, r=radius, height=abs(hi - lo), axis=axis).c(color).alpha(0.10))
                elif style == 'cone':
                    dim = params['dim']
                    c1, c2 = params['c1'], params['c2']
                    radlo, radhi = params['radlo'], params['radhi']
                    lo, hi = params['lo'], params['hi']
                    avg_r = max((radlo + radhi) / 2.0, 1e-6)
                    if dim == 'x':
                        pos = ((lo + hi) / 2.0, c1, c2)
                        axis = (1, 0, 0)
                    elif dim == 'y':
                        pos = (c1, (lo + hi) / 2.0, c2)
                        axis = (0, 1, 0)
                    else:
                        pos = (c1, c2, (lo + hi) / 2.0)
                        axis = (0, 0, 1)
                    actors.append(Cone(pos=pos, axis=axis, r=avg_r, height=abs(hi - lo)).c(color).alpha(0.10))
                elif style == 'plane':
                    px, py, pz = params['px'], params['py'], params['pz']
                    nx, ny, nz = params['nx'], params['ny'], params['nz']
                    actors.append(Plane(pos=(px, py, pz), normal=(nx, ny, nz), s=(plane_size, plane_size)).c(color).alpha(0.10))
            except Exception:
                continue

        return actors

    def plot_chain_data(self, chain_data):
        if chain_data.empty:
            return []

        pts = chain_data[['x','y','z']].values
    
        r = 0.005
        if 'diameter' in chain_data.columns:
            try:
                r = (chain_data['diameter'].values / 2.0)
            except Exception:
                pass

        # Create spheres (GPU instanced → fast)
        spheres = Spheres(pts, r=r, c='blue', alpha=0.6)

        actors = [spheres]

        # Optional: draw chains
        if getattr(self, 'draw_chains_var', None) and self.draw_chains_var.get():
            if len(pts) > 1:
                lines = Lines(pts, c='black', lw=1)
                actors.append(lines)

        return actors
    


    def _apply_highlight(self):
        mode = self.highlight_mode_var.get()
        id_str = self.highlight_id_var.get().strip()
        if not id_str: return
        try: n = int(id_str)
        except ValueError:
            self.highlight_status.config(text='Enter an integer ID.', fg='red')
            return
        ok, res = self.hl_ctrl.apply(mode, n)
        self.highlight_status.config(text=res, fg='black' if ok else 'red')
        self.plotter.render()

    def _clear_highlight(self):
        self.hl_ctrl.clear()
        self.highlight_status.config(text='—', fg='gray')
        self.plotter.render()


    # ── Time-series window openers ────────────────────────────────

    def _open_bond_win(self):
        if self.data_ctrl.df_mi is None:
            from tkinter import messagebox as _mb
            _mb.showinfo('No data', 'Load a simulation first.')
            return
        if self._bond_win is None or not self._bond_win.winfo_exists():
            self._bond_win = BondPlotWindow(
                self, self.data_ctrl.df_bonds, self.chain_size_var, self._bond_range_spec)
        else:
            self._bond_win.lift()

    def _open_angle_win(self):
        if self.data_ctrl.df_mi is None:
            from tkinter import messagebox as _mb
            _mb.showinfo('No data', 'Load a simulation first.')
            return
        if self._angle_win is None or not self._angle_win.winfo_exists():
            self._angle_win = AnglePlotWindow(
                self, self.data_ctrl.df_angles, self.chain_size_var, self._angle_range_spec)
        else:
            self._angle_win.lift()

    def _open_atom_win(self):
        if self.data_ctrl.df_mi is None:
            from tkinter import messagebox as _mb
            _mb.showinfo('No data', 'Load a simulation first.')
            return
        if self._atom_win is None or not self._atom_win.winfo_exists():
            self._atom_win = AtomPlotWindow(
                self, self.data_ctrl.df_mi, self.chain_size_var, self._atom_range_spec)
        else:
            self._atom_win.lift()

    def _cleanup(self):
        """Cleanup resources before exiting."""
        if hasattr(self, 'playback_ctrl'):
            self.playback_ctrl.pause()
        try:
            import matplotlib.pyplot as _plt
            _plt.close('all')
        except Exception:
            pass

    def _on_close(self):
        # Called when window is closed via window manager
        self._cleanup()
        try:
            self.destroy()
        except Exception:
            pass
        # Force exit to ensure background threads (if any) don't block terminal
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
