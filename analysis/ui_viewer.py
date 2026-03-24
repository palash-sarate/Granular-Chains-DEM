import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional
import pandas as pd
import numpy as np
from vedo import Plotter, Spheres, Lines, Axes, Box, Cylinder, Cone, Plane
from utilities import validate_chain_spacing
from data_manager import SimulationData, parse_simple_data_file, load_lammps_geometry

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
        # Let the window size to its content, then center.
        # (Avoid a fixed width that leaves empty space to the right.)
        # Ensure clean shutdown when window is closed
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        # Left controls
        ctrl = tk.Frame(self)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        tk.Button(ctrl, text='Open Dump Folder', command=self.open_dump_folder).pack(fill=tk.X)
        tk.Button(ctrl, text='Open Dump Files...', command=self.open_dump_files).pack(fill=tk.X, pady=(4,0))
        tk.Button(ctrl, text='Open Data File...', command=self.open_data_file).pack(fill=tk.X, pady=(4,0))
        self.refresh_button = tk.Button(ctrl, text='Refresh', command=self.refresh_current_source, state=tk.DISABLED)
        self.refresh_button.pack(fill=tk.X, pady=(8,0))

        # Drag-and-drop area
        drop_text = 'Drop dump/.data files here' if DND_AVAILABLE else 'Drag-and-drop disabled (install tkinterdnd2)'
        self.drop_label = tk.Label(ctrl, text=drop_text, relief='ridge', width=30, height=4)
        self.drop_label.pack(fill=tk.X, pady=(8,4))
        if DND_AVAILABLE:
            try:
                self.drop_label.drop_target_register(DND_CONST)
                self.drop_label.dnd_bind('<<Drop>>', self.on_drop)
            except Exception:
                # if registration fails, keep graceful fallback
                self.drop_label.config(text='Drag-and-drop not available on this platform')

        tk.Label(ctrl, text='Timesteps:').pack(anchor='w', pady=(8,0))
        self.ts_listbox = tk.Listbox(ctrl, width=30, height=20)
        self.ts_listbox.pack(fill=tk.Y)
        self.ts_listbox.bind('<<ListboxSelect>>', self.on_ts_select)

        tk.Label(ctrl, text='Frame').pack(anchor='w', pady=(8,0))
        self.frame_slider = tk.Scale(ctrl, from_=0, to=0, orient=tk.HORIZONTAL, command=self.on_slider)
        self.frame_slider.pack(fill=tk.X)

        tk.Button(ctrl, text='Fit View', command=self.fit_view).pack(fill=tk.X, pady=(8,0))
        # Option to draw chain lines between consecutive particles
        self.draw_chains_var = tk.BooleanVar(value=False)
        tk.Checkbutton(ctrl, text='Draw Chains', variable=self.draw_chains_var, command=self._on_draw_toggle).pack(fill=tk.X, pady=(4,0))
        self.show_geometry_var = tk.BooleanVar(value=False)
        tk.Checkbutton(ctrl, text='Show Geometry', variable=self.show_geometry_var, command=self._on_draw_toggle).pack(fill=tk.X, pady=(4,0))
        tk.Button(ctrl, text='Reset View', command=self.reset_view).pack(fill=tk.X, pady=(4,0))
        tk.Button(ctrl, text='Close', command=self.destroy).pack(fill=tk.X, pady=(20,0))

        self.plotter = Plotter(
            bg='white',
            interactive=True
        )

        # Size-to-content and center on screen after widgets are laid out
        self._autosize_and_center()

        # Data holders
        self.df = None
        self.timesteps = []
        self.current_timestep = None
        self.current_sim_folder = None
        self.geometry_data = None
        self.geometry_script_path = None
        # Store initial axis limits and default view angles
        self._init_limits = None
        self._default_view = (30, -60)
        try:
            self.ax.view_init(elev=self._default_view[0], azim=self._default_view[1])
        except Exception:
            pass

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

    def _load_simulation_folder(self, folder: str) -> None:
        """Load a simulation folder using SimulationData (expects folder/chain/*.dump)."""
        prev_ts = self.current_timestep
        try:
            sim = SimulationData(folder)
            df = sim.load_data(force_reload=True)
            if df.empty:
                messagebox.showerror('No data', 'No dump files found or parsing failed in selected folder')
                return
            # convert MultiIndex to columns for easier handling
            if isinstance(df.index, pd.MultiIndex):
                df = df.reset_index()
            self.load_dataframe(df)
            self.current_sim_folder = folder
            self.geometry_data, self.geometry_script_path = load_lammps_geometry(folder)
            try:
                self.refresh_button.config(state=tk.NORMAL)
            except Exception:
                pass

            # Try to keep the same timestep selected after refresh/reload
            if prev_ts is not None and prev_ts in self.timesteps:
                try:
                    idx = self.timesteps.index(prev_ts)
                    self.frame_slider.set(idx)
                    self.ts_listbox.selection_clear(0, tk.END)
                    self.ts_listbox.selection_set(idx)
                    self.ts_listbox.activate(idx)
                    self.show_timestep(prev_ts)
                except Exception:
                    pass
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load dumps: {e}')

    def refresh_current_source(self) -> None:
        """Reload the last opened/dropped simulation folder from disk."""
        if not self.current_sim_folder:
            messagebox.showinfo('Refresh', 'No simulation folder loaded yet. Use Open Dump Folder or drop a folder first.')
            return
        self._load_simulation_folder(self.current_sim_folder)

    def _normalize_dropped_path(self, p: str) -> str:
        p = (p or "").strip()
        if p.startswith("{") and p.endswith("}"):
            p = p[1:-1]
        # Keep OS-native separators for filesystem checks
        try:
            return os.path.normpath(p)
        except Exception:
            return p

    def _resolve_sim_root_from_drop_dir(self, dropped_dir: str) -> Optional[str]:
        """Given a dropped directory, resolve the simulation root directory.

        Accepts either:
        - the simulation root (contains chain/)
        - the chain/ directory itself (parent is the simulation root)
        """
        d = self._normalize_dropped_path(dropped_dir)
        if not d or not os.path.isdir(d):
            return None

        chain_dir = os.path.join(d, "chain")
        if os.path.isdir(chain_dir):
            return d

        base = os.path.basename(d).lower()
        if base == "chain":
            parent = os.path.dirname(d)
            if parent and os.path.isdir(os.path.join(parent, "chain")):
                return parent

        return None

    def open_dump_files(self):
        files = filedialog.askopenfilenames(title='Select dump files', filetypes=[('Dump files','*.dump'),('All','*.*')])
        if not files:
            return
        # Create a temporary in-memory concatenated DataFrame
        frames = []
        for f in files:
            try:
                # Use SimulationData parsing helper by emulating a chain folder
                # Here we reuse the private parser by instantiating SimulationData pointing to a temp dir is complex,
                # so parse single dump by reading ITEM: TIMESTEP and ATOMS sections
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
                        frames.append(df)
                        break
            except Exception as e:
                print(f'Failed to parse {f}: {e}')

        if not frames:
            messagebox.showerror('Error', 'No valid dump frames parsed')
            return
        full = pd.concat(frames)
        full.set_index(['timestep','id'], inplace=True)
        full = full.reset_index()
        self.load_dataframe(full)
        self.geometry_data = None
        self.geometry_script_path = None

    def on_drop(self, event):
        # event.data may be a list-like string; use tk splitlist for safety
        try:
            files = list(self.tk.splitlist(event.data))
        except Exception:
            # Fallback parsing
            raw = event.data.strip()
            files = [p.strip('{}') for p in raw.split()]

        # Normalize paths
        files = [self._normalize_dropped_path(f) for f in files]
        self.handle_dropped_files(files)

    def handle_dropped_files(self, files):
        if not files:
            return

        # Enforce single dropped item (folder or file)
        if len(files) != 1:
            messagebox.showerror('Drop one item', 'Please drop a single simulation folder (or chain/ folder) or a single file.')
            return

        p = files[0]

        # Dropped directory: treat as simulation folder (root or chain/)
        if os.path.isdir(p):
            sim_root = self._resolve_sim_root_from_drop_dir(p)
            if sim_root is None:
                messagebox.showerror(
                    'Invalid folder',
                    'Dropped folder must be the simulation run folder containing a chain/ subfolder, '
                    'or the chain/ subfolder itself.'
                )
                return
            self._load_simulation_folder(sim_root)
            return

        # If single .data file -> open as data file
        if p.lower().endswith('.data'):
            try:
                df = parse_simple_data_file(p)
                if df.empty:
                    messagebox.showerror('Parse failed', 'Could not parse the dropped data file')
                    return
                df = df.reset_index()
                self.load_dataframe(df)
                self.geometry_data = None
                self.geometry_script_path = None
            except Exception as e:
                messagebox.showerror('Error', f'Failed to load data file: {e}')
            return

        # Otherwise treat as dump files (one or many)
        frames = []
        for f in files:
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
                        frames.append(df)
                        break
            except Exception as e:
                print(f'Failed to parse dropped file {f}: {e}')

        if not frames:
            messagebox.showerror('Error', 'No valid dump frames parsed from dropped files')
            return
        full = pd.concat(frames)
        full.set_index(['timestep','id'], inplace=True)
        full = full.reset_index()
        self.load_dataframe(full)
        self.geometry_data = None
        self.geometry_script_path = None

    def open_data_file(self):
        path = filedialog.askopenfilename(title='Select LAMMPS data file', filetypes=[('Data files','*.data'),('All','*.*')])
        if not path:
            return
        df = parse_simple_data_file(path)
        if df.empty:
            messagebox.showerror('Parse failed', 'Could not parse the selected data file')
            return
        # reset index to columns
        df = df.reset_index()
        self.load_dataframe(df)
        self.geometry_data = None
        self.geometry_script_path = None

    def load_dataframe(self, df: pd.DataFrame):
        self.df = df.copy()
        # ensure columns present
        for c in ['x','y','z','id','timestep']:
            if c not in self.df.columns:
                messagebox.showerror('Invalid data', f'Missing column: {c}')
                return
            
        # compute and store initial/global axis limits for reset
        x_min, x_max = self.df['x'].min(), self.df['x'].max()
        y_min, y_max = self.df['y'].min(), self.df['y'].max()
        z_min, z_max = self.df['z'].min(), self.df['z'].max()
        self._init_limits = ((x_min, x_max), (y_min, y_max), (z_min, z_max))

        self.timesteps = sorted(self.df['timestep'].unique())
        self.ts_listbox.delete(0, tk.END)
        for t in self.timesteps:
            self.ts_listbox.insert(tk.END, str(t))

        self.frame_slider.config(from_=0, to=max(0, len(self.timesteps)-1))
        self.frame_slider.set(0)
        if self.timesteps:
            self.show_timestep(self.timesteps[0])

    def on_ts_select(self, event):
        sel = self.ts_listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        t = self.timesteps[idx]
        self.frame_slider.set(idx)
        self.show_timestep(t)

    def on_slider(self, val):
        idx = int(float(val))
        if idx < 0 or idx >= len(self.timesteps):
            return
        t = self.timesteps[idx]
        # update listbox selection without triggering event
        self.ts_listbox.selection_clear(0, tk.END)
        self.ts_listbox.selection_set(idx)
        self.ts_listbox.activate(idx)
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
            xs = pts['x'].values
            ys = pts['y'].values
            zs = pts['z'].values
            if xs.size > 1:
                self.ax.plot(xs, ys, zs, linewidth=0.8, alpha=0.8)

    def show_timestep(self, timestep):
        if self.df is None:
            return

        subset = self.df[self.df['timestep'] == timestep]
        self.current_timestep = timestep

        self.plotter.clear()

        actors = self.plot_chain_data(subset)

        # ---- Compute equal bounding box ----
        xs = subset['x'].values
        ys = subset['y'].values
        zs = subset['z'].values

        xmin, xmax = xs.min(), xs.max()
        ymin, ymax = ys.min(), ys.max()
        zmin, zmax = zs.min(), zs.max()

        # center
        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        cz = 0.5 * (zmin + zmax)

        # max range → enforce cube
        max_range = max(xmax - xmin, ymax - ymin, zmax - zmin) / 2.0

        bounds = [
            cx - max_range, cx + max_range,
            cy - max_range, cy + max_range,
            cz - max_range, cz + max_range
        ]

        axes = Axes(
            xrange=(bounds[0], bounds[1]),
            yrange=(bounds[2], bounds[3]),
            zrange=(bounds[4], bounds[5]),
            xtitle='X',
            ytitle='Y',
            ztitle='Z',
            c='black'
        )

        # ---- Add everything ----
        self.plotter.add(*actors)
        if getattr(self, 'show_geometry_var', None) and self.show_geometry_var.get():
            self.plotter.add(*self._build_geometry_actors(bounds))
        self.plotter.add(axes)

        # ---- Force camera to respect bounds ----
        self.plotter.reset_camera()

        self.plotter.render()

    def _build_geometry_actors(self, bounds):
        if not self.geometry_data:
            return []

        regions = self.geometry_data.get('regions', {})
        box_id = self.geometry_data.get('box_region')
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
    
        r = (chain_data['diameter'].values / 2.0)

        # Create spheres (GPU instanced → fast)
        spheres = Spheres(pts, r=r, c='blue', alpha=0.6)

        actors = [spheres]

        # Optional: draw chains
        if getattr(self, 'draw_chains_var', None) and self.draw_chains_var.get():
            if len(pts) > 1:
                lines = Lines(pts, c='black', lw=1)
                actors.append(lines)

        return actors
    
    def _cleanup(self):
        """Cleanup resources before exiting."""
        try:
            plt.close('all')
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
