import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
from vedo import Plotter, Spheres, Lines, Axes
from utilities import validate_chain_spacing

from data_manager import SimulationData

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

def parse_simple_data_file(path: str) -> pd.DataFrame:
    """Attempt a best-effort parse of a LAMMPS data file to extract atom positions.
    Returns a DataFrame with columns ['id','type','x','y','z','vx','vy','vz','diameter'] when available.
    This is intentionally permissive and will return an empty DataFrame if parsing fails.
    """
    if not os.path.exists(path):
        return pd.DataFrame()

    with open(path, 'r') as f:
        lines = f.readlines()

    # Find the 'Atoms' section
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith('Atoms'):
            start = i + 1
            break

    if start is None:
        return pd.DataFrame()

    # Advance past any blank/comment lines after the 'Atoms' header
    idx = start
    while idx < len(lines) and (lines[idx].strip() == '' or lines[idx].strip().startswith('#')):
        idx += 1

    # Read numeric lines until next blank or a new section header (e.g., Velocities, Bonds, Angles, Bonds, Masses)
    SECTION_HEADERS = set(['Velocities', 'Bonds', 'Angles', 'Masses', 'PairIJ', 'Velocities', 'Bonds', 'Angles'])
    data_lines = []
    for line in lines[idx:]:
        s = line.strip()
        if s == '':
            break
        # Stop if line looks like a section header (word with no numbers)
        first_tok = s.split()[0]
        if first_tok in SECTION_HEADERS:
            break
        # skip comments
        if s.startswith('#'):
            continue
        parts = line.split()
        # require at least id and x y z (3 coords + id -> 4)
        if len(parts) < 4:
            # if a short non-data line appears, stop parsing atoms
            break
        data_lines.append(parts)

    if not data_lines:
        return pd.DataFrame()

    # Convert to DataFrame trying to map common formats.
    # Typical atom styles: id mol type x y z ... or id type x y z ...
    # We'll try to detect whether second token is integer (mol) or float (x)
    first = data_lines[0]
    df = None
    try:
        # try id type x y z
        arr = np.array(data_lines, dtype=float)
        # If successful, map columns
        if arr.shape[1] >= 5:
            # id,type,x,y,z
            ids = arr[:,0].astype(int)
            types = arr[:,1].astype(int)
            x = arr[:,2]
            y = arr[:,3]
            z = arr[:,4]
            dia = arr[:,5]
            df = pd.DataFrame({'id': ids, 'type': types, 'x': x, 'y': y, 'z': z, 'diameter': dia})
    except Exception:
        # fallback: attempt to parse as mixed tokens
        rows = []
        for parts in data_lines:
            try:
                # assume first token id, last three are x y z
                if len(parts) >= 4:
                    idv = int(parts[0])
                    x = float(parts[-3])
                    y = float(parts[-2])
                    z = float(parts[-1])
                    dia = float(parts[4]) if len(parts) > 5 else 0.0
                    rows.append((idv, x, y, z, dia))
            except Exception:
                continue
        if rows:
            df = pd.DataFrame(rows, columns=['id','x','y','z','diameter'])

    if df is None:
        return pd.DataFrame()

    # Ensure columns exist
    for c in ['vx','vy','vz','diameter']:
        if c not in df.columns:
            df[c] = 0.0

    # attach timestep 0 for consistency with Animator expectations
    df['timestep'] = 0
    df.set_index(['timestep','id'], inplace=True)
    return df


class ViewerApp(BaseTk):
    def __init__(self):
        super().__init__()
        self.title('Chains Simulation Viewer')
        self.geometry('300x700')
        # Ensure clean shutdown when window is closed
        self.protocol("WM_DELETE_WINDOW", self._on_close)

        # Left controls
        ctrl = tk.Frame(self)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

        tk.Button(ctrl, text='Open Dump Folder', command=self.open_dump_folder).pack(fill=tk.X)
        tk.Button(ctrl, text='Open Dump Files...', command=self.open_dump_files).pack(fill=tk.X, pady=(4,0))
        tk.Button(ctrl, text='Open Data File...', command=self.open_data_file).pack(fill=tk.X, pady=(4,0))

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
        tk.Button(ctrl, text='Reset View', command=self.reset_view).pack(fill=tk.X, pady=(4,0))
        tk.Button(ctrl, text='Close', command=self.destroy).pack(fill=tk.X, pady=(20,0))

        self.plotter = Plotter(
            bg='white',
            interactive=True
        )

        # Data holders
        self.df = None
        self.timesteps = []
        self.current_timestep = None
        # Store initial axis limits and default view angles
        self._init_limits = None
        self._default_view = (30, -60)
        try:
            self.ax.view_init(elev=self._default_view[0], azim=self._default_view[1])
        except Exception:
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
        except Exception as e:
            messagebox.showerror('Error', f'Failed to load dumps: {e}')

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

    def on_drop(self, event):
        # event.data may be a list-like string; use tk splitlist for safety
        try:
            files = list(self.tk.splitlist(event.data))
        except Exception:
            # Fallback parsing
            raw = event.data.strip()
            files = [p.strip('{}') for p in raw.split()]

        # Normalize paths
        files = [f.replace('\\', '/') for f in files]
        self.handle_dropped_files(files)

    def handle_dropped_files(self, files):
        if not files:
            return
        # If single .data file -> open as data file
        if len(files) == 1 and files[0].lower().endswith('.data'):
            try:
                df = parse_simple_data_file(files[0])
                if df.empty:
                    messagebox.showerror('Parse failed', 'Could not parse the dropped data file')
                    return
                df = df.reset_index()
                self.load_dataframe(df)
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
        self.plotter.add(axes)

        # ---- Force camera to respect bounds ----
        self.plotter.reset_camera()

        self.plotter.render()

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
