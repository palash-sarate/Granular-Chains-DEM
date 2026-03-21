import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
            df = pd.DataFrame({'id': ids, 'type': types, 'x': x, 'y': y, 'z': z})
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
                    rows.append((idv, x, y, z))
            except Exception:
                continue
        if rows:
            df = pd.DataFrame(rows, columns=['id','x','y','z'])

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
        self.geometry('1000x700')
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
        tk.Button(ctrl, text='Close', command=self.destroy).pack(fill=tk.X, pady=(20,0))

        # Right: Matplotlib canvas
        self.fig = plt.figure(figsize=(7,6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)

        # Data holders
        self.df = None
        self.timesteps = []
        self.current_timestep = None

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

        self.timesteps = sorted(self.df['timestep'].unique())
        self.ts_listbox.delete(0, tk.END)
        for t in self.timesteps:
            self.ts_listbox.insert(tk.END, str(t))

        # setup slider
        if self.timesteps:
            self.frame_slider.config(from_=0, to=max(0, len(self.timesteps)-1))
            self.frame_slider.set(0)
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
        # recompute limits based on current dataframe
        if self.df is None:
            return
        x_min, x_max = self.df['x'].min(), self.df['x'].max()
        y_min, y_max = self.df['y'].min(), self.df['y'].max()
        z_min, z_max = self.df['z'].min(), self.df['z'].max()
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_zlim(z_min, z_max)
        self.canvas.draw()

    def show_timestep(self, timestep):
        if self.df is None:
            return
        subset = self.df[self.df['timestep'] == timestep]
        if subset.empty:
            return

        self.ax.cla()
        x = subset['x'].values
        y = subset['y'].values
        z = subset['z'].values

        # size
        if 'diameter' in subset.columns:
            sizes = (subset['diameter'].values / subset['diameter'].max()) * 100
        else:
            sizes = np.full_like(x, 30)

        self.ax.scatter(x, y, z, s=sizes, c='C0', depthshade=True)

        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.set_title(f'Timestep: {timestep}  (particles: {len(subset)})')

        # set equal aspect (best-effort)
        try:
            self.ax.set_box_aspect((1,1,1))
        except Exception:
            pass

        self.canvas.draw()

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
