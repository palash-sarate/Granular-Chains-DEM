import os
import sys
import threading
from queue import Queue
from typing import Optional, List
import pandas as pd
from analysis.data_manager import SimulationData, load_lammps_geometry, parse_simple_data_file

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
        self.timesteps: List[int] = []
        self._init_limits: Optional[tuple] = None
        self.sim_source: Optional[SimulationData] = None
        
        self.on_data_loaded_cb = None
        self.on_batch_ready_cb = None
        
        # Async loading state
        self.load_queue = Queue()
        self.loading_batches = set()
        self._worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._worker_thread.start()

    def load_folder(self, folder: str, force_reload: bool = False):
        try:
            self.sim_source = SimulationData(folder)
            data_dict = self.sim_source.load_data(force_reload=force_reload)
            df_atoms = data_dict['atoms']
            
            if df_atoms is None or df_atoms.empty:
                return False, "No dump files found or parsing failed."
            
            self.current_sim_folder = folder
            self.geometry_data, self.geometry_script_path = load_lammps_geometry(folder)
            
            # Metadata: Store ALL timesteps even if data isn't loaded yet
            self.timesteps = self.sim_source.timesteps
            
            self.df_bonds = data_dict['bonds']
            self.df_angles = data_dict['angles']
            
            self.load_dataframe(df_atoms, update_timesteps=False)
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

    def load_dataframe(self, df: pd.DataFrame, update_timesteps=True):
        if df is None: return
        self.df = df.copy()
        if 'timestep' in self.df.index.names:
            self.df = self.df.reset_index()

        # Ensure critical columns are numeric
        for c in ['x', 'y', 'z', 'id', 'timestep']:
            if c in self.df.columns:
                self.df[c] = pd.to_numeric(self.df[c], errors='coerce')
        
        if 'diameter' not in self.df.columns:
            self.df['diameter'] = 0.01
        else:
            self.df['diameter'] = pd.to_numeric(self.df['diameter'], errors='coerce').fillna(0.01)

        self.df.dropna(subset=['x', 'y', 'z', 'id', 'timestep'], inplace=True)

        if not self.df.empty:
            x_min, x_max = self.df['x'].min(), self.df['x'].max()
            y_min, y_max = self.df['y'].min(), self.df['y'].max()
            z_min, z_max = self.df['z'].min(), self.df['z'].max()
            self._init_limits = ((x_min, x_max), (y_min, y_max), (z_min, z_max))

        try:
            self.df_mi = self.df.set_index(['timestep', 'id']).sort_index()
        except Exception:
            self.df_mi = None

        if update_timesteps:
            self.timesteps = sorted(self.df['timestep'].unique())
            
        if self.on_data_loaded_cb:
            self.on_data_loaded_cb()

    def request_batch_for_timestep(self, ts):
        if not self.sim_source: return
        batch_idx = self.sim_source.get_batch_index_for_timestep(ts)
        if batch_idx not in self.sim_source.loaded_batches and batch_idx not in self.loading_batches:
            self.loading_batches.add(batch_idx)
            self.load_queue.put(batch_idx)

    def _worker_loop(self):
        while True:
            batch_idx = self.load_queue.get()
            if batch_idx is None: break
            try:
                print(f"Background loading batch {batch_idx}...")
                data_dict = self.sim_source.load_batch(batch_idx)
                
                # Merge new data into main dataframes
                if not data_dict['atoms'].empty:
                    new_atoms = data_dict['atoms'].reset_index()
                    self.df = pd.concat([self.df, new_atoms])
                    self.df = self.df[~self.df.duplicated(subset=['timestep', 'id'], keep='last')]
                    self.df_mi = self.df.set_index(['timestep', 'id']).sort_index()
                
                if not data_dict['bonds'].empty:
                    self.df_bonds = pd.concat([self.df_bonds, data_dict['bonds']])
                    self.df_bonds = self.df_bonds[~self.df_bonds.index.duplicated(keep='last')].sort_index()
                
                if not data_dict['angles'].empty:
                    self.df_angles = pd.concat([self.df_angles, data_dict['angles']])
                    self.df_angles = self.df_angles[~self.df_angles.index.duplicated(keep='last')].sort_index()

                if self.on_batch_ready_cb:
                    self.on_batch_ready_cb(batch_idx)
            except Exception as e:
                print(f"Error loading batch {batch_idx}: {e}")
            finally:
                self.loading_batches.remove(batch_idx)
                self.load_queue.task_done()

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

    def normalize_dropped_path(self, p: str) -> str:
        p = (p or "").strip()
        if p.startswith("{") and p.endswith("}"):
            p = p[1:-1]
        try:
            return os.path.normpath(p)
        except Exception:
            return p

    def resolve_sim_root(self, dropped_dir: str) -> Optional[str]:
        d = self.normalize_dropped_path(dropped_dir)
        if not d or not os.path.isdir(d):
            return None
        if os.path.isdir(os.path.join(d, "chain")):
            return d
        if os.path.basename(d).lower() == "chain":
            parent = os.path.dirname(d)
            if parent and os.path.isdir(os.path.join(parent, "chain")):
                return parent
        return None
