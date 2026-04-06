import os
import sys
import threading
from queue import Queue
from typing import Optional, List
import pandas as pd
from analysis.data_manager import SimulationData, load_lammps_geometry, parse_simple_data_file

class SimDataController:
    """Manages loading, parsing, and caching of simulation data."""
    @property
    def df(self) -> pd.DataFrame:
        """Lazy-fused monolithic atoms DataFrame."""
        if self.is_preview_mode and self.preview_df is not None:
            return self.preview_df.reset_index()
            
        if not self.sim_source:
            return self._df_manual if self._df_manual is not None else pd.DataFrame()
        
        # Merge all currently loaded batches for global analysis
        loaded_dfs = list(self.sim_source.atom_batches.values())
        if not loaded_dfs: return pd.DataFrame()
        # Ensure we return a reset index version (as expected by older code)
        return pd.concat(loaded_dfs).reset_index()

    @property
    def df_mi(self) -> Optional[pd.DataFrame]:
        """Lazy-fused MultiIndex atoms DataFrame."""
        if self.is_preview_mode and self.preview_df is not None:
            return self.preview_df
            
        if not self.sim_source:
            return self._df_mi_manual
        
        loaded_dfs = list(self.sim_source.atom_batches.values())
        if not loaded_dfs: return None
        return pd.concat(loaded_dfs)

    @property
    def df_bonds(self) -> pd.DataFrame:
        """Lazy-fused monolithic bonds DataFrame."""
        if self.is_preview_mode:
            # Restart files don't typically export bonds to the preview dump easily
            return pd.DataFrame()
            
        if not self.sim_source:
            return self._df_bonds_manual if self._df_bonds_manual is not None else pd.DataFrame()
            
        loaded_dfs = list(self.sim_source.bond_batches.values())
        if not loaded_dfs: return pd.DataFrame()
        return pd.concat(loaded_dfs)

    @property
    def df_angles(self) -> pd.DataFrame:
        """Lazy-fused monolithic angles DataFrame."""
        if self.is_preview_mode:
            return pd.DataFrame()
            
        if not self.sim_source:
            return self._df_angles_manual if self._df_angles_manual is not None else pd.DataFrame()
            
        loaded_dfs = list(self.sim_source.angle_batches.values())
        if not loaded_dfs: return pd.DataFrame()
        return pd.concat(loaded_dfs)

    def __init__(self):
        # Internal storage for manual/legacy loads
        self._df_manual: Optional[pd.DataFrame] = None
        self._df_mi_manual: Optional[pd.DataFrame] = None
        self._df_bonds_manual: Optional[pd.DataFrame] = None
        self._df_angles_manual: Optional[pd.DataFrame] = None
        self.current_sim_folder: Optional[str] = None
        self.geometry_data: Optional[dict] = None
        self.geometry_script_path: Optional[str] = None
        
        # Ghost/Preview Mode state for Restart Editor
        self.preview_df: Optional[pd.DataFrame] = None
        self.is_preview_mode: bool = False
        
        # New Simulation Mode state
        self.sim_source: Optional[SimulationData] = None
        
        self.timesteps: List[int] = []
        self._init_limits: Optional[tuple] = None
        
        self.on_data_loaded_cb = None
        self.on_batch_ready_cb = None
        
        # Async loading state
        self.load_queue = Queue()
        self.loading_batches = set()
        self._worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self._worker_thread.start()

    def load_folder(self, folder: str, force_reload: bool = False, enable_preloading: bool = True):
        try:
            new_source = SimulationData(folder)
            data_dict = new_source.load_data(force_reload=force_reload)
            df_atoms = data_dict['atoms']
            
            if df_atoms is None or df_atoms.empty:
                return False, "No dump files found or parsing failed.", 0
            
            # Commit changes only after success
            self.sim_source = new_source
            self.current_sim_folder = folder
            self.geometry_data, self.geometry_script_path = load_lammps_geometry(folder)
            
            # Metadata: Store ALL timesteps even if data isn't loaded yet
            self.timesteps = self.sim_source.timesteps
            
            # Identify cached batches for preloading if enabled
            cached_batches = []
            if enable_preloading:
                cached_batches = self.sim_source.get_available_cached_batches()
                
                # Batch 0 is already loaded by sim_source.load_data(), 
                # remove it from the background queue if it was in the cached list
                if 0 in cached_batches:
                    cached_batches.remove(0)
                
                # Queue remaining cached batches for background loading
                for b_idx in cached_batches:
                    if b_idx not in self.loading_batches:
                        self.loading_batches.add(b_idx)
                        self.load_queue.put(b_idx)
            
            self._df_bonds_manual = data_dict['bonds']
            self._df_angles_manual = data_dict['angles']
            
            self.load_dataframe(df_atoms, update_timesteps=False)
            
            # Return result and the number of batches we just queued
            return True, None, len(cached_batches)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return False, str(e), 0

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
        self.sim_source = None
        self.current_sim_folder = None
        self._dynamic_actors = []
        self._persistent_view_bounds = None
        self._axes = None
        self.load_dataframe(full)
        return True, None, 0

    def load_data_file(self, path: str):
        try:
            df = parse_simple_data_file(path)
            if df.empty:
                return False, "Could not parse the selected data file", 0
            df = df.reset_index()
            self.sim_source = None
            self.current_sim_folder = None
            self.geometry_data = None
            self.geometry_script_path = None
            self.load_dataframe(df)
            return True, None, 0
        except Exception as e:
            return False, str(e), 0

    def load_dataframe(self, df: pd.DataFrame, update_timesteps=True):
        if df is None: return
        self._df_manual = df.copy()
        if 'timestep' in self._df_manual.index.names:
            self._df_manual = self._df_manual.reset_index()

        # Ensure critical columns are numeric
        for c in ['x', 'y', 'z', 'id', 'timestep']:
            if c in self._df_manual.columns:
                self._df_manual[c] = pd.to_numeric(self._df_manual[c], errors='coerce')
        
        if 'diameter' not in self._df_manual.columns:
            self._df_manual['diameter'] = 0.01
        else:
            self._df_manual['diameter'] = pd.to_numeric(self._df_manual['diameter'], errors='coerce').fillna(0.01)
        
        self._df_manual.dropna(subset=['x', 'y', 'z', 'id', 'timestep'], inplace=True)

        if not self._df_manual.empty:
            x_min, x_max = self._df_manual['x'].min(), self._df_manual['x'].max()
            y_min, y_max = self._df_manual['y'].min(), self._df_manual['y'].max()
            z_min, z_max = self._df_manual['z'].min(), self._df_manual['z'].max()
            self._init_limits = ((x_min, x_max), (y_min, y_max), (z_min, z_max))

        try:
            self._df_mi_manual = self._df_manual.set_index(['timestep', 'id']).sort_index()
        except Exception:
            self._df_mi_manual = None

        if update_timesteps:
            self.timesteps = sorted(self._df_manual['timestep'].unique())
            
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
                
                # Data is now stored directly in self.sim_source by batch index.
                # We no longer perform expensive monolithic merges here.
                
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

    def get_atom_data_at_timestep(self, ts):
        if self.is_preview_mode and self.preview_df is not None:
            try:
                return self.preview_df.xs(ts, level='timestep')
            except (KeyError, TypeError):
                return self.preview_df.reset_index().set_index('id')

        if self.sim_source:
            # Trigger background load if batch is missing
            self.request_batch_for_timestep(ts)
            return self.sim_source.get_atoms_at_timestep(ts)
        
        # Fallback for manual dataframe loads (e.g. data file or manual dumps)
        if self.df_mi is None: return pd.DataFrame()
        try:
            return self.df_mi.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

    # Compatibility alias for older renderer versions
    def get_frame_data(self, ts):
        return self.get_atom_data_at_timestep(ts)

    def get_bond_data_at_timestep(self, ts):
        if self.is_preview_mode:
            return pd.DataFrame()

        if self.sim_source:
            self.request_batch_for_timestep(ts)
            return self.sim_source.get_bonds_at_timestep(ts)
            
        if self.df_bonds.empty: return pd.DataFrame()
        try:
            return self.df_bonds.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

    def get_angle_data_at_timestep(self, ts):
        if self.is_preview_mode:
            return pd.DataFrame()

        if self.sim_source:
            self.request_batch_for_timestep(ts)
            return self.sim_source.get_angles_at_timestep(ts)
            
        if self.df_angles.empty: return pd.DataFrame()
        try:
            return self.df_angles.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

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
