import pandas as pd
import glob
import os
import re
import pickle
import numpy as np
import json
from typing import Optional, Tuple, Any, Dict

from concurrent.futures import ProcessPoolExecutor
from analysis.lammps_parser import LammpsParser


def find_simulation_lineage(sim_dir: str) -> list:
    """Finds the full lineage (ancestors and descendants) for a simulation."""
    # 1. Walk backward to find the root
    ancestors = []
    curr = os.path.abspath(sim_dir)
    visited = set()
    while curr and os.path.isdir(curr) and curr not in visited:
        ancestors.append(curr)
        visited.add(curr)
        
        # Try to find metadata (prioritize unified name)
        meta_names = ["metadata.json", "grid_metadata.json", "sim_metadata.json"]
        metadata_path = None
        for name in meta_names:
            p = os.path.join(curr, name)
            if os.path.exists(p):
                metadata_path = p
                break
            
        if not metadata_path:
            break
            
        try:
            with open(metadata_path, 'r') as f:
                meta = json.load(f)
            
            parent = meta.get("source_dir") or meta.get("parent_run")
            if not parent and "restart_path" in meta:
                from pathlib import Path
                restart_p = Path(meta["restart_path"])
                # If it's a file path, we want the folder containing it
                if os.path.isfile(str(restart_p)):
                    parent = str(restart_p.parent)
                else:
                    # If it's a folder, it's the parent
                    parent = str(restart_p)
                    
            if parent:
                parent_abs = os.path.abspath(parent)
                # If it points to the 'chain' subfolder, go up one
                if os.path.basename(parent_abs) == 'chain':
                    parent_abs = os.path.dirname(parent_abs)
                    
                if os.path.isdir(parent_abs) and parent_abs not in visited:
                    curr = parent_abs
                else:
                    break
            else:
                break
        except:
            break
            
    # ancestors is currently [current, parent, grandparent...]
    # root is the last one
    root = ancestors[-1]
    all_related = set(ancestors)
    
    # 2. Walk forward to find descendants.
    # We check the siblings of the root AND scan dumping_yard subdirectories.
    search_dirs = []
    root_parent = os.path.dirname(root)
    if os.path.isdir(root_parent):
        search_dirs.append(root_parent)
    
    dy = os.path.abspath("dumping_yard")
    if os.path.isdir(dy) and dy != root_parent:
        search_dirs.append(dy)

    # Collect all potential candidate folders
    candidates = []
    for s_dir in search_dirs:
        try:
            for entry in os.scandir(s_dir):
                if entry.is_dir():
                    candidates.append(entry.path)
                    # Also check one level deep (e.g. dumping_yard/SimName/RunName)
                    try:
                        for sub in os.scandir(entry.path):
                            if sub.is_dir():
                                candidates.append(sub.path)
                    except: pass
        except: pass

    # Iteratively find children using a work-queue
    found_new = True
    processed_candidates = set()
    while found_new:
        found_new = False
        for cand in candidates:
            if cand in all_related or cand in processed_candidates: continue
            
            meta_names = ["metadata.json", "grid_metadata.json", "sim_metadata.json"]
            meta_path = None
            for name in meta_names:
                p = os.path.join(cand, name)
                if os.path.exists(p):
                    meta_path = p
                    break
                    
            if meta_path:
                try:
                    with open(meta_path, 'r') as f:
                        meta = json.load(f)
                    src = meta.get("source_dir") or meta.get("parent_run")
                    if src:
                        src_abs = os.path.abspath(src)
                        # Normalize common restart directory patterns
                        if os.path.basename(src_abs) == 'chain':
                            src_abs = os.path.dirname(src_abs)
                            
                        if src_abs in all_related:
                            all_related.add(os.path.abspath(cand))
                            found_new = True
                except:
                    pass
            processed_candidates.add(cand)

    # Sort chronologically by creation time
    lineage = sorted(list(all_related), key=lambda x: os.path.getctime(x))
    return lineage

def aggregate_lifecycle_metadata(lineage: list) -> dict:
    """Merges metadata from all runs in the lineage, oldest to newest."""
    aggregated = {}
    for d in lineage:
        meta_names = ["metadata.json", "grid_metadata.json", "sim_metadata.json"]
        metadata_path = None
        for name in meta_names:
            p = os.path.join(d, name)
            if os.path.exists(p):
                metadata_path = p
                break
        
        if metadata_path:
            try:
                with open(metadata_path, 'r') as f:
                    meta = json.load(f)
                aggregated.update(meta)
            except:
                pass
    return aggregated

def detect_lammps_input_script(sim_dir: str) -> Optional[str]:
    """Auto-detect a likely LAMMPS input script in a simulation folder."""
    if not sim_dir or not os.path.isdir(sim_dir):
        return None

    in_dot = sorted(glob.glob(os.path.join(sim_dir, "in.*")))
    for path in in_dot:
        if os.path.isfile(path):
            return path

    dot_in = sorted(glob.glob(os.path.join(sim_dir, "*.in")))
    for path in dot_in:
        if os.path.isfile(path):
            return path

    return None

def load_lammps_geometry(sim_dir: str) -> Tuple[Optional[dict], Optional[str]]:
    """Load geometry from an auto-detected LAMMPS input script."""
    script_path = detect_lammps_input_script(sim_dir)
    if not script_path:
        return None, None

    try:
        return parse_lammps_geometry_script(script_path), script_path
    except Exception as e:
        print(f"Failed to parse LAMMPS geometry from {script_path}: {e}")
        return None, script_path

def parse_lammps_geometry_script(script_path: str) -> Optional[dict]:
    """Parse geometry directly from a specific LAMMPS input script path."""
    if not script_path or not os.path.isfile(script_path):
        return None
    parser = LammpsParser(script_path)
    return parser.get_geometry()

def parse_simple_data_file(path: str, manual_map: dict = None) -> pd.DataFrame:
    """Parses a single LAMMPS data file, identifying columns by heuristics or manual map."""
    if not os.path.exists(path):
        return pd.DataFrame()

    with open(path, 'r') as f:
        lines = f.readlines()

    # Find the Atoms section
    atom_start = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atom_start = i + 1
            break
            
    if atom_start == -1:
        return pd.DataFrame()

    # Read atom data until blank line or next section
    rows = []
    SECTION_HEADERS = set(['Velocities', 'Bonds', 'Angles', 'Masses', 'PairIJ'])
    for i in range(atom_start, len(lines)):
        line = lines[i].strip()
        if not line:
            if rows: break
            continue
        if line.split()[0] in SECTION_HEADERS:
            break
        parts = line.split()
        if len(parts) < 3: continue
        rows.append(parts)

    if not rows:
        return pd.DataFrame()

    raw_rows = np.array(rows, dtype=object)
    num_cols = raw_rows.shape[1]

    if manual_map:
        # manual_map looks like {'x': 2, 'y': 3, 'z': 4, ...}
        df_cols = {}
        for name, idx in manual_map.items():
            if idx < num_cols:
                try:
                    df_cols[name] = pd.to_numeric(raw_rows[:, idx], errors='coerce')
                except Exception: pass
        df = pd.DataFrame(df_cols)
    else:
        # Use heuristics
        try:
            data_numeric = raw_rows.astype(float)
            is_int_col = [np.all(data_numeric[:,i] == data_numeric[:,i].astype(int)) for i in range(num_cols)]
            
            col_names = []
            first_float_idx = -1
            for i in range(2, min(5, num_cols)):
                if not is_int_col[i]:
                    first_float_idx = i
                    break
            
            if first_float_idx == 2:
                col_names = ['id', 'type', 'x', 'y', 'z']
            elif first_float_idx == 3:
                col_names = ['id', 'mol', 'type', 'x', 'y', 'z']
            elif first_float_idx == -1 and num_cols >= 5:
                if is_int_col[2] and is_int_col[3] and is_int_col[4] and num_cols >= 6:
                    col_names = ['id', 'mol', 'type', 'x', 'y', 'z']
                else:
                    col_names = ['id', 'type', 'x', 'y', 'z']
            else:
                col_names = ['id', 'type', 'x', 'y', 'z']

            remaining = list(range(len(col_names), num_cols))
            for idx in remaining:
                if np.all(data_numeric[:, idx] > 0) and np.all(data_numeric[:, idx] < 0.5) and np.unique(data_numeric[:, idx]).size == 1:
                    col_names.append('diameter')
                else:
                    col_names.append(f'v{idx}')
                
            df = pd.DataFrame(data_numeric[:, :len(col_names)], columns=col_names)
            for c in ['id', 'mol', 'type']:
                if c in df.columns:
                    df[c] = df[c].astype(int)
        except Exception:
            return pd.DataFrame()

    if df.empty:
        return pd.DataFrame()

    for c in ['vx','vy','vz','diameter']:
        if c not in df.columns:
            df[c] = 0.0 if c != 'diameter' else 0.01

    df['timestep'] = 0
    if 'id' in df.columns:
        df.set_index(['timestep','id'], inplace=True)
    else:
        df['id'] = range(len(df))
        df.set_index(['timestep','id'], inplace=True)
    return df

def parse_dump_file(filepath: str) -> pd.DataFrame:
    """Public wrapper for parsing a single LAMMPS dump file."""
    df = _parse_single_dump_fast(filepath)
    if df is not None and not df.empty:
        id_col = 'id' if 'id' in df.columns else 'index'
        if id_col in df.columns:
            df.set_index(['timestep', id_col], inplace=True)
            df.index.names = ['timestep', 'id']
            df.sort_index(inplace=True)
    return df

def _parse_single_dump_fast(filepath):
    """Optimized parsing of a single LAMMPS dump file using the C engine."""
    try:
        with open(filepath, 'r') as f:
            header_lines = [f.readline() for _ in range(15)]
        
        timestep = 0
        count = 0
        data_start_line = 0
        cols = []

        for i, line in enumerate(header_lines):
            if "ITEM: TIMESTEP" in line:
                try: timestep = int(header_lines[i+1].strip())
                except: pass
            elif "ITEM: NUMBER OF" in line:
                try: count = int(header_lines[i+1].strip())
                except: pass
            elif "ITEM: ATOMS" in line or "ITEM: ENTRIES" in line:
                cols = line.split()[2:]
                data_start_line = i + 1
                break
        
        if not cols: return pd.DataFrame()

        df = pd.read_csv(filepath, 
                         skiprows=data_start_line,
                         names=cols, 
                         nrows=count if count > 0 else None,
                         sep=r'\s+', 
                         engine='c',
                         memory_map=True)
        
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        id_col = 'id' if 'id' in df.columns else 'index'
        if id_col in df.columns:
            df.dropna(subset=[id_col], inplace=True)
        
        df['timestep'] = timestep
        return df
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
    return pd.DataFrame()

class SimulationMetadata:
    """Manages persistent metadata for a simulation folder."""
    def __init__(self, data_dir: str, filename: str = "metadata.json"):
        self.data_dir = data_dir
        self.path = os.path.join(data_dir, filename)
        self._data: Dict[str, Any] = {}
        self.load()

    def load(self):
        """Loads and merges metadata, supporting migration from legacy filenames."""
        merged = {}
        # Order of priority: primary file > grid_metadata > sim_metadata
        legacy_names = ["metadata.json", "grid_metadata.json", "sim_metadata.json"]
        
        for name in reversed(legacy_names):
            p = os.path.join(self.data_dir, name)
            if os.path.exists(p):
                try:
                    with open(p, 'r') as f:
                        data = json.load(f)
                    merged.update(data)
                except Exception as e:
                    print(f"Failed to load metadata from {p}: {e}")
        
        self._data = merged

    def save(self):
        try:
            os.makedirs(os.path.dirname(self.path), exist_ok=True)
            with open(self.path, 'w') as f:
                json.dump(self._data, f, indent=4)
        except Exception as e:
            print(f"Failed to save metadata to {self.path}: {e}")

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    def set(self, key: str, value: Any):
        self._data[key] = value
        self.save()

    def update(self, delta: Dict[str, Any]):
        self._data.update(delta)
        self.save()

    def __getitem__(self, key):
        return self._data.get(key)

    def __setitem__(self, key, value):
        self.set(key, value)

class SimulationData:
    def __init__(self, data_dirs, cache_name="sim_cache", batch_size=100):
        if isinstance(data_dirs, (str, bytes)):
            data_dirs = [data_dirs]
        self.data_dirs = [os.path.abspath(d) for d in data_dirs]
        self.data_dir = self.data_dirs[-1] # Primary directory for writing/active metadata
        
        self.lineage = find_simulation_lineage(self.data_dir)
        self.full_metadata = aggregate_lifecycle_metadata(self.lineage)
        
        # Determine cache location
        if len(self.data_dirs) > 1:
            import hashlib
            paths_str = "|".join(sorted(self.data_dirs))
            h = hashlib.md5(paths_str.encode()).hexdigest()[:12]
            # Place composite cache in a sibling directory to the root of the lineage
            self.cache_dir = os.path.join(os.path.dirname(self.data_dirs[0]), ".composite_cache", h)
        else:
            self.cache_dir = os.path.join(self.data_dir, cache_name)
            
        self.batch_size = batch_size
        
        if os.path.isdir(self.data_dir):
            os.makedirs(self.cache_dir, exist_ok=True)
            self.metadata = SimulationMetadata(self.data_dir)
        else:
            self.metadata = None
        
        self.atom_batches = {}
        self.bond_batches = {}
        self.angle_batches = {}
        
        self.timesteps = []
        self.atom_files = []
        self.bond_files = []
        self.angle_files = []
        self.loaded_batches = set()
        self.cached_batches = set()

    def _get_cache_file(self, batch_idx):
        return os.path.join(self.cache_dir, f"batch_{batch_idx}.pkl")

    def load_metadata(self):
        """Builds a unified timestep list by merging all directories (latest wins on overlap)."""
        step_to_atom = {}
        step_to_bond = {}
        step_to_angle = {}
        
        # Process chronologically - later directories in list overwrite earlier ones for same timestep
        for d in self.data_dirs:
            # Atom dumps
            ad = os.path.join(d, "chain")
            if os.path.isdir(ad):
                files = glob.glob(os.path.join(ad, "*.dump"))
                for f in files:
                    step_to_atom[self._get_step(f)] = f
            # Bond dumps
            bd = os.path.join(d, "bond")
            if os.path.isdir(bd):
                files = glob.glob(os.path.join(bd, "*.dump"))
                for f in files:
                    step_to_bond[self._get_step(f)] = f
            # Angle dumps
            andir = os.path.join(d, "angle")
            if os.path.isdir(andir):
                files = glob.glob(os.path.join(andir, "*.dump"))
                for f in files:
                    step_to_angle[self._get_step(f)] = f
                    
        # Synchronize timesteps (usually they align, but we take atoms as the master clock)
        self.timesteps = sorted(step_to_atom.keys())
        self.atom_files = [step_to_atom[ts] for ts in self.timesteps]
        self.bond_files = [step_to_bond.get(ts) for ts in self.timesteps]
        self.angle_files = [step_to_angle.get(ts) for ts in self.timesteps]
        
        return self.timesteps

    def get_available_cached_batches(self):
        """Returns a list of batch indices that have valid caches on disk."""
        available = []
        batch_count = self.get_batch_count()
        for i in range(batch_count):
            start_idx = i * self.batch_size
            end_idx = start_idx + self.batch_size
            batch_files = self.atom_files[start_idx:end_idx] + \
                          self.bond_files[start_idx:end_idx] + \
                          self.angle_files[start_idx:end_idx]
            if self._is_batch_cache_valid(i, batch_files):
                available.append(i)
        self.cached_batches = set(available)
        return available

    def get_batch_count(self):
        return max(1, (len(self.timesteps) + self.batch_size - 1) // self.batch_size) if self.timesteps else 0

    def get_batch_index_for_timestep(self, ts):
        if not self.timesteps: return 0
        try:
            idx = self.timesteps.index(ts)
            return idx // self.batch_size
        except ValueError:
            return 0

    def load_batch(self, batch_idx, force_reload=False):
        if batch_idx in self.loaded_batches and not force_reload:
            return {
                'atoms': self.atom_batches.get(batch_idx, pd.DataFrame()),
                'bonds': self.bond_batches.get(batch_idx, pd.DataFrame()),
                'angles': self.angle_batches.get(batch_idx, pd.DataFrame())
            }

        cache_file = self._get_cache_file(batch_idx)
        start_idx = batch_idx * self.batch_size
        end_idx = start_idx + self.batch_size

        batch_atoms_files = self.atom_files[start_idx:end_idx]
        batch_bonds_files = self.bond_files[start_idx:end_idx]
        batch_angles_files = self.angle_files[start_idx:end_idx]

        loaded_data = None
        if not force_reload and self._is_batch_cache_valid(batch_idx, batch_atoms_files + batch_bonds_files + batch_angles_files):
            try:
                with open(cache_file, 'rb') as f:
                    loaded_data = pickle.load(f)
            except Exception: pass

        if loaded_data is None:
            batch_atoms = self._parse_dump_list(batch_atoms_files)
            batch_bonds = self._parse_dump_list(batch_bonds_files)
            batch_angles = self._parse_dump_list(batch_angles_files)
            
            loaded_data = {'atoms': batch_atoms, 'bonds': batch_bonds, 'angles': batch_angles}
            try:
                with open(cache_file, 'wb') as f:
                    pickle.dump(loaded_data, f)
            except Exception: pass

        if not loaded_data['atoms'].empty: self.atom_batches[batch_idx] = loaded_data['atoms']
        if not loaded_data['bonds'].empty: self.bond_batches[batch_idx] = loaded_data['bonds']
        if not loaded_data['angles'].empty: self.angle_batches[batch_idx] = loaded_data['angles']

        self.loaded_batches.add(batch_idx)
        return loaded_data

    def get_atoms_at_timestep(self, ts):
        batch_idx = self.get_batch_index_for_timestep(ts)
        batch = self.atom_batches.get(batch_idx)
        if batch is None or batch.empty: return pd.DataFrame()
        try: return batch.xs(ts, level='timestep')
        except: return pd.DataFrame()

    def load_data(self, force_reload=False):
        self.load_metadata()
        if self.get_batch_count() > 0:
            return self.load_batch(0, force_reload=force_reload)
        return {'atoms': pd.DataFrame(), 'bonds': pd.DataFrame(), 'angles': pd.DataFrame()}

    def _is_batch_cache_valid(self, batch_idx, batch_files):
        cache_file = self._get_cache_file(batch_idx)
        valid_files = [f for f in batch_files if f and os.path.exists(f)]
        if not os.path.exists(cache_file) or not valid_files: return False
        latest_dump_mtime = max((os.path.getmtime(f) for f in valid_files), default=0)
        return os.path.getmtime(cache_file) > latest_dump_mtime

    def _get_step(self, filename):
        match = re.search(r'_(\d+)\.dump', filename)
        return int(match.group(1)) if match else 0
    
    def _parse_dump_list(self, dump_files):
        valid_files = [f for f in dump_files if f]
        if not valid_files: return pd.DataFrame()
        with ProcessPoolExecutor() as executor:
            all_frames = list(executor.map(_parse_single_dump_fast, valid_files))
        all_frames = [f for f in all_frames if f is not None and not f.empty]
        if not all_frames: return pd.DataFrame()
        try:
            full_df = pd.concat(all_frames)
            id_col = 'id' if 'id' in full_df.columns else 'index'
            if id_col not in full_df.columns: return pd.DataFrame()
            full_df.set_index(['timestep', id_col], inplace=True)
            full_df.index.names = ['timestep', 'id']
            full_df.sort_index(inplace=True)
            return full_df
        except: return pd.DataFrame()