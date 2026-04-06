import pandas as pd
import glob
import os
import re
import pickle
import numpy as np
from typing import Optional, Tuple

from concurrent.futures import ProcessPoolExecutor
from analysis.lammps_parser import LammpsParser


def detect_lammps_input_script(sim_dir: str) -> Optional[str]:
    """Auto-detect a likely LAMMPS input script in a simulation folder.

    Preference order:
    1) Files starting with ``in.`` (e.g. ``in.hopper_fill``)
    2) Files ending with ``.in``
    """
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
    """Load geometry from an auto-detected LAMMPS input script.

    Returns:
        (geometry_dict_or_none, script_path_or_none)
    """
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
    try:
        # Typical atom styles: id mol type x y z ... or id type x y z ...
        # Check first row to build column mapping
        raw_rows = np.array(data_lines, dtype=float)
        num_cols = raw_rows.shape[1]
        
        # Heuristic to detect mol vs type
        is_int_col = [np.all(raw_rows[:,i] == raw_rows[:,i].astype(int)) if i < num_cols else False for i in range(num_cols)]
        
        col_names = []
        if num_cols >= 6 and is_int_col[1] and is_int_col[2]:
             # Likely: id, mol, type, x, y, z ...
             col_names = ['id', 'mol', 'type', 'x', 'y', 'z']
        elif num_cols >= 5 and is_int_col[1]:
             # Likely: id, type, x, y, z ...
             col_names = ['id', 'type', 'x', 'y', 'z']
        else:
             # Fallback: assume minimal id, type, x, y, z
             col_names = ['id', 'type', 'x', 'y', 'z']
        
        # Append extra cols like vx, vy, vz, diameter, mass if they exist
        std_extras = ['vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'diameter', 'mass']
        for i in range(len(col_names), num_cols):
            idx = i - len(col_names)
            if idx < len(std_extras):
                col_names.append(std_extras[idx])
            else:
                col_names.append(f'v{i}')

        df = pd.DataFrame(raw_rows[:, :len(col_names)], columns=col_names)
        
        # Ensure ID, mol, type are integers
        for c in ['id', 'mol', 'type']:
            if c in df.columns:
                df[c] = df[c].astype(int)
                
    except Exception:
        # fallback: semi-brute force ID and coords
        rows = []
        for parts in data_lines:
            try:
                if len(parts) >= 4:
                    idv = int(parts[0])
                    # identify if 2nd token is likely mol
                    molv = int(parts[1]) if len(parts) > 5 else 0
                    x = float(parts[-3])
                    y = float(parts[-2])
                    z = float(parts[-1])
                    dia = float(parts[4]) if len(parts) > 5 and 'id' not in parts else 0.01 
                    rows.append((idv, molv, x, y, z, dia))
            except Exception:
                continue
        if rows:
            df = pd.DataFrame(rows, columns=['id', 'mol', 'x', 'y', 'z', 'diameter'])

    if df is None or df.empty:
        return pd.DataFrame()

    # Ensure columns exist
    for c in ['vx','vy','vz','diameter']:
        if c not in df.columns:
            df[c] = 0.0

    # attach timestep 0 for consistency with Animator expectations
    df['timestep'] = 0
    df.set_index(['timestep','id'], inplace=True)
    return df

def parse_dump_file(filepath: str) -> pd.DataFrame:
    """Public wrapper for parsing a single LAMMPS dump file.
    Returns a DataFrame with ['timestep', 'id'] as MultiIndex.
    """
    df = _parse_single_dump_fast(filepath)
    if df is not None and not df.empty:
        id_col = 'id' if 'id' in df.columns else 'index'
        if id_col in df.columns:
            # We enforce a MultiIndex to match the rest of the app's expectations
            df.set_index(['timestep', id_col], inplace=True)
            df.index.names = ['timestep', 'id']
            df.sort_index(inplace=True)
    return df

def _parse_single_dump_fast(filepath):
    """Optimized parsing of a single LAMMPS dump file using the C engine.
    Top-level function for multiprocessing compatibility.
    """
    try:
        with open(filepath, 'r') as f:
            # Quickly grab the first several lines to find structure
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

        # Use faster C engine with memory mapping
        df = pd.read_csv(filepath, 
                         skiprows=data_start_line,
                         names=cols, 
                         nrows=count if count > 0 else None,
                         sep=r'\s+', 
                         engine='c',
                         memory_map=True)
        
        # Numeric conversion
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Drop malformed rows
        id_col = 'id' if 'id' in df.columns else 'index'
        if id_col in df.columns:
            df.dropna(subset=[id_col], inplace=True)
        
        df['timestep'] = timestep
        return df
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
    return pd.DataFrame()

class SimulationData:
    def __init__(self, data_dir, cache_name="sim_cache", batch_size=100):
        self.data_dir = data_dir
        self.chain_dump_dir = os.path.join(data_dir, "chain")
        self.bond_dump_dir = os.path.join(data_dir, "bond")
        self.angle_dump_dir = os.path.join(data_dir, "angle")
        self.cache_dir = os.path.join(data_dir, cache_name)
        self.batch_size = batch_size
        
        # Ensure cache directory exists
        if os.path.isdir(self.data_dir):
            os.makedirs(self.cache_dir, exist_ok=True)
        
        self.atom_batches = {}  # Indexed by batch number
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
        """Scans directories, sorts dump files by timestep, and extracts all timesteps."""
        if os.path.isdir(self.chain_dump_dir):
            self.atom_files = sorted(glob.glob(os.path.join(self.chain_dump_dir, "*.dump")), key=self._get_step)
        if os.path.isdir(self.bond_dump_dir):
            self.bond_files = sorted(glob.glob(os.path.join(self.bond_dump_dir, "*.dump")), key=self._get_step)
        if os.path.isdir(self.angle_dump_dir):
            self.angle_files = sorted(glob.glob(os.path.join(self.angle_dump_dir, "*.dump")), key=self._get_step)

        # Build list of timesteps from atom files (or others if atom is empty)
        base_files = self.atom_files or self.bond_files or self.angle_files
        self.timesteps = [self._get_step(f) for f in base_files]
        return self.timesteps

    def get_batch_count(self):
        return max(1, (len(self.timesteps) + self.batch_size - 1) // self.batch_size) if self.timesteps else 0

    def get_available_cached_batches(self):
        """Returns a list of batch indices that have valid cache files."""
        available = []
        batch_count = self.get_batch_count()
        for i in range(batch_count):
            start_idx = i * self.batch_size
            end_idx = start_idx + self.batch_size
            
            # Form the list of source files for this batch
            batch_files = self.atom_files[start_idx:end_idx] + \
                          self.bond_files[start_idx:end_idx] + \
                          self.angle_files[start_idx:end_idx]
            
            if self._is_batch_cache_valid(i, batch_files):
                available.append(i)
        
        self.cached_batches = set(available)
        return available

    def get_batch_index_for_timestep(self, ts):
        if not self.timesteps: return 0
        try:
            # Find index in the metadata timesteps list
            idx = self.timesteps.index(ts)
            return idx // self.batch_size
        except ValueError:
            return 0

    def load_batch(self, batch_idx, force_reload=False):
        """Loads a specific batch of data."""
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
                print(f"Loading batch {batch_idx} from cache: {cache_file}")
                with open(cache_file, 'rb') as f:
                    loaded_data = pickle.load(f)
            except Exception as e:
                print(f"Cache loading failed for batch {batch_idx} ({e}), re-parsing...")

        if loaded_data is None:
            print(f"Parsing dump files for batch {batch_idx}...")
            batch_atoms = self._parse_dump_list(batch_atoms_files)
            batch_bonds = self._parse_dump_list(batch_bonds_files)
            batch_angles = self._parse_dump_list(batch_angles_files)
            
            loaded_data = {
                'atoms': batch_atoms,
                'bonds': batch_bonds,
                'angles': batch_angles
            }
            print(f"Saving batch {batch_idx} to cache...")
            try:
                with open(cache_file, 'wb') as f:
                    pickle.dump(loaded_data, f)
            except Exception as e:
                print(f"Could not save cache: {e}")

        # Store loaded data by batch index
        if not loaded_data['atoms'].empty:
            self.atom_batches[batch_idx] = loaded_data['atoms']
            
        if not loaded_data['bonds'].empty:
            self.bond_batches[batch_idx] = loaded_data['bonds']
            
        if not loaded_data['angles'].empty:
            self.angle_batches[batch_idx] = loaded_data['angles']

        self.loaded_batches.add(batch_idx)
        return loaded_data

    def get_atoms_at_timestep(self, ts):
        batch_idx = self.get_batch_index_for_timestep(ts)
        batch = self.atom_batches.get(batch_idx)
        if batch is None or batch.empty: return pd.DataFrame()
        try:
            return batch.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

    def get_bonds_at_timestep(self, ts):
        batch_idx = self.get_batch_index_for_timestep(ts)
        batch = self.bond_batches.get(batch_idx)
        if batch is None or batch.empty: return pd.DataFrame()
        try:
            return batch.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

    def get_angles_at_timestep(self, ts):
        batch_idx = self.get_batch_index_for_timestep(ts)
        batch = self.angle_batches.get(batch_idx)
        if batch is None or batch.empty: return pd.DataFrame()
        try:
            return batch.xs(ts, level='timestep')
        except (KeyError, TypeError):
            return pd.DataFrame()

    def load_data(self, force_reload=False):
        """Metadata-first loading. Returns only the first batch by default."""
        self.load_metadata()
        if self.get_batch_count() > 0:
            return self.load_batch(0, force_reload=force_reload)
        return {'atoms': pd.DataFrame(), 'bonds': pd.DataFrame(), 'angles': pd.DataFrame()}

    def _is_batch_cache_valid(self, batch_idx, batch_files):
        """Checks if batch cache exists and is newer than the dump files in it."""
        cache_file = self._get_cache_file(batch_idx)
        if not os.path.exists(cache_file):
            return False
        
        if not batch_files:
            return False
            
        latest_dump_mtime = max((os.path.getmtime(f) for f in batch_files if os.path.exists(f)), default=0)
        if latest_dump_mtime == 0:
            return False
            
        cache_mtime = os.path.getmtime(cache_file)
        return cache_mtime > latest_dump_mtime


    def _get_step(self, filename):
        match = re.search(r'_(\d+)\.dump', filename)
        return int(match.group(1)) if match else 0
    
    def _parse_dump_list(self, dump_files):
        """Parses a list of dump files in parallel and returns a MultiIndex DataFrame."""
        if not dump_files:
            return pd.DataFrame()
            
        print(f"Parsing {len(dump_files)} dump files in parallel...")

        # Use ProcessPoolExecutor for true parallelism in Python
        # Windows requires protection but since this is called from within the app 
        # it should be safe as long as we're not spawning recursive processes.
        with ProcessPoolExecutor() as executor:
            all_frames = list(executor.map(_parse_single_dump_fast, dump_files))
        
        all_frames = [f for f in all_frames if f is not None and not f.empty]
        
        if not all_frames:
            return pd.DataFrame()

        try:
            full_df = pd.concat(all_frames)
            id_col = 'id' if 'id' in full_df.columns else 'index'
            if id_col not in full_df.columns:
                return pd.DataFrame()
                
            full_df.set_index(['timestep', id_col], inplace=True)
            full_df.index.names = ['timestep', 'id']
            full_df.sort_index(inplace=True)
            return full_df
        except Exception as e:
            print(f"Failed to assemble DataFrames: {e}")
            return pd.DataFrame()