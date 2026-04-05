import pandas as pd
import glob
import os
import re
import pickle
import numpy as np
from typing import Optional, Tuple

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

class SimulationData:
    def __init__(self, data_dir, cache_prefix="sim_cache", batch_size=1000):
        self.data_dir = data_dir
        self.chain_dump_dir = os.path.join(data_dir, "chain")
        self.bond_dump_dir = os.path.join(data_dir, "bond")
        self.angle_dump_dir = os.path.join(data_dir, "angle")
        self.cache_prefix = os.path.join(data_dir, cache_prefix)
        self.batch_size = batch_size
        
        self.df_atoms = pd.DataFrame()
        self.df_bonds = pd.DataFrame()
        self.df_angles = pd.DataFrame()
        
        self.timesteps = []
        self.atom_files = []
        self.bond_files = []
        self.angle_files = []
        self.loaded_batches = set()

    def _get_cache_file(self, batch_idx):
        return f"{self.cache_prefix}_batch_{batch_idx}.pkl"

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
            return {'atoms': self.df_atoms, 'bonds': self.df_bonds, 'angles': self.df_angles}

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

        # Accumulate loaded data
        if not loaded_data['atoms'].empty:
            self.df_atoms = pd.concat([self.df_atoms, loaded_data['atoms']]) if not self.df_atoms.empty else loaded_data['atoms']
            # Re-sort to maintain clean multiindex
            self.df_atoms = self.df_atoms[~self.df_atoms.index.duplicated(keep='last')].sort_index()

        if not loaded_data['bonds'].empty:
            self.df_bonds = pd.concat([self.df_bonds, loaded_data['bonds']]) if not self.df_bonds.empty else loaded_data['bonds']
            self.df_bonds = self.df_bonds[~self.df_bonds.index.duplicated(keep='last')].sort_index()

        if not loaded_data['angles'].empty:
            self.df_angles = pd.concat([self.df_angles, loaded_data['angles']]) if not self.df_angles.empty else loaded_data['angles']
            self.df_angles = self.df_angles[~self.df_angles.index.duplicated(keep='last')].sort_index()

        self.loaded_batches.add(batch_idx)
        return {'atoms': self.df_atoms, 'bonds': self.df_bonds, 'angles': self.df_angles}

    def load_data(self, force_reload=False):
        """Metadata-first loading. Returns only the first batch by default."""
        self.load_metadata()
        if self.get_batch_count() > 0:
            return self.load_batch(0, force_reload=force_reload)
        return {'atoms': self.df_atoms, 'bonds': self.df_bonds, 'angles': self.df_angles}

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

    def _parse_single_dump(self, filepath):
        """Parses a single LAMMPS dump file (atoms, bonds, or angles)."""
        meta = {'timestep': 0}
        count = 0
        try:
            with open(filepath, 'r') as f:
                curr_line_idx = -1
                while True:
                    line = f.readline()
                    if not line:
                        break
                    curr_line_idx += 1
                    
                    if "ITEM: TIMESTEP" in line:
                        step_line = f.readline()
                        curr_line_idx += 1
                        if step_line:
                            try:
                                meta['timestep'] = int(step_line.strip())
                            except ValueError:
                                pass
                    elif "ITEM: NUMBER OF" in line:
                        count_line = f.readline()
                        curr_line_idx += 1
                        if count_line:
                            try:
                                count = int(count_line.strip())
                            except ValueError:
                                count = 0
                    elif "ITEM: ATOMS" in line or "ITEM: ENTRIES" in line:
                        columns = line.split()[2:]
                        # SKIP the headers we already read + the ITEM: ATOMS/ENTRIES header itself
                        df = pd.read_csv(filepath, skiprows=curr_line_idx + 1,
                                         names=columns, 
                                         nrows=count if count > 0 else None,
                                         sep=r'\s+', engine='python')
                        
                        # Defensive: drop any rows that failed to parse
                        for col in df.columns:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                        df.dropna(subset=df.columns.intersection(['x','y','z','id','index','dist','theta']), 
                                  inplace=True)
                        
                        df['timestep'] = meta['timestep']
                        return df
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
        return pd.DataFrame()

    def _get_step(self, filename):
        match = re.search(r'_(\d+)\.dump', filename)
        return int(match.group(1)) if match else 0
    
    def _parse_dump_list(self, dump_files):
        """Parses a specific list of dump files and returns a MultiIndex DataFrame."""
        if not dump_files:
            return pd.DataFrame()
            
        print(f"Parsing {len(dump_files)} dump files...")

        all_frames = [self._parse_single_dump(f) for f in dump_files]
        all_frames = [f for f in all_frames if not f.empty]
        
        if not all_frames:
            return pd.DataFrame()

        try:
            full_df = pd.concat(all_frames)
            # Create MultiIndex (Timestep, ID). 
            # For atoms it is 'id', for bonds/angles it is 'index'.
            id_col = 'id' if 'id' in full_df.columns else 'index'
            if id_col not in full_df.columns:
                print(f"Warning: Neither 'id' nor 'index' found in dumps.")
                return pd.DataFrame()
                
            full_df.set_index(['timestep', id_col], inplace=True)
            full_df.index.names = ['timestep', 'id']  # Unified index names
            full_df.sort_index(inplace=True)
            return full_df
        except Exception as e:
            print(f"Failed to assemble DataFrames: {e}")
            return pd.DataFrame()