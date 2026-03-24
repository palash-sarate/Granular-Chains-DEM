import pandas as pd
import glob
import os
import re
import pickle
import numpy as np
from typing import Optional, Tuple

from lammps_parser import LammpsParser


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
    def __init__(self, data_dir, cache_file="sim_cache.pkl"):
        self.data_dir = data_dir
        self.chain_dump_dir = os.path.join(data_dir, "chain")
        self.cache_file = os.path.join(data_dir, cache_file)
        self.df = None

    def load_data(self, force_reload=False):
        """
        Loads data from cache if available and fresh; otherwise parses dump files.
        """
        if not force_reload and self._is_cache_valid():
            print("Loading from cache...")
            self.df = pd.read_pickle(self.cache_file)
        else:
            print("Parsing dump files (this may take a moment)...")
            self.df = self._parse_all_dumps()
            print("Saving to cache...")
            self.df.to_pickle(self.cache_file)
        
        return self.df

    def _is_cache_valid(self):
        """Checks if cache exists and is newer than the latest dump file."""
        if not os.path.exists(self.cache_file):
            return False
        print(f"Checking cache validity against dumps in: {self.chain_dump_dir}")
        dump_files = glob.glob(os.path.join(self.chain_dump_dir, "*.dump"))
        if not dump_files:
            return False
            
        latest_dump_mtime = max(os.path.getmtime(f) for f in dump_files)
        cache_mtime = os.path.getmtime(self.cache_file)
        
        return cache_mtime > latest_dump_mtime

    def _parse_single_dump(self, filepath):
        """Parses a single LAMMPS dump file."""
        # (Your existing parsing logic, optimized slightly)
        meta = {'timestep': 0}
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        for i, line in enumerate(lines):
            if "ITEM: TIMESTEP" in line:
                meta['timestep'] = int(lines[i+1])
            elif "ITEM: ATOMS" in line:
                columns = line.split()[2:]
                df = pd.read_csv(filepath, skiprows=i+1, names=columns, sep=r'\s+', engine='python')
                df['timestep'] = meta['timestep']
                return df
        return pd.DataFrame()

    def _get_step(self, filename):
        match = re.search(r'_(\d+)\.dump', filename)
        return int(match.group(1)) if match else 0
    
    def _parse_all_dumps(self):
        """Parses all dump files in directory and returns a MultiIndex DataFrame."""
        dump_files = glob.glob(os.path.join(self.chain_dump_dir, "*.dump"))
        print(f"Found {len(dump_files)} dump files to parse in {self.chain_dump_dir}")
        # Sort by timestep
        dump_files.sort(key=self._get_step)

        all_frames = [self._parse_single_dump(f) for f in dump_files]
        
        if not all_frames:
            return pd.DataFrame()

        full_df = pd.concat(all_frames)
        # Create MultiIndex (Timestep, Atom ID)
        full_df.set_index(['timestep', 'id'], inplace=True)
        full_df.sort_index(inplace=True)
        return full_df