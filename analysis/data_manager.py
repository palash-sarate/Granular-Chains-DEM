import pandas as pd
import glob
import os
import re
import pickle

class SimulationData:
    def __init__(self, data_dir, cache_file="sim_cache.pkl"):
        self.data_dir = data_dir
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
        
        dump_files = glob.glob(os.path.join(self.data_dir, "chain_*.dump"))
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
                df = pd.read_csv(filepath, skiprows=i+1, names=columns, sep='\s+', engine='python')
                df['timestep'] = meta['timestep']
                return df
        return pd.DataFrame()

    def _parse_all_dumps(self):
        """Parses all dump files in directory and returns a MultiIndex DataFrame."""
        dump_files = glob.glob(os.path.join(self.data_dir, "*.dump"))
        
        # Sort by timestep
        def get_step(filename):
            match = re.search(r'_(\d+)\.dump', filename)
            return int(match.group(1)) if match else 0
        dump_files.sort(key=get_step)

        all_frames = [self._parse_single_dump(f) for f in dump_files]
        
        if not all_frames:
            return pd.DataFrame()

        full_df = pd.concat(all_frames)
        # Create MultiIndex (Timestep, Atom ID)
        full_df.set_index(['timestep', 'id'], inplace=True)
        full_df.sort_index(inplace=True)
        return full_df