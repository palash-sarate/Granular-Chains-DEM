import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time
from typing import List, Dict, Any, Callable, Tuple, Optional, Union

# --- Data Structures ---

class SimulationFrame:
    """
    Represents a single snapshot of the simulation.
    """
    def __init__(self, timestep: int, atoms: Dict[int, Dict[str, float]] = None, 
                 bonds: Dict[int, Dict[str, float]] = None, 
                 angles: Dict[int, Dict[str, float]] = None):
        self.timestep = timestep
        self.atoms = atoms or {}
        self.bonds = bonds or {}
        self.angles = angles or {}

class Visualizer:
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the Visualizer with a configuration dictionary.
        
        Config Structure:
        {
            "layout": (rows, cols),
            "figsize": (width, height),
            "title": "Window Title",
            "files": {
                "atoms": "path/to/atom.dump",
                "bonds": "path/to/bond.dump",
                "angles": "path/to/angle.dump"
            },
            "plots": [
                {
                    "position": (row, col), # 0-indexed
                    "title": "Plot Title",
                    "xlabel": "Time",
                    "ylabel": "Value",
                    "series": [
                        {
                            "label": "Series Label",
                            "type": "atom" | "bond" | "angle" | "function",
                            "id": int, # ID of atom/bond/angle
                            "var": "x" | "y" | "z" | "length" | "angle" | ...,
                            "func": Callable[[SimulationFrame], float], # For type='function'
                            "style": {"color": "red", "linestyle": "--"} # Optional matplotlib kwargs
                        },
                        ...
                    ]
                },
                ...
            ]
        }
        """
        self.config = config
        self.layout = config.get("layout", (1, 1))
        self.figsize = config.get("figsize", (12, 8))
        self.title = config.get("title", "Simulation Visualizer")
        
        # Setup Matplotlib
        plt.ion()
        self.fig, self.axes = plt.subplots(self.layout[0], self.layout[1], figsize=self.figsize)
        self.fig.suptitle(self.title)
        
        # Normalize axes to 1D array
        if self.layout[0] * self.layout[1] == 1:
            self.axes = np.array([self.axes])
        self.axes = self.axes.flatten()
        
        # Internal state
        self.lines = {} # Map key -> Line2D
        self.data_history = {} # Map key -> list of values
        self.time_history = []
        self.file_pointers = {} # Map file_type -> file position
        
        self._setup_plots()
        
    def _setup_plots(self):
        """Initialize plot axes and lines based on config."""
        plots_config = self.config.get("plots", [])
        
        for plot_cfg in plots_config:
            row, col = plot_cfg.get("position", (0, 0))
            ax_idx = row * self.layout[1] + col
            
            if ax_idx >= len(self.axes):
                print(f"Warning: Plot position ({row}, {col}) is out of layout bounds.")
                continue
                
            ax = self.axes[ax_idx]
            ax.set_title(plot_cfg.get("title", ""))
            ax.set_xlabel(plot_cfg.get("xlabel", "Time"))
            ax.set_ylabel(plot_cfg.get("ylabel", "Value"))
            ax.grid(True, alpha=0.3)
            
            for i, series in enumerate(plot_cfg.get("series", [])):
                label = series.get("label", f"Series {i+1}")
                # Unique key for this series
                key = f"{ax_idx}_{i}_{label}"
                
                style = series.get("style", {})
                line, = ax.plot([], [], label=label, **style)
                
                self.lines[key] = line
                self.data_history[key] = []
                
                # Store series config with the key for easy lookup during update
                series['_key'] = key
            
            ax.legend()

    def update(self, frame: SimulationFrame):
        """
        Update all plots with data from the given SimulationFrame.
        """
        self.time_history.append(frame.timestep)
        
        plots_config = self.config.get("plots", [])
        
        for plot_cfg in plots_config:
            row, col = plot_cfg.get("position", (0, 0))
            ax_idx = row * self.layout[1] + col
            if ax_idx >= len(self.axes): continue
            
            ax = self.axes[ax_idx]
            updated_any = False
            
            for series in plot_cfg.get("series", []):
                key = series.get('_key')
                if not key: continue
                
                val = self._evaluate_series(series, frame)
                
                # Append data
                self.data_history[key].append(val)
                
                # Update line
                # Handle cases where history lengths might mismatch if update failed previously (unlikely here)
                self.lines[key].set_data(self.time_history, self.data_history[key])
                updated_any = True
            
            if updated_any:
                ax.relim()
                ax.autoscale_view()
        
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def _evaluate_series(self, series_cfg: Dict, frame: SimulationFrame) -> float:
        """Extracts value for a single series from the frame."""
        s_type = series_cfg.get("type", "atom")
        
        try:
            if s_type == "function":
                func = series_cfg.get("func")
                if callable(func):
                    return func(frame)
                return np.nan
                
            # Data lookup types
            obj_id = series_cfg.get("id")
            var_name = series_cfg.get("var")
            
            data_source = None
            if s_type == "atom":
                data_source = frame.atoms
            elif s_type == "bond":
                data_source = frame.bonds
            elif s_type == "angle":
                data_source = frame.angles
            
            if data_source and obj_id in data_source:
                return data_source[obj_id].get(var_name, np.nan)
            
            return np.nan
            
        except Exception as e:
            # print(f"Error evaluating series {series_cfg.get('label')}: {e}")
            return np.nan

    def sync_from_files(self):
        """
        Reads the latest available frame from configured files and updates the plots.
        """
        files_cfg = self.config.get("files", {})
        
        # We need at least one file to get a timestep
        if not files_cfg:
            return

        # Parse latest frame from each file
        # Note: This assumes all files are synchronized to the same timestep
        # or we just take the latest from each.
        
        atoms_data, timestep = self._read_latest_dump(files_cfg.get("atoms"), "atoms")
        bonds_data, _ = self._read_latest_dump(files_cfg.get("bonds"), "bonds")
        angles_data, _ = self._read_latest_dump(files_cfg.get("angles"), "angles")
        
        if timestep is not None:
            # Only update if we have a valid timestep (from atoms usually)
            frame = SimulationFrame(timestep, atoms=atoms_data, bonds=bonds_data, angles=angles_data)
            
            # Check if this timestep is new
            if not self.time_history or timestep > self.time_history[-1]:
                self.update(frame)

    def _read_latest_dump(self, filepath: str, file_type: str) -> Tuple[Dict, Optional[int]]:
        """
        Reads the last complete frame from a LAMMPS dump file.
        Returns (data_dict, timestep).
        """
        if not filepath or not os.path.exists(filepath):
            return {}, None

        # Naive implementation: Read whole file (slow for large files)
        # Optimized: Read last N bytes and parse backwards
        
        try:
            with open(filepath, 'rb') as f:
                # Seek to end
                f.seek(0, os.SEEK_END)
                file_size = f.tell()
                
                # Read last 50KB (adjust based on system size)
                read_size = min(file_size, 100 * 1024) 
                f.seek(-read_size, os.SEEK_END)
                content = f.read().decode('utf-8', errors='ignore')
                
            # Find last "ITEM: TIMESTEP"
            # We split by "ITEM: TIMESTEP" and take the last chunk
            chunks = content.split("ITEM: TIMESTEP")
            if len(chunks) < 2:
                return {}, None
            
            last_chunk = chunks[-1].strip().split('\n')
            if len(last_chunk) < 4: # Header needs at least timestep, num atoms, box, header
                # Try previous chunk if last one is incomplete
                if len(chunks) >= 3:
                    last_chunk = chunks[-2].strip().split('\n')
                else:
                    return {}, None

            # Parse Timestep
            try:
                timestep = int(last_chunk[0])
            except ValueError:
                return {}, None

            # Parse Body
            # Look for "ITEM: ATOMS" or "ITEM: ENTRIES"
            header_idx = -1
            for i, line in enumerate(last_chunk):
                if line.startswith("ITEM: ATOMS") or line.startswith("ITEM: ENTRIES"):
                    header_idx = i
                    break
            
            if header_idx == -1:
                return {}, timestep

            # Parse columns
            # Line format: ITEM: ATOMS id type x y z ...
            headers = last_chunk[header_idx].split()[2:]
            
            data = {}
            for line in last_chunk[header_idx+1:]:
                parts = line.split()
                if len(parts) != len(headers): continue
                
                # Parse row
                row_data = {}
                try:
                    # Assume first column is ID if 'id' is in headers, else use index?
                    # LAMMPS usually puts ID first or we look for 'id' column
                    
                    # Convert all to float
                    values = [float(x) for x in parts]
                    
                    # Map headers to values
                    row_data = dict(zip(headers, values))
                    
                    # Determine ID
                    obj_id = int(row_data.get('id', 0))
                    if obj_id == 0:
                        # Fallback if no ID column, use first column as ID?
                        # Or just skip
                        pass
                    
                    data[obj_id] = row_data
                except ValueError:
                    continue
                    
            return data, timestep

        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            return {}, None

    def monitor(self, interval: float = 1.0):
        """
        Continuously monitor files and update plots.
        Blocking call.
        """
        print(f"Monitoring files every {interval}s...")
        try:
            while True:
                self.sync_from_files()
                plt.pause(interval)
        except KeyboardInterrupt:
            print("Monitoring stopped.")
