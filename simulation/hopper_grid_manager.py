import os
import math
import shutil
from pathlib import Path
from typing import List, Dict, Any, Tuple
import numpy as np

from .runner import SimulationRunner
from .config import SimulationConfig

class HopperGridManager:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def create_hopper_grid_data(self, hopper_template_path: Path, n_hoppers: int, 
                                 spacing: float = 2.0) -> Tuple[Path, Dict[int, Any]]:
        """
        Clones a base hopper geometry into a 3D grid.
        
        Returns:
            Path to the new grid data file.
            Metadata about hopper offsets.
        """
        # 1. Parse the base hopper
        h_atoms, h_bonds, h_angles = self._parse_data_file(hopper_template_path)
        
        grid_dim = math.ceil(n_hoppers**(1/3))
        all_atoms = []
        hopper_metadata = {}
        
        for i in range(n_hoppers):
            ix = i % grid_dim
            iy = (i // grid_dim) % grid_dim
            iz = i // (grid_dim * grid_dim)
            
            offset = np.array([ix * spacing, iy * spacing, iz * spacing])
            
            # Shift atoms and assign to a new molecule ID representing the hopper number
            # (Though hopper atoms are usually fixed/type 2, we still offset them)
            atom_offset = len(all_atoms)
            for atom in h_atoms:
                new_atom = atom.copy()
                new_atom['id'] += atom_offset
                new_atom['x'] += offset[0]
                new_atom['y'] += offset[1]
                new_atom['z'] += offset[2]
                # We can use molecule ID to distinguish which hopper wall it belongs to
                new_atom['mol'] = i + 1 
                all_atoms.append(new_atom)
                
            hopper_metadata[i + 1] = {
                "offset": offset,
                "index": i
            }
            
        # Write the grid data
        grid_path = Path("chain_data/hopper_grid_layout.data")
        self._write_data(grid_path, all_atoms, h_bonds, h_angles) # Hoppers usually don't have bonds/angles
        
        return grid_path, hopper_metadata

    def _parse_data_file(self, path: Path):
        atoms = []
        bonds = []
        angles = []
        with open(path, 'r') as f:
            section = None
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                if line.startswith("Atoms"): section = "atoms"; continue
                if line.startswith("Bonds"): section = "bonds"; continue
                if line.startswith("Angles"): section = "angles"; continue
                
                parts = line.split()
                if section == "atoms" and len(parts) >= 8:
                    atoms.append({
                        'id': int(parts[0]), 'type': int(parts[1]),
                        'x': float(parts[2]), 'y': float(parts[3]), 'z': float(parts[4]),
                        'diam': float(parts[5]), 'dens': float(parts[6]), 'mol': int(parts[7])
                    })
        return atoms, bonds, angles

    def _write_data(self, path: Path, atoms: List[Dict], bonds: List[Dict], angles: List[Dict]):
        with open(path, 'w') as f:
            f.write("# LAMMPS Hopper Grid Data\n\n")
            f.write(f"{len(atoms)} atoms\n\n")
            
            # Determine types (hoppers usually type 2, chains type 1)
            # For the base layout, it might just be the hopper atoms
            types = set(a['type'] for a in atoms)
            f.write(f"{max(types) if types else 2} atom types\n\n")
            
            xs, ys, zs = [a['x'] for a in atoms], [a['y'] for a in atoms], [a['z'] for a in atoms]
            pad = 2.0
            f.write(f"{min(xs)-pad} {max(xs)+pad} xlo xhi\n")
            f.write(f"{min(ys)-pad} {max(ys)+pad} ylo yhi\n")
            f.write(f"{min(zs)-pad} {max(zs)+pad} zlo zhi\n\n")
            
            f.write("Masses\n\n1 1100.0\n2 1100.0\n\n")
            
            f.write("Atoms # hybrid sphere molecular\n\n")
            for a in sorted(atoms, key=lambda x: x['id']):
                f.write(f"{a['id']} {a['type']} {a['x']} {a['y']} {a['z']} {a['diam']} {a['dens']} {a['mol']} 0 0 0\n")
