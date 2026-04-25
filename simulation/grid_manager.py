import os
import math
import shutil
from pathlib import Path
from typing import List, Dict, Any, Tuple
import numpy as np

from .runner import SimulationRunner
from .config import SimulationConfig
from .chain_generator import ChainConfig, write_chain_data

class GridRelaxManager:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def generate_grid_data(self, n_beads_list: List[int], n_states_per_n: int, 
                           spacing: float = 0.5) -> Tuple[Path, Dict[int, Dict[str, Any]]]:
        """
        Creates a single LAMMPS data file with all chains placed in a 3D grid.
        
        Returns:
            Path to the generated data file.
            Dictionary mapping molecule_id to metadata (n_beads, state_index, original_offset).
        """
        all_atoms = []
        all_bonds = []
        all_angles = []
        
        molecule_metadata = {}
        mol_id = 1
        
        total_chains = len(n_beads_list) * n_states_per_n
        grid_dim = math.ceil(total_chains**(1/3))
        
        # We'll use a simple temporary directory for building
        temp_dir = Path("chain_data/grid_build_temp")
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        chain_idx = 0
        for n_beads in n_beads_list:
            # Generate one prototype chain for this N
            # We use a standard vertical linear chain as the starting point
            proto_cfg = ChainConfig(
                beads=n_beads,
                spacing=0.0025,
                mode="linear",
                orientation="vert",
                output_dir=temp_dir
            )
            proto_path = write_chain_data(proto_cfg)
            
            # Read the prototype atoms/bonds
            proto_atoms, proto_bonds, proto_angles = self._parse_data_file(proto_path)
            
            for s_idx in range(n_states_per_n):
                # Calculate grid position
                ix = chain_idx % grid_dim
                iy = (chain_idx // grid_dim) % grid_dim
                iz = chain_idx // (grid_dim * grid_dim)
                
                offset = np.array([ix * spacing, iy * spacing, iz * spacing])
                
                # Shift atoms and assign new IDs
                atom_offset = len(all_atoms)
                for atom in proto_atoms:
                    new_atom = atom.copy()
                    new_atom['id'] += atom_offset
                    new_atom['mol'] = mol_id
                    new_atom['x'] += offset[0]
                    new_atom['y'] += offset[1]
                    new_atom['z'] += offset[2]
                    all_atoms.append(new_atom)
                
                for bond in proto_bonds:
                    new_bond = bond.copy()
                    new_bond['id'] += len(all_bonds)
                    new_bond['a1'] += atom_offset
                    new_bond['a2'] += atom_offset
                    all_bonds.append(new_bond)
                    
                for angle in proto_angles:
                    new_angle = angle.copy()
                    new_angle['id'] += len(all_angles)
                    new_angle['a1'] += atom_offset
                    new_angle['a2'] += atom_offset
                    new_angle['a3'] += atom_offset
                    all_angles.append(new_angle)
                
                molecule_metadata[mol_id] = {
                    "n_beads": n_beads,
                    "state_index": s_idx,
                    "offset": offset
                }
                
                mol_id += 1
                chain_idx += 1
        
        # Write the master grid file
        grid_data_path = Path("chain_data/relaxed_grid_input.data")
        self._write_master_data(grid_data_path, all_atoms, all_bonds, all_angles)
        
        return grid_data_path, molecule_metadata

    def split_grid_results(self, big_data_path: Path, metadata: Dict[int, Dict[str, Any]], 
                           output_base_dir: str):
        """
        Parses the final relaxed grid data file and splits it into individual state files.
        """
        print(f"Splitting {big_data_path} into individual files...")
        all_atoms, all_bonds, all_angles = self._parse_data_file(big_data_path)
        
        # Group by molecule ID
        mol_atoms = {}
        for atom in all_atoms:
            m_id = atom['mol']
            if m_id not in mol_atoms: mol_atoms[m_id] = []
            mol_atoms[m_id].append(atom)
            
        mol_bonds = {}
        for bond in all_bonds:
            # We need to know which molecule this bond belongs to.
            # We'll check the first atom of the bond.
            # (Assuming atoms in a bond belong to the same molecule)
            atom_id = bond['a1']
            # Fast lookup: atom_id to mol_id (could optimize with a map)
            m_id = next(a['mol'] for a in all_atoms if a['id'] == atom_id)
            if m_id not in mol_bonds: mol_bonds[m_id] = []
            mol_bonds[m_id].append(bond)
            
        mol_angles = {}
        for angle in all_angles:
            atom_id = angle['a1']
            m_id = next(a['mol'] for a in all_atoms if a['id'] == atom_id)
            if m_id not in mol_angles: mol_angles[m_id] = []
            mol_angles[m_id].append(angle)
            
        for mol_id, meta in metadata.items():
            n_beads = meta['n_beads']
            s_idx = meta['state_index']
            offset = meta['offset']
            
            target_dir = Path(output_base_dir) / f"N{n_beads}"
            target_dir.mkdir(parents=True, exist_ok=True)
            out_path = target_dir / f"state_{s_idx}.data"
            
            # Re-center atoms for this molecule
            atoms = mol_atoms[mol_id]
            bonds = mol_bonds.get(mol_id, [])
            angles = mol_angles.get(mol_id, [])
            
            # 1-index the atom IDs within the new file
            id_map = {}
            for i, atom in enumerate(sorted(atoms, key=lambda x: x['id'])):
                old_id = atom['id']
                id_map[old_id] = i + 1
                atom['id'] = i + 1
                atom['x'] -= offset[0]
                atom['y'] -= offset[1]
                atom['z'] -= offset[2]
            
            for i, bond in enumerate(bonds):
                bond['id'] = i + 1
                bond['a1'] = id_map[bond['a1']]
                bond['a2'] = id_map[bond['a2']]
                
            for i, angle in enumerate(angles):
                angle['id'] = i + 1
                angle['a1'] = id_map[angle['a1']]
                angle['a2'] = id_map[angle['a2']]
                angle['a3'] = id_map[angle['a3']]
                
            self._write_master_data(out_path, atoms, bonds, angles)
            
        print(f"Successfully split into {len(metadata)} individual files in {output_base_dir}")

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
                if line.startswith("Masses") or line.startswith("Pair Coeffs") or \
                   line.startswith("Bond Coeffs") or line.startswith("Angle Coeffs"):
                    section = "skip"; continue
                
                parts = line.split()
                if section == "atoms" and len(parts) >= 8:
                    # Match chain_generator.py: id type x y z diam dens mol ...
                    atoms.append({
                        'id': int(parts[0]), 
                        'type': int(parts[1]),
                        'x': float(parts[2]), 
                        'y': float(parts[3]), 
                        'z': float(parts[4]),
                        'diam': float(parts[5]), 
                        'dens': float(parts[6]),
                        'mol': int(parts[7])
                    })
                elif section == "bonds" and len(parts) == 4:
                    # id type a1 a2
                    bonds.append({
                        'id': int(parts[0]), 'type': int(parts[1]),
                        'a1': int(parts[2]), 'a2': int(parts[3])
                    })
                elif section == "angles" and len(parts) == 5:
                    # id type a1 a2 a3
                    angles.append({
                        'id': int(parts[0]), 'type': int(parts[1]),
                        'a1': int(parts[2]), 'a2': int(parts[3]), 'a3': int(parts[4])
                    })
                    
        return atoms, bonds, angles

    def _write_master_data(self, path: Path, atoms: List[Dict], bonds: List[Dict], angles: List[Dict]):
        with open(path, 'w') as f:
            f.write("# LAMMPS Grid Relaxation Data File\n\n")
            f.write(f"{len(atoms)} atoms\n")
            f.write(f"{len(bonds)} bonds\n")
            f.write(f"{len(angles)} angles\n\n")
            
            f.write("1 atom types\n")
            f.write("1 bond types\n")
            f.write("1 angle types\n\n")
            
            # Determine box size
            xs = [a['x'] for a in atoms]
            ys = [a['y'] for a in atoms]
            zs = [a['z'] for a in atoms]
            pad = 1.0
            f.write(f"{min(xs)-pad} {max(xs)+pad} xlo xhi\n")
            f.write(f"{min(ys)-pad} {max(ys)+pad} ylo yhi\n")
            f.write(f"{min(zs)-pad} {max(zs)+pad} zlo zhi\n\n")
            
            f.write("Masses\n\n1 1100.0\n\n") 
            
            f.write("Atoms # hybrid sphere molecular\n\n")
            for a in sorted(atoms, key=lambda x: x['id']):
                # id type x y z diam dens mol 0 0 0
                f.write(f"{a['id']} {a['type']} {a['x']} {a['y']} {a['z']} {a['diam']} {a['dens']} {a['mol']} 0 0 0\n")
                
            if bonds:
                f.write("\nBonds\n\n")
                for b in sorted(bonds, key=lambda x: x['id']):
                    f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']}\n")
                    
            if angles:
                f.write("\nAngles\n\n")
                for an in sorted(angles, key=lambda x: x['id']):
                    f.write(f"{an['id']} {an['type']} {an['a1']} {an['a2']} {an['a3']}\n")
