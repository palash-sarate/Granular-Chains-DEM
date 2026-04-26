import os
import math
import shutil
import random
from pathlib import Path
from typing import List, Dict, Any, Tuple
import numpy as np

from .runner import SimulationRunner
from .config import SimulationConfig

class GridHopperManager:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def run_grid_filling(self, n_hoppers: int, n_fill_per_hopper: int, N: int, 
                         hopper_template_data: str = "simulation_geometries/2D_hopper_setup.inc",
                         source_dir: str = None,
                         spacing: float = 1.0,
                         dt: float = 1e-6,
                         relax_steps: int = 500000,
                         seed: int = 12345,
                         output_dir: str = "chain_data/grid_filled",
                         lepton_file: str = "simulation_templates/lepton.inc",
                         dump_file: str = "simulation_templates/quiet_dump.inc",
                         viscosity: float = 0.001,
                         num_procs: int = 1,
                         num_threads: int = 1,
                         use_kokkos: bool = True,
                         mode: str = "2D_stacked",
                         simulation: str = "Grid_Hopper_Filling",
                         template: str = "in.grid_hopper_fill"):
        """
        Main entry point for grid-based batch hopper filling using analytical regions.
        """
        if source_dir is None:
            source_dir = f"chain_data/relaxed/N{N}"
            
        # 1. Prepare molecules
        from .hopper_manager import HopperManager
        h_manager = HopperManager(self.runner)
        mol_dir = "chain_data/molecules_temp"
        inc_file = h_manager.prepare_molecules(source_dir, mol_dir)
        n_templates = len(list(Path(source_dir).glob("*.data")))

        # 2. Setup Grid Geometry (Analytical Regions)
        run_name = f"Grid_Fill_{n_hoppers}H_N{N}_S{seed}"
        job_dir = Path(f"dumping_yard/Grid_Hopper_Filling/{run_name}")
        job_dir.mkdir(parents=True, exist_ok=True)
        
        geometry_inc, metadata = self._generate_replicated_geometry(Path(hopper_template_data), n_hoppers, spacing, job_dir)
        
        # 3. Create Consolidated Insertion File (The Chains)
        insertion_file, z_max = self._generate_grid_insertion_file(
            job_dir, n_hoppers, n_fill_per_hopper, n_templates, seed, N, spacing, metadata, mode=mode
        )

        # 4. Configure & Run Simulation
        config = SimulationConfig(
            template=template,
            simulation=simulation,
            run=run_name,
            data_file=None, # No data file, we use create_box
            lepton_file=lepton_file,
            dump_file=dump_file,
            extra_vars={
                "geometry_inc": geometry_inc,
                "mol_include_file": inc_file,
                "insertion_file": insertion_file,
                "z_max": z_max,
                "relax_steps": relax_steps,
                "dt": dt,
                "seed": seed,
                "N": N,
                "viscosity": viscosity,
                "n_hoppers": n_hoppers,
                "spacing": spacing
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos
        )
        
        print(f"--- Starting Analytical Grid Hopper Filling: {n_hoppers} hoppers ---")
        # Save metadata for later splitting or resuming
        import json
        with open(job_dir / "grid_metadata.json", 'w') as f:
            json.dump({
                "metadata": {str(k): {"offset": v["offset"].tolist()} for k, v in metadata.items()}, 
                "N": N,
                "simulation": simulation
            }, f)

        self.runner.run(config)
        
        # 5. Split Results
        final_grid_data = Path(config.output_dir) / "final_grid.data"
        if final_grid_data.exists():
            self.split_grid_results(final_grid_data, metadata, output_dir)
        else:
            print(f"Error: Final grid data not found at {final_grid_data}")

    def resume_grid_filling(self, restart_path: str, relax_steps: int = 500000,
                           dt: float = 1e-6, output_dir: str = "chain_data/grid_filled",
                           lepton_file: str = "simulation_templates/lepton.inc",
                           dump_file: str = "simulation_templates/quiet_dump.inc",
                           viscosity: float = 0.001,
                           num_procs: int = 1, num_threads: int = 1, use_kokkos: bool = True):
        """
        Resumes a grid filling simulation from a restart file.
        """
        import json
        restart_p = Path(restart_path)
        job_dir = restart_p.parent
        metadata_path = job_dir / "grid_metadata.json"
        
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found at {metadata_path}. Cannot resume/split.")
            
        with open(metadata_path, 'r') as f:
            meta_raw = json.load(f)
            N = meta_raw["N"]
            simulation = meta_raw.get("simulation", "Grid_Hopper_Filling")
            geometry_inc = meta_raw.get("geometry_inc")
            metadata = {int(k): {"offset": np.array(v["offset"])} for k, v in meta_raw["metadata"].items()}

        run_name = f"Resume_{job_dir.name}"
        
        config = SimulationConfig(
            template="in.grid_hopper_fill_resume",
            simulation=simulation,
            run=run_name,
            data_file=str(restart_p), # Passed as restart_path in template
            lepton_file=lepton_file,
            dump_file=dump_file,
            extra_vars={
                "restart_path": str(restart_p),
                "geometry_inc": geometry_inc,
                "relax_steps": relax_steps,
                "dt": dt,
                "viscosity": viscosity
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos
        )
        
        print(f"--- Resuming Grid Hopper Filling from {restart_p.name} ---")
        self.runner.run(config)
        
        final_grid_data = Path(config.output_dir) / "final_grid.data"
        if final_grid_data.exists():
            self.split_grid_results(final_grid_data, metadata, output_dir)

    def _generate_replicated_geometry(self, setup_path: Path, n_hoppers: int, spacing: float, job_dir: Path):
        import re
        def resolve_includes(path: Path):
            lines = []
            with open(path, 'r') as f:
                for line in f:
                    if line.strip().startswith("include"):
                        inc_path = Path(line.split()[1])
                        if not inc_path.is_absolute():
                            inc_path = path.parent / inc_path
                        lines.extend(resolve_includes(inc_path))
                    else:
                        lines.append(line)
            return lines

        original_commands = resolve_includes(setup_path)
        
        # 1. Parse all variables into a Python dictionary for resolution
        vars_dict = {"pi": math.pi}
        for line in original_commands:
            stripped = line.strip()
            if stripped.startswith("variable"):
                parts = stripped.split()
                v_name = parts[1]
                v_expr = " ".join(parts[3:])
                # Remove comments
                v_expr = v_expr.split("#")[0].strip()
                # Simple substitution of already known variables
                for k, v in vars_dict.items():
                    v_expr = v_expr.replace(f"${{{k}}}", str(v))
                
                try:
                    # Handle basic math (sin, cos, tan, etc.)
                    safe_expr = v_expr.replace("sin(", "math.sin(").replace("cos(", "math.cos(").replace("tan(", "math.tan(")
                    vars_dict[v_name] = eval(safe_expr, {"math": math, "__builtins__": None}, vars_dict)
                except:
                    # If it's a complex formula LAMMPS handles, we might fail here, but let's try
                    pass

        # Reserved words that should never be suffixed
        reserved = {"side", "out", "in", "units", "move", "rotate", "open", "all"}
        # Extract all original region names to ensure we only suffix valid references
        original_region_names = {l.strip().split()[1] for l in original_commands if l.strip().startswith("region")}
        
        grid_dim = math.ceil(n_hoppers**(1/3))
        replicated_lines = []
        
        # 2. Write global variables/comments once
        replicated_lines.append("# --- Global Variables (Resolved in Python) ---\n")
        for k, v in vars_dict.items():
            if k != "pi":
                replicated_lines.append(f"variable {k} equal {v}\n")
        
        metadata = {}
        for i in range(n_hoppers):
            ix, iy, iz = i % grid_dim, (i // grid_dim) % grid_dim, i // (grid_dim * grid_dim)
            offset = np.array([ix * spacing, iy * spacing, iz * spacing])
            metadata[i] = {"offset": offset}
            
            replicated_lines.append(f"\n# --- Replicated Hopper {i} at {offset} ---\n")
            
            for line in original_commands:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith("variable"):
                    continue
                
                # Replicate Regions
                if stripped.startswith("region"):
                    parts = stripped.split()
                    name = parts[1]
                    style = parts[2]
                    rest = parts[3:]
                    new_name = f"{name}_{i}"
                    
                    def resolve_val(val, off):
                        # Try to resolve val using vars_dict
                        expr = val
                        for k, v in vars_dict.items():
                            expr = expr.replace(f"${{{k}}}", str(v))
                        try:
                            # Handle leading minus
                            clean_expr = expr.replace("math.", "") # in case it was already processed
                            return float(eval(clean_expr, {"math": math, "__builtins__": None})) + off
                        except:
                            return f"{val}+{off}" # Fallback if we can't resolve

                    if style == "block":
                        coords = [str(resolve_val(c, offset[j//2])) for j, c in enumerate(rest[:6])]
                        final_rest = [f"{p}_{i}" if p in original_region_names else p for p in rest[6:]]
                        new_line = f"region {new_name} block {' '.join(coords)} {' '.join(final_rest)}\n"
                    elif style == "plane":
                        coords = [str(resolve_val(c, offset[j])) for j, c in enumerate(rest[:3])]
                        final_rest = [f"{p}_{i}" if p in original_region_names else p for p in rest[3:]]
                        new_line = f"region {new_name} plane {' '.join(coords)} {' '.join(final_rest)}\n"
                    elif style in ["intersect", "union"]:
                        num = rest[0]
                        regs = [f"{r}_{i}" if r in original_region_names else r for r in rest[1:]]
                        new_line = f"region {new_name} {style} {num} {' '.join(regs)}\n"
                    else:
                        new_line = f"region {new_name} {style} {' '.join([f'{p}_{i}' if p in original_region_names else p for p in rest])}\n"
                    replicated_lines.append(new_line)
                
                # Replicate Fixes
                elif stripped.startswith("fix"):
                    parts = stripped.split()
                    f_id = parts[1]
                    f_group = parts[2]
                    f_style = parts[3]
                    new_id = f"{f_id}_{i}"
                    
                    new_parts = []
                    for p in parts:
                        if p in original_region_names:
                            new_parts.append(f"{p}_{i}")
                        else:
                            new_parts.append(p)
                    new_parts[1] = new_id
                    replicated_lines.append(f"{' '.join(new_parts)}\n")
                    
        out_path = job_dir / "replicated_geometry.inc"
        with open(out_path, 'w') as f:
            f.writelines(replicated_lines)
        return str(out_path).replace("\\", "/"), metadata

    def _generate_grid_insertion_file(self, outdir: Path, n_hoppers: int, n_fill: int, 
                                      n_templates: int, seed: int, N: int, 
                                      spacing: float, metadata: Dict[int, Any],
                                      mode: str = "2D_stacked"):
        rng = random.Random(seed)
        insertion_lines = []
        
        # Internal pouring grid relative to hopper center
        # y_half from 2D_hopper.inc is 0.155
        hopper_y_half = 0.155
        
        # If chains are oriented horizontally or can tumble, 
        # we need a buffer of at least N*bead_radius/2.
        # Even for vertical stacking, we want a small safety margin.
        bead_diam = 0.003
        safe_buffer = max(0.005, (N * bead_diam / 2.0) + 0.002)
        
        y_max = hopper_y_half - safe_buffer
        y_min = -hopper_y_half + safe_buffer
        
        # Ensure we didn't shrink the region to nothing for very long chains
        if y_max <= y_min:
            y_max = 0.001
            y_min = -0.001
            
        z_start = 0.51
        y_width = y_max - y_min
        
        # Exclusion distance for 2D stacking
        # Horizontal (Y) gap: beads are ~2mm, so 2.5mm spacing is safe
        # Vertical (Z) gap: length of chain + small buffer
        dy_gap = 0.003 
        dz_gap = (N * 0.0025) + 0.005 

        z_max_global = z_start
        total_inserted = 0
        
        for h_idx in range(n_hoppers):
            offset = metadata[h_idx]['offset']
            count = 0
            
            if mode == "2D_stacked":
                ny = max(1, int(y_width / dy_gap))
                nz = math.ceil(n_fill / ny)
                
                for iz in range(nz):
                    if count >= n_fill: break
                    pz = z_start + iz * dz_gap
                    z_max_global = max(z_max_global, pz + offset[2])
                    
                    for iy in range(ny):
                        if count >= n_fill: break
                        px = offset[0] # Fixed X at center of grid cell
                        py = y_min + (iy + 0.5) * (y_width/ny) + offset[1]
                        
                        mol_id = rng.randint(1, n_templates)
                        insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz+offset[2]:.6f} mol m{mol_id} 12345 rotate 0.0 0.0 0.0 1.0")
                        count += 1
                        total_inserted += 1
            else:
                # Original 3D insertion logic
                x_min, x_max = -0.08, 0.08
                x_width = x_max - x_min
                ex_diam = 3 * N * 0.001
                nx = max(1, int(x_width / ex_diam))
                ny = max(1, int(y_width / ex_diam))
                grid_2d = nx * ny
                nz = math.ceil(n_fill / grid_2d)
                dz = ex_diam
                
                for iz in range(nz):
                    if count >= n_fill: break
                    ox, oy = rng.uniform(0, x_width/nx), rng.uniform(0, y_width/ny)
                    pz = z_start + iz * dz
                    z_max_global = max(z_max_global, pz + offset[2])
                    for ix in range(nx):
                        if count >= n_fill: break
                        for iy in range(ny):
                            if count >= n_fill: break
                            px = x_min + (ix + 0.5) * (x_width/nx) + ox + offset[0]
                            py = y_min + (iy + 0.5) * (y_width/ny) + oy + offset[1]
                            
                            mol_id = rng.randint(1, n_templates)
                            mol_seed = seed + total_inserted
                            insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz+offset[2]:.6f} mol m{mol_id} {mol_seed}")
                            count += 1
                            total_inserted += 1
                        
        insertion_path = outdir / "insertions.inc"
        with open(insertion_path, 'w') as f:
            f.write("\n".join(insertion_lines))
        return str(insertion_path).replace("\\", "/"), z_max_global

    def split_grid_results(self, big_data_path: Path, metadata: Dict[int, Any], output_base_dir: str):
        # [Splitting logic similar to GridRelaxManager but aware of hopper boundaries]
        print(f"Splitting grid filled results into {len(metadata)} hoppers...")
        all_atoms, all_bonds, all_angles = self._parse_data_file(big_data_path)
        
        for h_idx, meta in metadata.items():
            offset = meta['offset']
            # Find atoms belonging to this hopper (within a bounding box around the offset)
            # A hopper is roughly 0.5m wide, so we'll use a +/- 0.45m filter
            h_atoms = [a for a in all_atoms if abs(a['x'] - offset[0]) < 0.45 and abs(a['y'] - offset[1]) < 0.45]
            
            target_dir = Path(output_base_dir) / f"state_{h_idx}"
            target_dir.mkdir(parents=True, exist_ok=True)
            out_path = target_dir / "final_hopper.data"
            
            # Re-center and re-ID
            id_map = {}
            for i, a in enumerate(sorted(h_atoms, key=lambda x: x['id'])):
                id_map[a['id']] = i + 1
                a['id'] = i + 1
                a['x'] -= offset[0]
                a['y'] -= offset[1]
                a['z'] -= offset[2]
                if a['mol'] >= 1000000: a['mol'] = 0 # Mark as hopper wall if needed
            
            # Filter and re-ID bonds/angles
            h_atom_ids = set(id_map.keys())
            bonds = [b for b in all_bonds if b['a1'] in h_atom_ids]
            for i, b in enumerate(bonds):
                b['id'] = i + 1
                b['a1'] = id_map[b['a1']]
                b['a2'] = id_map[b['a2']]
            
            angles = [an for an in all_angles if an['a1'] in h_atom_ids]
            for i, an in enumerate(angles):
                an['id'] = i + 1
                an['a1'] = id_map[an['a1']]
                an['a2'] = id_map[an['a2']]
                an['a3'] = id_map[an['a3']]
                
            self._write_data(out_path, h_atoms, bonds, angles)
        print("Splitting complete.")

    def _parse_data_file(self, path: Path):
        atoms, bonds, angles = [], [], []
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
                    atoms.append({'id': int(parts[0]), 'type': int(parts[1]), 'x': float(parts[2]), 'y': float(parts[3]), 'z': float(parts[4]), 'diam': float(parts[5]), 'dens': float(parts[6]), 'mol': int(parts[7])})
                elif section == "bonds" and len(parts) == 4:
                    bonds.append({'id': int(parts[0]), 'type': int(parts[1]), 'a1': int(parts[2]), 'a2': int(parts[3])})
                elif section == "angles" and len(parts) == 5:
                    angles.append({'id': int(parts[0]), 'type': int(parts[1]), 'a1': int(parts[2]), 'a2': int(parts[3]), 'a3': int(parts[4])})
        return atoms, bonds, angles

    def _write_data(self, path: Path, atoms: List[Dict], bonds: List[Dict], angles: List[Dict]):
        with open(path, 'w') as f:
            f.write("# LAMMPS Data\n\n")
            f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n\n")
            f.write(f"2 atom types\n1 bond types\n1 angle types\n\n")
            xs, ys, zs = [a['x'] for a in atoms], [a['y'] for a in atoms], [a['z'] for a in atoms]
            pad = 0.5
            f.write(f"{min(xs)-pad} {max(xs)+pad} xlo xhi\n{min(ys)-pad} {max(ys)+pad} ylo yhi\n{min(zs)-pad} {max(zs)+pad} zlo zhi\n\n")
            f.write("Masses\n\n1 1100.0\n2 1100.0\n\nAtoms # hybrid sphere molecular\n\n")
            for a in sorted(atoms, key=lambda x: x['id']):
                f.write(f"{a['id']} {a['type']} {a['x']} {a['y']} {a['z']} {a['diam']} {a['dens']} {a['mol']} 0 0 0\n")
            if bonds:
                f.write("\nBonds\n\n")
                for b in bonds: f.write(f"{b['id']} {b['type']} {b['a1']} {b['a2']}\n")
            if angles:
                f.write("\nAngles\n\n")
                for an in angles: f.write(f"{an['id']} {an['type']} {an['a1']} {an['a2']} {an['a3']}\n")
