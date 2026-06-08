import os
import math
import shutil
import random
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import numpy as np
import time
import json

from .runner import SimulationRunner
from .config import SimulationConfig

class GridHopperManager:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def _ensure_local_directory(self, target_path: str):
        """Helper to ensure a directory is present locally, restoring it from Drive if synced."""
        from pathlib import Path
        import os
        
        path_p = Path(target_path)
        path_str = str(path_p.absolute()).replace("\\", "/")
        
        # 1. Resolve to the main run directory by finding the segment that matches the run name pattern
        if "/restart" in path_str:
            dir_str = path_str.split("/restart")[0]
            dir_p = Path(dir_str)
        elif "/results" in path_str:
            dir_str = path_str.split("/results")[0]
            dir_p = Path(dir_str)
        else:
            if path_p.suffix:
                dir_p = path_p.parent
            else:
                dir_p = path_p
                
        # 2. Check if the local directory actually exists and is complete
        metadata_exists = (dir_p / "metadata.json").exists() or (dir_p / "grid_metadata.json").exists()
        if not dir_p.exists() or not metadata_exists:
            # 3. Find the exact matching lineage key by matching relative path
            target_rel = str(dir_p).replace("\\", "/")
            # Strip common root prefixes to isolate the relative path
            target_rel = target_rel.replace("/home/guest/palash/Granular-Chains-DEM/", "")
            target_rel = target_rel.replace("/Data/palash_data/", "")
            target_rel = target_rel.lstrip("/")
            
            try:
                from Pulse.pulse_core import PBSManager, SyncManager
                lineage = PBSManager.load_lineage()
                
                matching_key = None
                for k, info in lineage.items():
                    k_rel = info.get('rel_path', '').replace("\\", "/").lstrip("/")
                    if target_rel == k_rel or k == str(dir_p.absolute()).replace("\\", "/"):
                        matching_key = k
                        break
                
                if matching_key:
                    print(f"[RECOVERY] Local directory for {target_rel} is missing or incomplete. Restoring from cloud...")
                    success, msg = SyncManager.restore_run(matching_key)
                    if success:
                        print(f"[RECOVERY] {msg}")
                    else:
                        print(f"[RECOVERY WARNING] Could not restore: {msg}")
                else:
                    print(f"[RECOVERY WARNING] Could not find a lineage match for {target_rel}. Cannot restore from cloud.")
            except Exception as e:
                print(f"[RECOVERY ERROR] Failed during cloud restoration attempt: {e}")



    def _resolve_includes(self, path: Path):
        lines = []
        with open(path, 'r') as f:
            for line in f:
                if line.strip().startswith("include"):
                    inc_path = Path(line.split()[1])
                    if not inc_path.is_absolute():
                        inc_path = path.parent / inc_path
                    lines.extend(self._resolve_includes(inc_path))
                else:
                    lines.append(line)
        return lines

    def _generate_resume_path(self, job_dir: Path, seed: int, simulation: Optional[str] = None) -> Tuple[Path, str]:
        """
        Generates a standardized resume path: OriginalName_S<seed1>_S<seed2>
        """
        original_name = job_dir.name
        # Clean up any legacy _ResN suffixes if present
        import re
        base_name = re.sub(r'_Res\d+', '', original_name)
        run_name = f"{base_name}_S{seed}"
        
        # Determine the simulation type from the job_dir path if not provided
        # dumping_yard/<SimulationType>/<RunName>
        if simulation is None:
            simulation = job_dir.parent.name
        new_path = Path(f"dumping_yard/{simulation}/{run_name}")
        return new_path, run_name

    def run_grid_filling(self, n_hoppers: int, n_fill_per_hopper: Any, N: Any, 
                         hopper_template_data: str = "simulation_geometries/2D_hopper_setup.inc",
                         source_dir: str = None,
                         spacing: float = 1.0,
                         dt: float = 1e-6,
                         relax_steps: int = 500000,
                         seed: int = 12345,
                         lepton_file: str = "simulation_templates/lepton.inc",
                         dump_file: str = "simulation_templates/default_dump.inc",
                         viscosity: float = 0.000001,
                         num_procs: int = 1,
                         num_threads: int = 1,
                         use_kokkos: bool = True,
                         mode: str = "2D_stacked",
                         simulation: str = "Hopper_Fill",
                         template: str = "in.grid_hopper_fill",
                         generate_vtk: bool = True,
                         geometry_vars: Dict[str, Any] = None):
        t_start = time.time()
        print(f"[INFO] Grid filling simulation started. VTK Generation: {generate_vtk}")
        # 1. Normalize N and n_fill into lists
        if isinstance(N, (int, float)):
            N_list = [int(N)] * n_hoppers
        else:
            N_list = list(N)
            if len(N_list) < n_hoppers:
                N_list = (N_list * (n_hoppers // len(N_list) + 1))[:n_hoppers]
        
        if isinstance(n_fill_per_hopper, (int, float)):
            n_fill_list = [int(n_fill_per_hopper)] * n_hoppers
        else:
            n_fill_list = list(n_fill_per_hopper)
            if len(n_fill_list) < n_hoppers:
                n_fill_list = (n_fill_list * (n_hoppers // len(n_fill_list) + 1))[:n_hoppers]
        
        # 1b. Normalize geometry variables
        normalized_geo_vars = [{} for _ in range(n_hoppers)]
        if geometry_vars:
            for var_name, values in geometry_vars.items():
                if isinstance(values, (list, tuple)):
                    # Cycle values if shorter than n_hoppers
                    expanded = (list(values) * (n_hoppers // len(values) + 1))[:n_hoppers]
                    for i in range(n_hoppers):
                        normalized_geo_vars[i][var_name] = expanded[i]
                else:
                    # Single value applied to all
                    for i in range(n_hoppers):
                        normalized_geo_vars[i][var_name] = values

        unique_Ns = sorted(list(set(N_list)))
        cols = math.ceil(math.sqrt(n_hoppers))

        # 2. Setup Run Name and Job Directory
        mixed_tag = "_MixedN" if len(unique_Ns) > 1 else ""
        run_name = f"Grid_Fill_{n_hoppers}H{mixed_tag}_S{seed}"
        job_dir = Path(f"dumping_yard/{simulation}/{run_name}")
        job_dir.mkdir(parents=True, exist_ok=True)

        # 3. Prepare molecules from multiple N sources in a job-specific directory
        from .hopper_manager import HopperManager
        from .molecule_converter import convert_data_to_molecule
        mol_dir = job_dir / "molecules"
        mol_dir.mkdir(parents=True, exist_ok=True)
        
        relaxed_sources = {} # N -> source_path
        mol_ranges = {} # N -> { 'start': int, 'count': int }
        mol_bboxes = {} # mol_id -> bbox
        combined_inc_lines = []
        current_id_offset = 0
        
        for n_val in unique_Ns:
            if source_dir:
                n_source = Path(source_dir) / f"N{n_val}"
            else:
                n_source = Path(f"chain_data/relaxed/N{n_val}")
            
            relaxed_sources[n_val] = str(n_source).replace("\\", "/")
            data_files = list(n_source.glob("*.data"))
            if not data_files:
                print(f"Warning: No templates found for N={n_val} in {n_source}")
                continue
            
            mol_ranges[n_val] = { 'start': current_id_offset + 1, 'count': len(data_files) }
            
            for i, data_file in enumerate(data_files):
                mol_id = current_id_offset + i + 1
                mol_filename = f"mol_{mol_id}.mol"
                output_mol = mol_dir / mol_filename
                bbox = convert_data_to_molecule(str(data_file), str(output_mol))
                mol_bboxes[mol_id] = bbox
                # Use path relative to project root for LAMMPS
                mol_rel_path = str(output_mol).replace("\\", "/")
                combined_inc_lines.append(f"molecule m{mol_id} {mol_rel_path}")
            
            current_id_offset += len(data_files)
            
        inc_file = mol_dir / "molecules.inc"
        with open(inc_file, 'w') as f:
            f.write("\n".join(combined_inc_lines))
        
        # 4. Setup Grid Geometry (Analytical Regions)
        geometry_inc, envelope = self._generate_replicated_geometry(
            Path(hopper_template_data), n_hoppers, spacing, job_dir, normalized_geo_vars
        )
        
        # geometry_inc, envelope = self._generate_replicated_geometry(
        #     Path(hopper_template_data), n_hoppers, spacing, job_dir, normalized_geo_vars
        # )

        # 3.2 Optional: Generate VTK meshes for UI visualization
        if generate_vtk:
            print("--- Generating VTK Mesh for UI Visualization ---")
            try:
                from analysis.geometry_extractor import GeometryExtractor
                extractor = GeometryExtractor(self.runner.lammps_exe)
                
                vtk_bounds = [envelope['total_bounds']['x'][0], envelope['total_bounds']['x'][1],
                              envelope['total_bounds']['y'][0], envelope['total_bounds']['y'][1],
                              envelope['total_bounds']['z'][0], envelope['total_bounds']['z'][1]]
                
                extractor.extract(
                    inc_file=geometry_inc,
                    outdir=job_dir / "Geometry_vtk",
                    auto_vis=True,
                    combined=False,
                    bounds=vtk_bounds,
                    spacing=0.002, # High resolution to resolve small features (e.g. 2.5mm gaps)
                    num_procs=num_procs,
                    num_threads=num_threads,
                    use_kokkos=False # Permanent -no-kokkos for extraction as it often conflicts with OVITO/sampling
                )
            except Exception as e:
                print(f"[WARNING] Could not generate VTK mesh: {e}")
                print("You can generate it manually later using analysis/geometry_extractor.py")
        
        # 4. Create Consolidated Insertion File (The Chains)
        z_start = 0.51
        try:
            original_commands = self._resolve_includes(Path(hopper_template_data))
            global_overrides = normalized_geo_vars[0] if normalized_geo_vars else {}
            vars_dict = self._evaluate_variables(original_commands, global_overrides)
            if "converging_sec_ht" in vars_dict:
                z_start = float(vars_dict["converging_sec_ht"])
                print(f"[INFO] Using dynamic z_start={z_start:.4f} from geometry variable converging_sec_ht")
        except Exception as e:
            print(f"[WARNING] Could not parse geometry for converging_sec_ht: {e}. Falling back to z_start=0.51")

        insertion_file, z_max = self._generate_grid_insertion_file(
            job_dir, n_hoppers, n_fill_list, mol_ranges, seed, N_list, spacing, envelope['hoppers'], 
            mol_bboxes=mol_bboxes, mode=mode, z_start=z_start
        )

        setup_duration = time.time() - t_start

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
                "N": N_list[0],
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
        metadata = {str(k): {"offset": v["offset"].tolist(), "N": N_list[k], "n_fill": n_fill_list[k], "geometry_vars": normalized_geo_vars[k]} for k, v in envelope['hoppers'].items()}
        
        import json
        with open(job_dir / "metadata.json", 'w') as f:
            json.dump({
                "metadata": metadata, 
                "total_bounds": envelope['total_bounds'],
                "N": N_list,
                "n_fill": n_fill_list,
                "spacing": spacing,
                "geometry_vars": normalized_geo_vars,
                "geometry_inc": geometry_inc,
                "simulation": simulation,
                "relaxed_sources": relaxed_sources,
                "n_hoppers": n_hoppers,
                "cols": cols,
                "hopper_template_data": str(hopper_template_data)
            }, f)

        self.runner.run(config, clean_dir=False)
        
        # 5. Split Results
        for k_str in metadata:
            k = int(k_str)
            metadata[k_str]['N'] = N_list[k]
            metadata[k_str]['n_fill'] = n_fill_list[k]
            metadata[k_str]['source_dir'] = relaxed_sources.get(N_list[k])

        final_grid_data = Path(config.output_dir) / "final_grid.data"
        if final_grid_data.exists():
            split_dir = Path(config.output_dir) / "split_states"
            self.split_grid_results(final_grid_data, metadata, str(split_dir), spacing=spacing)
        else:
            print(f"Error: Final grid data not found at {final_grid_data}")

        return setup_duration

    def resume_grid_filling(self, restart_path: str, relax_steps: int = 500000,
                           dt: float = 1e-6,
                           lepton_file: str = "simulation_templates/lepton.inc",
                           dump_file: str = "simulation_templates/default_dump.inc",
                           viscosity: float = 0.000001,
                           num_procs: int = 1, num_threads: int = 1, use_kokkos: bool = True,
                           template: str = "in.grid_hopper_fill_resume",
                           simulation: str = None,
                           seed: int = 42,
                           inplace: bool = False):
        """
        Resumes a grid filling simulation from a restart file.
        """
        import json
        self._ensure_local_directory(restart_path)
        restart_p = Path(restart_path)
        if restart_p.is_dir():
            prev_job_dir = restart_p
            # Find latest restart in directory
            r_file = restart_p / "restart" / "restart.final.bin"
            if not r_file.exists():
                restarts = list((restart_p / "restart").glob("restart.*.bin"))
                if restarts:
                    # Sort by step number in restart.STEP.bin
                    try:
                        r_file = max(restarts, key=lambda p: int(p.stem.split('.')[-1]) if '.' in p.stem else 0)
                    except:
                        r_file = restarts[-1]
                else:
                    raise FileNotFoundError(f"No restart files found in {restart_p}/restart")
            restart_p = r_file
        else:
            prev_job_dir = restart_p.parent
            
        metadata_path = prev_job_dir / "metadata.json"
        if not metadata_path.exists():
            # Try one level up if we were inside a 'restart' or 'results' subfolder
            metadata_path = prev_job_dir.parent / "metadata.json"
            if metadata_path.exists():
                prev_job_dir = prev_job_dir.parent
        
        # Backward compatibility check for legacy naming
        if not metadata_path.exists():
            legacy_names = ["grid_metadata.json", "sim_metadata.json"]
            for name in legacy_names:
                p = prev_job_dir / name
                if p.exists():
                    metadata_path = p
                    break
                    
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found at {metadata_path}. Cannot resume/split.")
            
        with open(metadata_path, 'r') as f:
            meta_raw = json.load(f)
            if simulation is None:
                simulation = meta_raw.get("simulation", "Grid_Hopper_Filling")
            geometry_inc = meta_raw.get("geometry_inc")
            spacing = meta_raw.get("spacing", 2.0)
            
            # Reconstruction parameters for self-contained geometry
            n_hoppers = meta_raw.get("n_hoppers", len(meta_raw["metadata"]))
            cols = meta_raw.get("cols", 1)
            geo_vars = meta_raw.get("geometry_vars", [])
            hopper_template = meta_raw.get("hopper_template_data", "simulation_geometries/2D_hopper_with_orifice_cover.inc")

            metadata = {int(k): {
                "offset": np.array(v["offset"]),
                "N": v.get("N"),
                "n_fill": v.get("n_fill"),
                "geometry_vars": v.get("geometry_vars")
            } for k, v in meta_raw["metadata"].items()}

        # Naming Strategy: Inplace vs New Folder
        if inplace:
            new_job_dir = prev_job_dir
            run_name = prev_job_dir.name
        else:
            run_name = f"{prev_job_dir.name}_S{seed}"
            new_job_dir = Path(f"dumping_yard/{simulation}/{run_name}")
            new_job_dir.mkdir(parents=True, exist_ok=True)
        
        # Regeneration of geometry to ensure self-contained jobs (fixes missing parent files after sync)
        local_geometry_inc, _ = self._generate_replicated_geometry(
            Path(hopper_template), n_hoppers, spacing, new_job_dir, geo_vars
        )
        
        # Metadata Inheritance: Copy and Update with Lineage
        new_meta = meta_raw.copy()
        if not inplace:
            new_meta["source_dir"] = str(prev_job_dir.absolute()).replace("\\", "/")
            new_meta["run_name"] = run_name
        
        # Track active seeds for dashboard color logic
        if "active_seeds" not in new_meta: new_meta["active_seeds"] = []
        if str(seed) not in new_meta["active_seeds"]: new_meta["active_seeds"].append(str(seed))
        
        # Reset sync status if resuming in-place
        if inplace:
            self._reset_sync_status(str(new_job_dir.absolute()))

        with open(new_job_dir / "grid_metadata.json", 'w') as f:
            json.dump(new_meta, f, indent=4)
        
        config = SimulationConfig(
            template="in.grid_hopper_fill_resume",
            simulation=simulation,
            run=run_name,
            data_file=str(restart_p), # Passed as restart_path in template
            lepton_file=lepton_file,
            dump_file=dump_file,
            outdir_override=str(new_job_dir),
            extra_vars={
                "restart_path": str(restart_p).replace("\\", "/"),
                "geometry_inc": local_geometry_inc.replace("\\", "/"),
                "relax_steps": relax_steps,
                "dt": dt,
                "viscosity": viscosity
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos
        )
        
        print(f"--- Resuming Grid Hopper Filling from {restart_p.name} ---")
        self.runner.run(config, clean_dir=False)
        
        final_grid_data = Path(config.output_dir) / "final_grid.data"
        if final_grid_data.exists():
            split_dir = Path(config.output_dir) / "split_states"
            self.split_grid_results(final_grid_data, metadata, str(split_dir), spacing=spacing)

    def run_grid_flow(self, source_dir: str, run_steps: int = 1000000,
                      freq: Any = 10.0, amp: Any = 0.01, osc_dir: str = 'z',
                      dt: float = 1e-6,
                      lepton_file: str = "simulation_templates/lepton.inc",
                      dump_file: str = "simulation_templates/default_dump.inc",
                      viscosity: float = 0.000001,
                      num_procs: int = 1, num_threads: int = 1, use_kokkos: bool = True,
                      simulation: str = "Grid_Hopper_Flow",
                      template: str = "in.grid_hopper_flow",
                      seed: int = 12345,
                      inplace: bool = False):
        """
        Takes a filled grid state and starts the flow simulation (opens orifices + oscillation).
        """
        import json
        self._ensure_local_directory(source_dir)
        source_p = Path(source_dir)
        metadata_path = source_p / "metadata.json"
        if not metadata_path.exists():
            metadata_path = source_p / "grid_metadata.json" # Legacy fallback
            
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found at {metadata_path}. Cannot start flow.")

        with open(metadata_path, 'r') as f:
            meta_raw = json.load(f)
            n_hoppers = len(meta_raw["metadata"])
            spacing = meta_raw.get("spacing", 2.0)
            geometry_vars = meta_raw.get("geometry_vars", [])
            # Reconstruct metadata for splitting later
            metadata_original = {int(k): {
                "offset": np.array(v["offset"]),
                "N": v.get("N"),
                "n_fill": v.get("n_fill"),
                "geometry_vars": v.get("geometry_vars")
            } for k, v in meta_raw["metadata"].items()}

        # 1. Normalize freq and amp
        if isinstance(freq, (int, float)):
            freq_list = [float(freq)] * n_hoppers
        else:
            freq_list = (list(freq) * (n_hoppers // len(freq) + 1))[:n_hoppers]
            
        if isinstance(amp, (int, float)):
            amp_list = [float(amp)] * n_hoppers
        else:
            amp_list = (list(amp) * (n_hoppers // len(amp) + 1))[:n_hoppers]

        # 2. Setup Flow Run Directory
        if inplace:
            job_dir = source_p
            run_name = source_p.name
        else:
            is_mixed = "_MixedN" in source_p.name
            run_name = f"Grid_Flow_{n_hoppers}H{'_MixedN' if is_mixed else ''}_S{seed}"
            job_dir = Path(f"dumping_yard/{simulation}/{run_name}")
            job_dir.mkdir(parents=True, exist_ok=True)

        # 3. Generate Flow Geometry (No lids, includes oscillation variables)
        hopper_template = "simulation_geometries/2D_hopper_with_orifice_cover.inc"
        
        osc_params = [{"freq": freq_list[i], "amp": amp_list[i], "dir": osc_dir} for i in range(n_hoppers)]
        
        geometry_flow_inc, _ = self._generate_replicated_geometry(
            Path(hopper_template), n_hoppers, spacing, job_dir, geometry_vars,
            exclude_regions=["orifice_cover"],
            oscillation_params=osc_params,
            inc_name="replicated_geometry_flow.inc"
        )

        # 4. Save New Metadata
        new_meta = meta_raw.copy()
        new_meta["freq"] = freq_list
        new_meta["amp"] = amp_list
        new_meta["osc_dir"] = osc_dir
        if not inplace:
            new_meta["source_dir"] = str(source_p).replace("\\", "/")
            new_meta["source_run"] = source_p.name
            new_meta["run_name"] = run_name
        new_meta["simulation"] = simulation
        new_meta["geometry_inc"] = geometry_flow_inc
        new_meta["hopper_template_data"] = str(hopper_template)
        
        # Track active seeds
        if "active_seeds" not in new_meta: new_meta["active_seeds"] = []
        if str(seed) not in new_meta["active_seeds"]: new_meta["active_seeds"].append(str(seed))

        if inplace:
            self._reset_sync_status(str(job_dir.absolute()))

        with open(job_dir / "grid_metadata.json", 'w') as f:
            json.dump(new_meta, f, indent=4)

        # 6. Configure & Run
        restart_path = source_p / "restart" / "restart.final.bin"
        if not restart_path.exists():
            restarts = list((source_p / "restart").glob("restart.*.bin"))
            if restarts:
                restart_path = sorted(restarts, key=os.path.getmtime)[-1]
            else:
                raise FileNotFoundError(f"No restart file found in {source_p / 'restart'}")

        config = SimulationConfig(
            template=template,
            simulation=simulation,
            run=run_name,
            data_file=None,
            lepton_file=lepton_file,
            dump_file=dump_file,
            outdir_override=str(job_dir),
            extra_vars={
                "restart_path": str(restart_path).replace("\\", "/"),
                "geometry_inc": geometry_flow_inc.replace("\\", "/"),
                "run_steps": run_steps,
                "dt": dt,
                "viscosity": viscosity,
                "seed": seed
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos
        )

        print(f"--- Starting Grid Hopper Flow: {n_hoppers} hoppers ---")
        self.runner.run(config, clean_dir=False)

        # 7. Split Results
        final_grid_data = Path(config.output_dir) / "final_grid_flow.data"
        if final_grid_data.exists():
            split_dir = Path(config.output_dir) / "split_states"
            self.split_grid_results(final_grid_data, metadata_original, str(split_dir), spacing=spacing)

    def resume_grid_flow(self, restart_path: str, run_steps: int = 1000000,
                         dt: float = 1e-6,
                         lepton_file: str = "simulation_templates/lepton.inc",
                         dump_file: str = "simulation_templates/default_dump.inc",
                         viscosity: float = 0.000001,
                         num_procs: int = 1, num_threads: int = 1, use_kokkos: bool = True,
                         template: str = "in.grid_hopper_flow_resume",
                         simulation: str = None,
                         seed: int = 42,
                         inplace: bool = False):
        """
        Resumes a grid flow simulation from a restart file.
        """
        import json
        self._ensure_local_directory(restart_path)
        restart_p = Path(restart_path)
        if restart_p.is_dir():
            job_dir = restart_p
            # Find latest restart in directory
            r_file = job_dir / "restart" / "restart.final.bin"
            if not r_file.exists():
                restarts = list((job_dir / "restart").glob("restart.*.bin"))
                if restarts:
                    r_file = sorted(restarts, key=os.path.getmtime)[-1]
                else:
                    raise FileNotFoundError(f"No restart file found in {job_dir}")
            restart_p = r_file
        else:
            job_dir = restart_p.parent
            if job_dir.name == "restart":
                job_dir = job_dir.parent
            
        metadata_path = job_dir / "grid_metadata.json"
        if not metadata_path.exists():
            metadata_path = job_dir.parent / "grid_metadata.json"
        
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found for {restart_path}. Cannot resume flow.")

        with open(metadata_path, 'r') as f:
            meta_raw = json.load(f)
            if simulation is None:
                simulation = meta_raw.get("simulation", "Grid_Hopper_Flow")
            geometry_inc = meta_raw.get("geometry_inc")
            spacing = meta_raw.get("spacing", 2.0)
            
            # Reconstruction parameters for self-contained geometry
            n_hoppers = len(meta_raw["metadata"])
            geo_vars = meta_raw.get("geometry_vars", [])
            hopper_template = meta_raw.get("hopper_template_data", "simulation_geometries/2D_hopper_with_orifice_cover.inc")
            
            freq_list = meta_raw.get("freq", [])
            amp_list = meta_raw.get("amp", [])
            osc_dir = meta_raw.get("osc_dir", "z")

            metadata_original = {int(k): {
                "offset": np.array(v["offset"]),
                "N": v.get("N"),
                "n_fill": v.get("n_fill"),
                "geometry_vars": v.get("geometry_vars")
            } for k, v in meta_raw["metadata"].items()}

        # Naming Strategy
        if inplace:
            new_job_dir = job_dir
            run_name = job_dir.name
        else:
            new_job_dir, run_name = self._generate_resume_path(job_dir, seed, simulation)
            new_job_dir.mkdir(parents=True, exist_ok=True)
        
        # Regeneration of geometry to ensure self-contained jobs
        osc_params = [{"freq": freq_list[i], "amp": amp_list[i], "dir": osc_dir} for i in range(n_hoppers)] if freq_list else None
        
        local_geometry_inc, _ = self._generate_replicated_geometry(
            Path(hopper_template), n_hoppers, spacing, new_job_dir, geo_vars,
            exclude_regions=["orifice_cover"] if "Flow" in simulation else [],
            oscillation_params=osc_params,
            inc_name="replicated_geometry_flow.inc" if "Flow" in simulation else "replicated_geometry.inc"
        )
        
        # Inherit Metadata and Update Lineage
        new_meta = meta_raw.copy()
        if not inplace:
            new_meta["source_dir"] = str(job_dir.absolute()).replace("\\", "/")
            new_meta["run_name"] = run_name
        
        # Track active seeds
        if "active_seeds" not in new_meta: new_meta["active_seeds"] = []
        if str(seed) not in new_meta["active_seeds"]: new_meta["active_seeds"].append(str(seed))

        if inplace:
            self._reset_sync_status(str(new_job_dir.absolute()))

        with open(new_job_dir / "grid_metadata.json", 'w') as f:
            json.dump(new_meta, f, indent=4)

        config = SimulationConfig(
            template=template,
            simulation=simulation,
            run=run_name,
            data_file=str(restart_p).replace("\\", "/"),
            lepton_file=lepton_file,
            dump_file=dump_file,
            outdir_override=str(new_job_dir),
            extra_vars={
                "restart_path": str(restart_p).replace("\\", "/"),
                "geometry_inc": local_geometry_inc.replace("\\", "/"),
                "run_steps": run_steps,
                "dt": dt,
                "viscosity": viscosity
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos
        )

        print(f"--- Resuming Grid Hopper Flow from {restart_p.name} ---")
        self.runner.run(config, clean_dir=False)

        final_grid_data = Path(config.output_dir) / "final_grid_flow.data"
        if final_grid_data.exists():
            split_dir = Path(config.output_dir) / "split_states"
            self.split_grid_results(final_grid_data, metadata_original, str(split_dir), spacing=spacing)

    def _generate_replicated_geometry(self, setup_path: Path, n_hoppers: int, spacing: float, job_dir: Path, 
                                      normalized_geo_vars: List[Dict], exclude_regions: List[str] = None,
                                      oscillation_params: List[Dict] = None, inc_name: str = "replicated_geometry.inc"):
        original_commands = self._resolve_includes(setup_path)
        # Use overrides from the first hopper for the global block if available
        global_overrides = normalized_geo_vars[0] if normalized_geo_vars else {}
        vars_dict = self._evaluate_variables(original_commands, global_overrides)
        
        if normalized_geo_vars:
            all_override_vars = set().union(*[g.keys() for g in normalized_geo_vars])
            missing_vars = all_override_vars - set(vars_dict.keys())
            if missing_vars:
                raise ValueError(f"Error: The following geometry variables were not found in {setup_path}: {missing_vars}")

        # Extract all original region names to ensure we only suffix valid references
        original_region_names = {l.strip().split()[1] for l in original_commands if l.strip().startswith("region")}
        
        cols = math.ceil(math.sqrt(n_hoppers))
        replicated_lines = []
        
        # 2. Write global variables/comments once
        replicated_lines.append("# --- Global Variables (Resolved in Python) ---\n")
        for k, v in vars_dict.items():
            if k != "pi":
                replicated_lines.append(f"variable {k} equal {v}\n")
        
        metadata = {}
        for i in range(n_hoppers):
            ix, iy, iz = i % cols, i // cols, 0
            offset = np.array([ix * spacing, iy * spacing, iz * spacing])
            metadata[i] = {"offset": offset}
            
            replicated_lines.append(f"\n# --- Replicated Hopper {i} at {offset} ---\n")
            
            # Apply per-hopper variable overrides and re-evaluate derived variables
            current_vars = self._evaluate_variables(original_commands, normalized_geo_vars[i])
            
            # 3. Define Oscillation Variables for this Hopper
            if oscillation_params:
                p = oscillation_params[i]
                freq = p.get('freq', 0.0)
                amp = p.get('amp', 0.0)
                if amp > 0:
                    replicated_lines.append(f"variable osc_{i} equal \"{amp}*sin(2*PI*{freq}*step*dt)\"\n")
                else:
                    replicated_lines.append(f"variable osc_{i} equal 0.0\n")

            for line in original_commands:
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or stripped.startswith("variable"):
                    continue
                
                parts = stripped.split()
                name = parts[1]
                
                # Check for region exclusion (e.g., skip orifice_cover for flow)
                if exclude_regions and any(ex in name for ex in exclude_regions):
                    continue

                # Replicate Regions
                if stripped.startswith("region"):
                    style = parts[2]
                    rest = parts[3:]
                    new_name = f"{name}_{i}"
                    
                    def resolve_val(val, off):
                        expr = val
                        for k, v in current_vars.items():
                            expr = expr.replace(f"${{{k}}}", str(v))
                        try:
                            clean_expr = expr.replace("math.", "")
                            return float(eval(clean_expr, {"math": math, "__builtins__": None})) + off
                        except:
                            return f"{val}+{off}"

                    move_str = ""
                    if oscillation_params and style in ["block", "plane", "cylinder", "sphere", "prism"]:
                        p = oscillation_params[i]
                        amp = p.get('amp', 0.0)
                        odir = p.get('dir', 'z').lower()
                        if amp > 0:
                            if odir == 'x': move_str = f" move v_osc_{i} NULL NULL"
                            elif odir == 'y': move_str = f" move NULL v_osc_{i} NULL"
                            else: move_str = f" move NULL NULL v_osc_{i}"

                    if style == "block":
                        coords = [str(resolve_val(c, offset[j//2])) for j, c in enumerate(rest[:6])]
                        final_rest = [f"{p}_{i}" if p in original_region_names else p for p in rest[6:]]
                        new_line = f"region {new_name} block {' '.join(coords)} {' '.join(final_rest)}{move_str}\n"
                    elif style == "plane":
                        coords = [str(resolve_val(c, offset[j])) for j, c in enumerate(rest[:3])]
                        # Resolve normal components without offset
                        normals = [str(resolve_val(c, 0)) for c in rest[3:6]]
                        final_rest = [f"{p}_{i}" if p in original_region_names else p for p in rest[6:]]
                        new_line = f"region {new_name} plane {' '.join(coords)} {' '.join(normals)} {' '.join(final_rest)}{move_str}\n"
                    elif style == "cylinder":
                        # region ID cylinder dim c1 c2 radius lo hi
                        dim = rest[0]
                        c1_idx = 1 if dim == 'x' else 0
                        c2_idx = 2 if dim == 'z' else 1
                        c1 = str(resolve_val(rest[1], offset[c1_idx]))
                        c2 = str(resolve_val(rest[2], offset[c2_idx]))
                        radius = str(resolve_val(rest[3], 0))
                        lo = str(resolve_val(rest[4], offset["xyz".find(dim)]))
                        hi = str(resolve_val(rest[5], offset["xyz".find(dim)]))
                        final_rest = [f"{p}_{i}" if p in original_region_names else p for p in rest[6:]]
                        new_line = f"region {new_name} cylinder {dim} {c1} {c2} {radius} {lo} {hi} {' '.join(final_rest)}{move_str}\n"
                    elif style in ["intersect", "union"]:
                        num = rest[0]
                        regs = [f"{r}_{i}" if r in original_region_names else r for r in rest[1:]]
                        new_line = f"region {new_name} {style} {num} {' '.join(regs)}{move_str}\n"
                    else:
                        # Fallback for other styles: resolve everything as best as possible
                        resolved_rest = []
                        for part in rest:
                            if part in original_region_names:
                                resolved_rest.append(f"{part}_{i}")
                            else:
                                resolved_rest.append(str(resolve_val(part, 0)))
                        new_line = f"region {new_name} {style} {' '.join(resolved_rest)}{move_str}\n"
                    replicated_lines.append(new_line)
                
                # Replicate Fixes
                elif stripped.startswith("fix"):
                    f_id = parts[1]
                    new_id = f"{f_id}_{i}"
                    
                    # Check if this fix references an excluded region
                    skip_fix = False
                    if exclude_regions:
                        for p in parts:
                            if any(ex in p for ex in exclude_regions):
                                skip_fix = True
                                break
                    if skip_fix:
                        continue

                    new_parts = []
                    for p in parts:
                        if p in original_region_names:
                            new_parts.append(f"{p}_{i}")
                        else:
                            new_parts.append(p)
                    new_parts[1] = new_id
                    replicated_lines.append(f"{' '.join(new_parts)}\n")
                    
        rows = math.ceil(n_hoppers / cols)
        total_bounds = {
            'x': [-spacing, cols * spacing],
            'y': [-spacing, rows * spacing],
            'z': [-0.1, 1.0]  # Default Z bounds
        }
        
        envelope = {
            "total_bounds": total_bounds,
            "hoppers": metadata
        }
        
        out_path = job_dir / inc_name
        with open(out_path, 'w') as f:
            f.writelines(replicated_lines)
        return str(out_path).replace("\\", "/"), envelope

    def _evaluate_variables(self, commands: List[str], overrides: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Parses and evaluates LAMMPS-style variables from a list of command lines.
        Supports basic math and variable substitution.
        """
        import math
        vars_dict = {"pi": math.pi}
        if overrides:
            vars_dict.update(overrides)
            
        for line in commands:
            stripped = line.strip()
            if stripped.startswith("variable"):
                parts = stripped.split()
                if len(parts) < 4: continue
                v_name = parts[1]
                
                # If this variable was explicitly overridden by the user, we skip the template definition
                if overrides and v_name in overrides:
                    continue
                    
                v_expr = " ".join(parts[3:])
                v_expr = v_expr.split("#")[0].strip()
                
                # Substitute known variables
                # Sort keys by length descending to prevent partial replacements (e.g., 'y_half' matching inside 'h_y_half')
                for k in sorted(vars_dict.keys(), key=len, reverse=True):
                    v = vars_dict[k]
                    v_expr = v_expr.replace(f"${{{k}}}", str(v))
                
                try:
                    # Map common LAMMPS math functions to Python's math module
                    safe_expr = v_expr.replace("sin(", "math.sin(").replace("cos(", "math.cos(").replace("tan(", "math.tan(")
                    vars_dict[v_name] = eval(safe_expr, {"math": math, "__builtins__": None}, vars_dict)
                except:
                    # Fallback for expressions we can't parse in Python
                    pass
        return vars_dict

    def _generate_grid_insertion_file(self, outdir: Path, n_hoppers: int, n_fill_list: List[int], 
                                      mol_ranges: Dict[int, Any], seed: int, N_list: List[int], 
                                      spacing: float, metadata: Dict[int, Any],
                                      mol_bboxes: Dict[int, Any] = None,
                                      mode: str = "2D_stacked",
                                      z_start: float = 0.51):
        rng = random.Random(seed)
        insertion_lines = []
        
        # Internal pouring grid relative to hopper center
        # y_half from 2D_hopper.inc is 0.155
        hopper_y_half = 0.155
        z_max_global = z_start
        total_inserted = 0
        
        for h_idx in range(n_hoppers):
            offset = metadata[h_idx]['offset']
            n_fill = n_fill_list[h_idx]
            N = N_list[h_idx]
            
            # Get molecule range for this N
            m_range = mol_ranges.get(N)
            if not m_range:
                print(f"Error: No molecules prepared for N={N} (Hopper {h_idx})")
                continue
            
            m_start = m_range['start']
            m_count = m_range['count']

            # Hopper-specific geometry based on N
            bead_diam = 0.003            
            count = 0
            
            if mode == "2D_worst_case":
                safe_buffer = max(0.005, (N * bead_diam / 2.0) + 0.002)
            
                y_max = hopper_y_half - safe_buffer
                y_min = -hopper_y_half + safe_buffer
                if y_max <= y_min:
                    y_max, y_min = 0.001, -0.001
                y_width = y_max - y_min

                dy_gap = N * 0.003 
                dz_gap = N * 0.003
                ny = max(1, int(y_width / dy_gap))
                nz = math.ceil(n_fill / ny)
                
                for iz in range(nz):
                    if count >= n_fill: break
                    pz = z_start + iz * dz_gap
                    z_max_global = max(z_max_global, pz + offset[2])
                    for iy in range(ny):
                        if count >= n_fill: break
                        px = offset[0]
                        py = y_min + (iy + 0.5) * (y_width/ny) + offset[1]
                        mol_id = m_start + rng.randint(0, m_count - 1)
                        insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz+offset[2]:.6f} mol m{mol_id} 12345 rotate 0.0 0.0 0.0 1.0")
                        count += 1
                        total_inserted += 1
            elif mode == "2D_stacked":
                safe_buffer = bead_diam 
                # safe width to pour in
                y_max = hopper_y_half - safe_buffer
                y_min = -hopper_y_half + safe_buffer
                if y_max <= y_min:
                    y_max, y_min = 0.001, -0.001
                y_width = y_max - y_min

                # Smart packing using actual molecule bounding boxes
                current_z = z_start
                buffer = 0.002 # 2mm safety gap between bounding boxes
                
                while count < n_fill:
                    current_y = y_min
                    max_h_in_row = 0
                    row_empty = True
                    
                    while count < n_fill:
                        # We need to peek at the next molecule's bbox
                        # To keep it simple, we'll pick the molecule first
                        mol_id = m_start + rng.randint(0, m_count - 1)
                        bbox = mol_bboxes[mol_id]
                        m_w = bbox['width'] + buffer
                        m_h = bbox['height'] + buffer
                        
                        # Check if it fits in current row
                        if current_y + m_w > y_max:
                            if row_empty:
                                # Even a single molecule doesn't fit? 
                                # This happens if safe_buffer is too large or hopper too narrow.
                                # Force it into center and move to next row.
                                px = offset[0] - (bbox['x'][0] + bbox['x'][1])/2
                                py = offset[1] - (bbox['y'][0] + bbox['y'][1])/2
                                pz = current_z + buffer/2 - bbox['z'][0] + offset[2]
                                insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz:.6f} mol m{mol_id} 12345 rotate 0.0 0.0 0.0 1.0")
                                z_max_global = max(z_max_global, pz + bbox['z'][1])
                                current_z += m_h
                                count += 1
                                total_inserted += 1
                                break
                            else:
                                # Row is full, move to next Z level
                                current_z += max_h_in_row
                                break # Exit inner loop to start new row
                        
                        # Place at current_y (centered in its allocated slot)
                        px = offset[0] - (bbox['x'][0] + bbox['x'][1])/2
                        py = current_y + buffer/2 - bbox['y'][0] + offset[1]
                        pz = current_z + buffer/2 - bbox['z'][0] + offset[2]
                        insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz:.6f} mol m{mol_id} 12345 rotate 0.0 0.0 0.0 1.0")
                        z_max_global = max(z_max_global, pz + bbox['z'][1])
                        
                        current_y += m_w
                        max_h_in_row = max(max_h_in_row, m_h)
                        row_empty = False
                        count += 1
                        total_inserted += 1
            else:
                # Default 3D grid-like pouring
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
                            mol_id = m_start + rng.randint(0, m_count - 1)
                            mol_seed = seed + total_inserted
                            insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz+offset[2]:.6f} mol m{mol_id} {mol_seed}")
                            count += 1
                            total_inserted += 1
                        
        insertion_path = outdir / "insertions.inc"
        with open(insertion_path, 'w') as f:
            f.write("\n".join(insertion_lines))
        return str(insertion_path).replace("\\", "/"), z_max_global

    def _reset_sync_status(self, run_id: str):
        """Resets the sync status in lineage.json to Local."""
        try:
            from Pulse.pulse_core import PBSManager
            lineage = PBSManager.load_lineage()
            if run_id in lineage:
                lineage[run_id]["sync_status"] = "Local"
                PBSManager.save_lineage(lineage)
                print(f"[INFO] Reset sync status for {run_id} to Local due to in-place resumption.")
        except Exception as e:
            print(f"[WARNING] Could not reset sync status: {e}")

    def split_grid_results(self, big_data_path: Path, metadata: Dict[int, Any], output_base_dir: str, spacing: float = 2.0):
        # [Splitting logic similar to GridRelaxManager but aware of hopper boundaries]
        print(f"Splitting grid filled results into {len(metadata)} hoppers...")
        all_atoms, all_bonds, all_angles = self._parse_data_file(big_data_path)
        
        for h_idx, meta in metadata.items():
            offset = meta['offset']
            n_val = meta.get('N', 'unknown')
            fill_val = meta.get('n_fill', 'unknown')
            
            # Find atoms belonging to this hopper (within a bounding box around the offset)
            # We use a filter slightly smaller than half-spacing to avoid boundary atoms
            spatial_filter = spacing * 0.48 
            h_atoms = [a for a in all_atoms if abs(a['x'] - offset[0]) < spatial_filter and abs(a['y'] - offset[1]) < spatial_filter]
            
            target_dir = Path(output_base_dir) / f"Hopper_{h_idx}"
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
                
            header = f"# LAMMPS Data - Split Hopper {h_idx} (N={n_val}, n_fill={fill_val})\n"
            self._write_data(out_path, h_atoms, bonds, angles, header=header)
            
            # Save a small properties file for easy analysis access
            import json
            with open(target_dir / "properties.json", 'w') as f:
                json.dump({
                    "N": n_val, 
                    "n_fill": fill_val, 
                    "offset": offset.tolist() if hasattr(offset, "tolist") else offset,
                    "geometry_overrides": meta.get('geometry_vars', {}),
                    "source_dir": meta.get('source_dir')
                }, f)
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

    def _write_data(self, path: Path, atoms: List[Dict], bonds: List[Dict], angles: List[Dict], header: str = "# LAMMPS Data\n"):
        with open(path, 'w') as f:
            f.write(header + "\n")
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
