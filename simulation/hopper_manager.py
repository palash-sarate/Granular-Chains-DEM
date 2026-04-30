import os
import glob
from pathlib import Path
from typing import Optional, Any
from .molecule_converter import convert_data_to_molecule
from .runner import SimulationRunner
from .config import SimulationConfig

class HopperManager:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def prepare_molecules(self, source_dir: str, mol_dir: str) -> str:
        """
        Converts all .data files in source_dir to .mol files in mol_dir.
        Generates a molecules.inc file.
        Returns the path to molecules.inc.
        """
        source_path = Path(source_dir)
        mol_path = Path(mol_dir)
        mol_path.mkdir(parents=True, exist_ok=True)
        
        data_files = list(source_path.glob("*.data"))
        if not data_files:
            raise ValueError(f"No .data files found in {source_dir}")
            
        inc_lines = []
        
        print(f"Converting {len(data_files)} chains from {source_dir} to molecules in {mol_dir}...")
        
        for i, data_file in enumerate(data_files):
            mol_id = i + 1
            mol_filename = f"mol_{mol_id}.mol"
            output_mol = mol_path / mol_filename
            
            # Convert
            convert_data_to_molecule(str(data_file), str(output_mol))
            
            # Add to include file
            # Use forward slashes for LAMMPS and relative path from project root
            mol_rel_path = str(output_mol).replace("\\", "/")
            inc_lines.append(f"molecule m{mol_id} {mol_rel_path}")
            
        inc_file = mol_path / "molecules.inc"
        with open(inc_file, 'w') as f:
            f.write("\n".join(inc_lines))
            
        print(f"Generated molecule index: {inc_file}")
        return str(inc_file).replace("\\", "/")

    def run_hopper_flow(self, 
                        chain_source_dir: str, 
                        n_fill: int, 
                        drop_steps: int = None,
                        freq: float = 5.0, 
                        amp: float = 0.005,
                        run_name: str = "hopper_test",
                        dt: float = 1e-6,
                        run_steps: int = 10000,
                        setup_inc: str = "simulation_geometries/2D_hopper_flow_setup.inc",
                        num_procs: int = None, 
                        num_threads: int = 1
                        ):
        
        if drop_steps is None:
            drop_steps = 5 * (1 / dt) * 1e-2
            # print(f"Setting drop_steps to {drop_steps} based on dt={dt}")
            
        # 1. Prepare Molecules
        # We'll store molecules in a subdir of the source or a temp dir
        mol_dir = "chain_data/molecules_temp"
        inc_file = self.prepare_molecules(chain_source_dir, mol_dir)
        
        # Count templates
        n_templates = len(list(Path(chain_source_dir).glob("*.data")))
        
        # 2. Configure Simulation
        sim_config = SimulationConfig(
            template="in.hopper_flow",
            simulation="Hopper_Flow",
            run=run_name,
            extra_vars={
                "setup_inc": setup_inc,
                "mol_include_file": inc_file,
                "n_templates": n_templates,
                "n_fill": n_fill,
                "freq": freq,
                "amp": amp,
                "seed": 12345,
                "run_steps": run_steps, # Post-fill run
                "dt": dt,
                "drop_steps": int(drop_steps)
            },
            num_procs=num_procs,
            num_threads=num_threads
        )
        
        # 3. Run
        print(f"Starting Hopper Flow simulation: {run_name}")
        # Enable directory cleaning to prevent mixing old and new data
        self.runner.run(sim_config, verbose=True, clean_dir=True)

    def generate_filled_state(self, source_dir: str, n_fill: int, relax_steps: int,
                              dt: float = 1e-6, run_name: str = None, seed: int = 12345,
                              mol_dir: Optional[str] = None, setup_inc: str = "",
                              dump_inc: str = "simulation_templates/default_dump.inc", 
                              viscosity: float = 0.001, N: int = 4, fill_template: str = "in.hopper_fill",
                              outdir: str = None, num_procs: int = None, num_threads: int = 1,
                              use_kokkos: bool = True, use_intel: bool = True,
                              clean_dir: bool = True) -> str:
        
        """Create a filled hopper state from relaxed chain files and save data+restart.
        Returns the path to the saved data file (forward-slashes).
        """
        if run_name is None:
            run_name = f"filled_N{n_fill}_s{seed}"

        # 1. Resolve output directory early to place molecules inside
        from .config import SimulationConfig
        actual_outdir = outdir if outdir else SimulationConfig.compute_output_dir("Hopper_Fill", run_name)
        
        if mol_dir is None or mol_dir == "chain_data/molecules_temp":
            mol_dir = os.path.join(actual_outdir, "molecules")

        # Prepare molecules and include file
        inc_file = self.prepare_molecules(source_dir, mol_dir)
        n_templates = len(list(Path(source_dir).glob("*.data")))

        # Determine default drop_steps similar to run_hopper_flow
        # drop_steps = int(5 * (1 / dt) * 1e-2)

        cfg = SimulationConfig(
            template=fill_template,
            simulation="Hopper_Flow",
            run=run_name,
            data_file=None,
            dump_file=dump_inc,
            outdir_override=outdir,
            extra_vars={
                "viscosity": viscosity,
                "setup_inc": setup_inc,
                "mol_include_file": inc_file,
                "n_templates": n_templates,
                "n_fill": n_fill,
                # "drop_steps": int(drop_steps),
                "relax_steps": relax_steps,
                "dt": dt,
                "seed": seed,
                "N": N,
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )

        # Prepare directories manually first so we can place the insertion file inside
        if clean_dir:
            self.runner._prepare_directories(cfg, clean=True)

        # If using the tall template, pre-generate the insertion file via Python
        if "tall" in fill_template:
            output_dir = cfg.output_dir
            insertion_file, z_max = self._generate_tall_insertion_file(output_dir, n_fill, n_templates, seed, N)
            cfg.extra_vars["insertion_file"] = insertion_file
            cfg.extra_vars["z_max"] = z_max

        # Tell runner not to clean/prep again (prep_dirs=False)
        self.runner.run(cfg, verbose=True, clean_dir=False, prep_dirs=False)

        saved_data = f"{cfg.output_dir}/final_hopper.data".replace("\\", "/")
        return saved_data

    def resume_filled_state(self, run_name: str, source_dir: str, n_fill: int, relax_steps: int,
                            dt: float = 1e-6, restart_path: str = None, seed: int = 12345,
                            mol_dir: Optional[str] = None, setup_inc: str = "",
                            dump_inc: str = "simulation_templates/default_dump.inc", 
                            viscosity: float = 0.001, fill_template: str = "in.hopper_fill_resume",
                            N: int = 4, outdir: str = None, num_procs: int = None, num_threads: int = 1,
                            use_kokkos: bool = True, use_intel: bool = True) -> str:
        """Resume a hopper fill simulation from a restart file and add more chains.
        Returns the path to the saved data file.
        """
        
        if restart_path is None:
            # Look for restart files in potential locations:
            # 1. The custom outdir (if provided)
            # 2. The default directory computed from run_name
            search_dirs = []
            if outdir:
                search_dirs.append(outdir)
            
            default_dir = SimulationConfig.compute_output_dir("Hopper_Flow", run_name)
            if default_dir not in search_dirs:
                search_dirs.append(default_dir)

            found_restart = False
            for d in search_dirs:
                restart_search_dir = os.path.join(d, "restart").replace("\\", "/")
                
                if not os.path.exists(restart_search_dir):
                    continue
                
                restart_files = [f for f in os.listdir(restart_search_dir) if f.endswith(".bin")]
                if not restart_files:
                    continue
                    
                restart_files.sort(key=lambda x: os.path.getmtime(os.path.join(restart_search_dir, x)))
                restart_path = os.path.join(restart_search_dir, restart_files[-1]).replace("\\", "/")
                print(f"Using latest restart file found in {d}: {restart_path}")
                found_restart = True
                break
                
            if not found_restart:
                 raise FileNotFoundError(f"Could not find any restart files (.bin) in searched locations: {search_dirs}")

        # Ensure unique mol_dir for resumption too
        if mol_dir is None or mol_dir == "chain_data/molecules_temp":
            # outdir is handled in SimulationConfig compute_output_dir if not provided
            from .config import SimulationConfig
            actual_outdir = outdir if outdir else SimulationConfig.compute_output_dir("Hopper_Fill", run_name)
            mol_dir = os.path.join(actual_outdir, "molecules")

        # Prepare molecules (ensure we have the .mol files for create_atoms)
        inc_file = self.prepare_molecules(source_dir, mol_dir)
        n_templates = len(list(Path(source_dir).glob("*.data")))

        # Determine drop_steps
        # drop_steps = int(5 * (1 / dt) * 1e-2)

        cfg = SimulationConfig(
            template=fill_template,
            simulation="Hopper_Flow",
            run=run_name,
            resume_file=restart_path,
            dump_file=dump_inc,
            outdir_override=outdir,
            extra_vars={
                "viscosity": viscosity,
                "setup_inc": setup_inc,
                "mol_include_file": inc_file,
                "n_templates": n_templates,
                "n_fill": n_fill,
                # "drop_steps": int(drop_steps),
                "relax_steps": relax_steps,
                "dt": dt,
                "seed": seed,
                "N": N,
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )

        print(f"Resuming filled hopper state: {run_name}")
        # Use runner.resume instead of runner.run
        self.runner.resume(cfg, verbose=True)

        saved_data = f"{cfg.output_dir}/final_hopper_resume.data".replace("\\", "/")
        return saved_data

    def run_flow_from_saved(self, saved_data_path: str = None, restart_path: str = None,
                            run_name: str = None, freq: float = 5.0, amp: float = 0.005,
                            dt: float = 1e-6, run_steps: int = 10000, num_procs: int = None, num_threads: int = 1) -> str:
        """Run an oscillating hopper flow starting from a saved data or restart file.
        Provide either `saved_data_path` or `restart_path` (restart preferred).
        Returns the output directory path.
        """
        if run_name is None:
            run_name = "flow_from_saved"

        extra = {
            "dt": dt,
            "v_freq": freq,
            "v_amp": amp,
            "run_steps": run_steps,
        }

        cfg = SimulationConfig(
            template="in.hopper_run_from_saved",
            simulation="Hopper_Flow",
            run=run_name,
            data_file=saved_data_path or "",
            resume_file=restart_path or None,
            extra_vars=extra,
            num_procs=num_procs,
            num_threads=num_threads
        )

        print(f"Running hopper flow from saved state: {run_name}")
        self.runner.run(cfg, verbose=True, clean_dir=True)
        return cfg.output_dir

    def _generate_tall_insertion_file(self, outdir: str, n_fill: int, n_templates: int, seed: int, N: int):
        """
        Generates a 3D grid of insertion coordinates in a tall column via Python.
        Writes 'create_atoms' commands to a file to be included by LAMMPS.
        Returns (insertion_file_path, z_max).
        """
        import math
        import random
        
        # Use a local random instance to avoid interfering with other parts of the program
        rng = random.Random(seed)
        
        # Grid bounds (matched to in.hopper_fill_uniform)
        x_min, x_max = -0.08, 0.08
        y_min, y_max = -0.15, 0.15
        z_start = 0.51
        
        x_width = x_max - x_min
        y_width = y_max - y_min
        
        # Exclusion diameter (3 * N * 1mm)
        ex_diam = 3 * N * 0.001
        
        nx = max(1, int(x_width / ex_diam))
        ny = max(1, int(y_width / ex_diam))
        grid_2d = nx * ny
        
        nz = math.ceil(n_fill / grid_2d)
        dz = ex_diam # Vertical spacing
        
        insertion_lines = []
        count = 0
        z_max = z_start
        
        for iz in range(nz):
            if count >= n_fill: break
            
            # Random offset for this layer to avoid perfect vertical alignment
            ox = rng.uniform(0, x_width/nx)
            oy = rng.uniform(0, y_width/ny)
            
            pz = z_start + iz * dz
            z_max = max(z_max, pz)
            
            for ix in range(nx):
                if count >= n_fill: break
                for iy in range(ny):
                    if count >= n_fill: break
                    
                    px = x_min + (ix + 0.5) * (x_width/nx) + ox
                    py = y_min + (iy + 0.5) * (y_width/ny) + oy
                    
                    # Wrap x/y if offset pushed them out
                    if px > x_max: px -= x_width
                    if py > y_max: py -= y_width
                    
                    mol_id = rng.randint(1, n_templates)
                    insertion_lines.append(f"create_atoms 0 single {px:.6f} {py:.6f} {pz:.6f} mol m{mol_id} 12345 rotate 0.0 0.0 0.0 1.0")
                    count += 1
        
        insertion_path = os.path.join(outdir, "insertions.inc").replace("\\", "/")
        os.makedirs(outdir, exist_ok=True)
        with open(insertion_path, 'w') as f:
            f.write("\n".join(insertion_lines))
            
        print(f"Generated 3D insertion file with {count} molecules: {insertion_path} (z_max={z_max:.3f})")
        return insertion_path, z_max
