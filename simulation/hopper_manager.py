import os
import glob
from pathlib import Path
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
        
        print(f"Converting {len(data_files)} chains from {source_dir} to molecules...")
        
        for i, data_file in enumerate(data_files):
            mol_id = i + 1
            mol_filename = f"mol_{mol_id}.mol"
            output_mol = mol_path / mol_filename
            
            # Convert
            convert_data_to_molecule(str(data_file), str(output_mol))
            
            # Add to include file
            # Use forward slashes for LAMMPS
            mol_rel_path = f"{mol_dir}/{mol_filename}".replace("\\", "/")
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
                        setup_inc: str = "simulation_geometries/2D_hopper_flow_setup.inc"
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
            }
        )
        
        # 3. Run
        print(f"Starting Hopper Flow simulation: {run_name}")
        # Enable directory cleaning to prevent mixing old and new data
        self.runner.run(sim_config, verbose=True, clean_dir=True)

    def generate_filled_state(self, source_dir: str, n_fill: int, relax_steps: int,
                              dt: float = 1e-6, run_name: str = None, seed: int = 12345,
                              mol_dir: str = "chain_data/molecules_temp", setup_inc: str = "",
                              dump_inc: str = "simulation_templates/default_dump.inc") -> str:
        
        """Create a filled hopper state from relaxed chain files and save data+restart.
        Returns the path to the saved data file (forward-slashes).
        """
        if run_name is None:
            run_name = f"filled_N{n_fill}_s{seed}"

        # Prepare molecules and include file
        inc_file = self.prepare_molecules(source_dir, mol_dir)
        n_templates = len(list(Path(source_dir).glob("*.data")))

        # Determine default drop_steps similar to run_hopper_flow
        drop_steps = int(5 * (1 / dt) * 1e-2)

        cfg = SimulationConfig(
            template="in.hopper_fill",
            simulation="Hopper_Flow",
            run=run_name,
            data_file=None,
            extra_vars={
                "viscosity": 0.001,
                "setup_inc": setup_inc,
                "dump_inc": dump_inc,
                "mol_include_file": inc_file,
                "n_templates": n_templates,
                "n_fill": n_fill,
                "drop_steps": int(drop_steps),
                "relax_steps": relax_steps,
                "dt": dt,
                "seed": seed,
            }
        )

        print(f"Generating filled hopper state: {run_name}")
        self.runner.run(cfg, verbose=True, clean_dir=True)

        saved_data = f"{cfg.output_dir}/final_hopper.data".replace("\\", "/")
        return saved_data

    def run_flow_from_saved(self, saved_data_path: str = None, restart_path: str = None,
                            run_name: str = None, freq: float = 5.0, amp: float = 0.005,
                            dt: float = 1e-6, run_steps: int = 10000) -> str:
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
            extra_vars=extra
        )

        print(f"Running hopper flow from saved state: {run_name}")
        self.runner.run(cfg, verbose=True, clean_dir=True)
        return cfg.output_dir

