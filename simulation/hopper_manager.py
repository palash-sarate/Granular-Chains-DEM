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
                        freq: float = 5.0, 
                        amp: float = 0.005,
                        run_name: str = "hopper_test",
                        dt: float = 1e-6,
                        run_steps: int = 10000,
                        temperature: float = 1e12):
        
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
                "mol_include_file": inc_file,
                "n_templates": n_templates,
                "n_fill": n_fill,
                "freq": freq,
                "amp": amp,
                "seed": 12345,
                "run_steps": run_steps, # Post-fill run
                "dt": dt,
                "temperature": temperature # Not used for creation but for thermostat if needed
            }
        )
        
        # 3. Run
        print(f"Starting Hopper Flow simulation: {run_name}")
        self.runner.run(sim_config, verbose=True)

