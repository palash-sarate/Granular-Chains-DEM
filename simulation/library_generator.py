import os
import shutil
import random
from pathlib import Path
from .runner import SimulationRunner
from .config import SimulationConfig
from .chain_generator import ChainConfig, write_chain_data

class LibraryGenerator:
    def __init__(self, runner: SimulationRunner):
        self.runner = runner

    def generate_library(self, n_beads: int, n_states: int, output_dir: str = "chain_data/relaxed"):
        """
        Generates a library of relaxed chain states.
        
        Args:
            n_beads: Number of beads in the chain.
            n_states: Number of independent states to generate.
            output_dir: Base directory to store the library.
        """
        # 1. Setup directories
        target_dir = Path(output_dir) / f"N{n_beads}"
        target_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Generating {n_states} relaxed states for N={n_beads} in {target_dir}...")

        # 2. Generate Base Linear Chain (Temporary)
        # We place it in chain_data/temp_lib_gen so runner.py can find it easily
        temp_chain_dir = Path("chain_data/temp_lib_gen")
        temp_chain_dir.mkdir(parents=True, exist_ok=True)
        
        chain_cfg = ChainConfig(
            beads=n_beads,
            spacing=0.0025,
            mode="linear",
            orientation="vert",
            output_dir=temp_chain_dir
        )
        base_chain_path = write_chain_data(chain_cfg)
        
        # Calculate relative path for SimulationConfig (relative to chain_data/)
        # base_chain_path is like "chain_data/temp_lib_gen/N...data"
        # we need "temp_lib_gen/N...data"
        rel_data_path = base_chain_path.relative_to("chain_data")
        
        # 3. Loop to generate states
        for i in range(n_states):
            seed = random.randint(1, 999999)
            run_name = f"relax_N{n_beads}_state_{i}"
            
            # Config for relaxation
            sim_config = SimulationConfig(
                template="in.relax_3d_gen_template",
                data_file=str(rel_data_path).replace("\\", "/"), # Ensure forward slashes for LAMMPS
                simulation="Relax_3d_Library_Gen",
                run=run_name,
                extra_vars={
                    "seed": seed,
                    "run_steps": 50000, # 0.05s of relaxation
                    "dt": 1e-6,
                    "temperature": 1e12
                }
            )
            
            print(f"  [State {i+1}/{n_states}] Running relaxation (Seed: {seed})...")
            try:
                self.runner.run(sim_config, verbose=False)
                
                # 4. Move and Rename Output
                # The output will be in dumping_yard/Library_Gen/relax_N.../relaxed_chain.data
                generated_file = Path(f"dumping_yard/Relax_3d_Library_Gen/{run_name}/relaxed_chain.data")
                
                if generated_file.exists():
                    final_dest = target_dir / f"state_{i}.data"
                    shutil.copy(generated_file, final_dest)
                    print(f"    -> Saved: {final_dest}")
                else:
                    print(f"    -> Error: Output file not found for state {i}")
            except Exception as e:
                print(f"    -> Simulation failed: {e}")

        # 5. Cleanup
        # shutil.rmtree(temp_chain_dir) 
        print(f"Library generation for N={n_beads} complete.")
