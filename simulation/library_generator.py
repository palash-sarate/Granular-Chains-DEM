import os
import shutil
import random
import concurrent.futures
from pathlib import Path
from .runner import SimulationRunner
from .config import SimulationConfig
from .chain_generator import ChainConfig, write_chain_data

class LibraryGenerator:
    def __init__(self, runner: SimulationRunner, forced: bool = True):
        self.runner = runner
        # If False, skip generation if final relaxed state already exists in target directory
        self.forced = forced

    def generate_library(self, n_beads: int, n_states: int,
                         output_dir: str = "chain_data/relaxed",
                         dump_inc: str = "simulation_templates/default_dump.inc",
                         n_parallel: int = 1,
                         num_procs: int = None,
                         num_threads: int = 1,
                         use_kokkos: bool = True,
                         use_intel: bool = True):
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
        temp_chain_dir = Path("chain_data/lib_gen_temp")
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
        
        # 3. Define single state generation task
        def run_single_state(i):
            seed = random.randint(1, 999999)
            run_name = f"relax_N{n_beads}_state_{i}"
            
            final_dest = target_dir / f"state_{i}.data"
            if not self.forced and final_dest.exists():
                print(f"  [State {i+1}/{n_states}] Skipping (already exists): {final_dest}")
                return
            
            # Config for relaxation
            sim_config = SimulationConfig(
                template="in.relax_3d_gen",
                data_file=str(rel_data_path).replace("\\", "/"),
                dump_file=dump_inc, 
                simulation="Relax_3d_Library_Gen",
                run=run_name,
                extra_vars={
                    "seed": seed,
                    "motion_steps": 500000, 
                    "explore_steps": 200000, 
                    "viscous_relax_steps": 300000, 
                    "dt": 1e-6,
                    "temperature": 1e15
                },
                num_procs=num_procs,
                num_threads=num_threads,
                use_kokkos=use_kokkos,
                use_intel=use_intel
            )
            
            print(f"  [State {i+1}/{n_states}] Starting relaxation (Seed: {seed})...")
            try:
                # We need a fresh runner or a thread-safe runner
                # Since SimulationRunner is lightweight and uses subprocess, 
                # we'll use the existing one but be mindful of its state if any.
                # Actually, creating a fresh runner for each thread is safer.
                from .runner import SimulationRunner
                runner = SimulationRunner(lammps_executable=self.runner.lammps_exe)
                runner.run(sim_config, verbose=False)
                
                # 4. Move and Rename Output
                generated_file = Path(f"dumping_yard/Relax_3d_Library_Gen/{run_name}/relaxed_chain.data")
                if generated_file.exists():
                    shutil.copy(generated_file, final_dest)
                    print(f"    -> Saved: {final_dest}")
                else:
                    print(f"    -> Error: Output file not found for state {i}")
            except Exception as e:
                print(f"    -> Simulation failed for state {i}: {e}")

        # 4. Execute in Parallel
        if n_parallel > 1:
            print(f"  Running in parallel with {n_parallel} workers...")
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
                executor.map(run_single_state, range(n_states))
        else:
            for i in range(n_states):
                run_single_state(i)

        # 5. Cleanup
        # shutil.rmtree(temp_chain_dir) 
        print(f"Library generation for N={n_beads} complete.")
