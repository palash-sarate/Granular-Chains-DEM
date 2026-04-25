import os
import shutil
import random
import concurrent.futures
from pathlib import Path
from typing import Optional, List
from .runner import SimulationRunner
from .config import SimulationConfig
from .chain_generator import ChainConfig, write_chain_data

class LibraryGenerator:
    def __init__(self, runner: SimulationRunner, forced: bool = True):
        self.runner = runner
        # If False, skip generation if final relaxed state already exists in target directory
        self.forced = forced

    def run_single_state_task(self, n_beads: int, state_index: int, 
                             target_dir: Path, rel_data_path: str,
                             template: str, simulation_name: str, 
                             run_name_prefix: Optional[str],
                             dump_inc: str, num_procs: int, num_threads: int,
                             use_kokkos: bool, use_intel: bool):
        """Task for a single simulation run, suitable for parallel execution."""
        seed = random.randint(1, 999999)
        prefix = run_name_prefix if run_name_prefix else f"relax_N{n_beads}"
        run_name = f"{prefix}_state_{state_index}"
        
        final_dest = target_dir / f"state_{state_index}.data"
        if not self.forced and final_dest.exists():
            return f"[N={n_beads}, State {state_index}] Skipping (already exists)"
        
        sim_config = SimulationConfig(
            template=template,
            data_file=str(rel_data_path).replace("\\", "/"),
            dump_file=dump_inc, 
            simulation=simulation_name,
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
        
        try:
            from .runner import SimulationRunner
            runner = SimulationRunner(lammps_executable=self.runner.lammps_exe)
            runner.run(sim_config, verbose=False)
            
            generated_file = Path(f"dumping_yard/{simulation_name}/{run_name}/relaxed_chain.data")
            if generated_file.exists():
                shutil.copy(generated_file, final_dest)
                return f"[N={n_beads}, State {state_index}] Saved: {final_dest}"
            else:
                return f"[N={n_beads}, State {state_index}] Error: Output not found"
        except Exception as e:
            return f"[N={n_beads}, State {state_index}] Failed: {e}"

    def generate_library(self, n_beads: int, n_states: int,
                         output_dir: str = "chain_data/relaxed",
                         template: str = "in.relax_3d_gen",
                         run_name_prefix: Optional[str] = None,
                         simulation_name: str = "Relax_Library_Gen",
                         dump_inc: str = "simulation_templates/default_dump.inc",
                         n_parallel: int = 1,
                         num_procs: int = None,
                         num_threads: int = 1,
                         use_kokkos: bool = True,
                         use_intel: bool = True):
        """Generates a library of relaxed chain states for a single N."""
        return self.generate_batch_library(
            n_beads_list=[n_beads],
            n_states_per_n=n_states,
            output_dir=output_dir,
            template=template,
            run_name_prefix=run_name_prefix,
            simulation_name=simulation_name,
            dump_inc=dump_inc,
            n_parallel=n_parallel,
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )

    def generate_batch_library(self, n_beads_list: List[int], n_states_per_n: int,
                               output_dir: str = "chain_data/relaxed",
                               template: str = "in.relax_3d_gen",
                               run_name_prefix: Optional[str] = None,
                               simulation_name: str = "Relax_Library_Gen",
                               dump_inc: str = "simulation_templates/default_dump.inc",
                               n_parallel: int = 1,
                               num_procs: int = None,
                               num_threads: int = 1,
                               use_kokkos: bool = True,
                               use_intel: bool = True):
        """Generates a library for multiple N values in a single parallel pool."""
        
        # 1. Prepare all base data files first
        temp_chain_dir = Path("chain_data/lib_gen_temp")
        temp_chain_dir.mkdir(parents=True, exist_ok=True)
        
        tasks = []
        for n_beads in n_beads_list:
            target_dir = Path(output_dir) / f"N{n_beads}"
            target_dir.mkdir(parents=True, exist_ok=True)
            
            chain_cfg = ChainConfig(beads=n_beads, spacing=0.0025, mode="linear", 
                                   orientation="vert", output_dir=temp_chain_dir)
            base_chain_path = write_chain_data(chain_cfg)
            rel_data_path = base_chain_path.relative_to("chain_data")
            
            prefix = run_name_prefix
            if prefix and "{N}" in prefix:
                prefix = prefix.replace("{N}", str(n_beads))
            elif not prefix:
                prefix = f"relax_N{n_beads}"

            for i in range(n_states_per_n):
                tasks.append({
                    "n_beads": n_beads,
                    "state_index": i,
                    "target_dir": target_dir,
                    "rel_data_path": rel_data_path,
                    "template": template,
                    "simulation_name": simulation_name,
                    "run_name_prefix": prefix,
                    "dump_inc": dump_inc,
                    "num_procs": num_procs,
                    "num_threads": num_threads,
                    "use_kokkos": use_kokkos,
                    "use_intel": use_intel
                })

        print(f"Batch generating {len(tasks)} states across N={n_beads_list}...")
        print(f"Using {n_parallel} workers.")

        if n_parallel > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as executor:
                future_to_task = {executor.submit(self.run_single_state_task, **task): task for task in tasks}
                for future in concurrent.futures.as_completed(future_to_task):
                    result = future.result()
                    print(f"  {result}")
        else:
            for task in tasks:
                result = self.run_single_state_task(**task)
                print(f"  {result}")

        print("Batch library generation complete.")
