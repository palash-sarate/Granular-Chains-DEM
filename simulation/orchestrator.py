import os
import sys
import time
import ast
from pathlib import Path
from typing import List, Optional, Dict, Any

from simulation import SimulationConfig, SimulationRunner
from simulation.chain_generator import ChainConfig, write_chain_data
from simulation.library_generator import LibraryGenerator
from simulation.hopper_manager import HopperManager
from analysis.utilities import get_dt_token, get_viscosity_token, ETAEstimator

class SimulationOrchestrator:
    def __init__(self, lammps_executable: str = "lmp"):
        self.lammps_executable = lammps_executable

    def run_hopper_fill(self, 
                        source_dir: Optional[str] = None, 
                        fill_template: str = "in.hopper_fill",
                        n_fill: int = 10, 
                        relax_steps: int = 100000, 
                        run_name: Optional[str] = None, 
                        seed: Optional[int] = None, 
                        dt: float = 1e-6,
                        mol_dir: str = "chain_data/molecules_temp",
                        setup_inc: str = "simulation_geometries/2D_hopper.inc",
                        dump_inc: str = "simulation_templates/default_dump.inc",
                        viscosity: float = 0.001, 
                        N: int = 4,
                        outdir: Optional[str] = None,
                        num_procs: int = 1,
                        num_threads: int = 1,
                        use_kokkos: bool = True,
                        use_intel: bool = True):
        """Pre-fill a hopper with relaxed molecular chains."""
        if source_dir is None:
            source_dir = f"chain_data/relaxed/N{N}"
            
        if seed is None:
            seed = int(time.time()) % 1000000

        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        manager = HopperManager(runner)
        manager.generate_filled_state(source_dir=source_dir, 
                                     fill_template=fill_template,
                                     n_fill=n_fill, 
                                     dt=dt, 
                                     relax_steps=relax_steps,
                                     run_name=run_name,
                                     seed=seed,
                                     mol_dir=mol_dir, 
                                     setup_inc=setup_inc,
                                     dump_inc=dump_inc,
                                     viscosity=viscosity,
                                     N=N,
                                     outdir=outdir,
                                     num_procs=num_procs,
                                     num_threads=num_threads,
                                     use_kokkos=use_kokkos,
                                     use_intel=use_intel)

    def resume_hopper_fill(self, 
                          restart_path: str = None, 
                          source_dir: str = "chain_data/relaxed/N4", 
                          fill_template: str = "in.hopper_fill_resume",
                          n_fill: int = 10, 
                          relax_steps: int = 100000, 
                          run_name: Optional[str] = None, 
                          seed: Optional[int] = None, 
                          dt: float = 1e-6,
                          mol_dir: str = "chain_data/molecules_temp",
                          setup_inc: str = "simulation_geometries/2D_hopper.inc",
                          dump_inc: str = "simulation_templates/default_dump.inc",
                          viscosity: float = 0.001, 
                          N: int = 4,
                          outdir: Optional[str] = None,
                          num_procs: int = None,
                          num_threads: int = 1,
                          use_kokkos: bool = True,
                          use_intel: bool = True):
        """Resume filling a hopper from a binary restart file."""
        if seed is None:
            seed = int(time.time()) % 1000000

        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        manager = HopperManager(runner)
        manager.resume_filled_state(restart_path=restart_path,
                                    fill_template=fill_template,
                                    source_dir=source_dir,
                                    n_fill=n_fill,
                                    relax_steps=relax_steps,
                                    run_name=run_name,
                                    dt=dt,
                                    seed=seed,
                                    mol_dir=mol_dir, 
                                    setup_inc=setup_inc,
                                    dump_inc=dump_inc,
                                    viscosity=viscosity,
                                    N=N,
                                    outdir=outdir,
                                    num_procs=num_procs,
                                    num_threads=num_threads,
                                    use_kokkos=use_kokkos,
                                    use_intel=use_intel)

    def run_flop_simulation(self, N: int = 4, run_steps: int = 50000, viscosity: float = 0.001, dt: float = 1e-6, num_procs: int = None, num_threads: int = 1, use_kokkos: bool = True, use_intel: bool = True):
        """Run a single chain-flop simulation to analyze mobility."""
        viscosity_token = get_viscosity_token(viscosity)
        dt_token = get_dt_token(dt)

        config = SimulationConfig(
            template="in.chain_flop_template",
            data_file=f"chains_linear_x/N{N}_chain_horz.data",
            lepton_file="simulation_templates/lepton.inc",
            simulation="Chain_flop",
            run=f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}",
            extra_vars={
                "viscosity": viscosity,
                "run_steps": run_steps,
                "dt": dt
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )
        
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        print(f"Running simulation: {config.simulation}=>{config.run}")
        runner.run(config)

    def resume_flop_simulation(self, N: int = 4, run_steps: int = 50000, viscosity: float = 0.001, dt: float = 1e-6, resume_token: str = "100000", num_procs: int = None, num_threads: int = 1, use_kokkos: bool = True, use_intel: bool = True):
        """Resume a chain-flop simulation from a restart point."""
        viscosity_token = get_viscosity_token(viscosity)
        dt_token = get_dt_token(dt)

        config = SimulationConfig(**{
                "template": "in.chain_flop_resume",
                "resume_file": f"./dumping_yard/Chain_flop/N{N}_Viscosity_{viscosity_token}_dt_{dt_token}/restart/restart.{resume_token}.bin",
                "simulation": "Chain_flop",
                "run": f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}",
                "extra_vars":{
                    "viscosity": viscosity,
                    "run_steps": run_steps,
                    "dt": dt
                },
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )
        
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        print(f"Resuming simulation: {config.simulation}=>{config.run}")
        runner.resume(config)

    def generate_relaxed_library(self, n_beads: Any = 4, n_states: int = 10, forced: bool = False, 
                                 output_dir: str = "chain_data/relaxed",
                                 template: str = "in.relax_3d_gen",
                                 run_name_prefix: Optional[str] = None,
                                 simulation_name: str = "Relax_Library_Gen",
                                 dump_inc: str = "simulation_templates/default_dump.inc", 
                                 n_parallel: int = 1, num_procs: int = None, num_threads: int = 1, 
                                 use_kokkos: bool = True, use_intel: bool = True):
        """Generate a library of relaxed chain states for future hopper insertions."""
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        lib_gen = LibraryGenerator(runner, forced)
        
        # Handle single N vs multiple Ns
        if isinstance(n_beads, str):
            n_beads_list = [int(n.strip()) for n in n_beads.split(",")]
        elif isinstance(n_beads, list):
            n_beads_list = n_beads
        else:
            n_beads_list = [int(n_beads)]

        lib_gen.generate_batch_library(
            n_beads_list=n_beads_list, 
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

    def generate_chains(self, Ns: str = "4,6,8", orientation: str = "x", output_dir_name: str = "linear_x"):
        """Generate initial linear chain data files."""
        # Convert comma-separated string to list of ints for CLI/UI convenience
        if isinstance(Ns, str):
            N_list = [int(n.strip()) for n in Ns.split(",")]
        else:
            N_list = Ns if isinstance(Ns, list) else [Ns]

        for N in N_list:
            print(f"Generating linear chain with {N} beads, orientation={orientation}...")
            linear_config = ChainConfig(
                beads=N,
                spacing=0.0025,
                mode="linear",
                orientation=orientation,
                output_dir=Path(f"chain_data/{output_dir_name}"),
            )
            path = write_chain_data(linear_config)
            print(f"Created: {path}")

    def run_flop_batch(self, Ns: List[int] = [4, 6], run_steps: List[int] = [50000], viscosities: List[float] = [0.001], dt: float = 1e-6, num_procs: int = None, num_threads: int = 1, use_kokkos: bool = True, use_intel: bool = True):
        """Run a batch of flop simulations across multiple parameters."""
        total = len(Ns) * len(run_steps) * len(viscosities)
        eta = ETAEstimator(total=total)
        eta.start()
        completed = 0

        for N in Ns:
            for run_step in run_steps:
                for viscosity in viscosities:
                    completed += 1
                    viscosity_token = get_viscosity_token(viscosity)
                    dt_token = get_dt_token(dt)
                    config = SimulationConfig(**{
                        "template": "in.chain_flop_template",
                        "data_file": f"chains_linear_x/N{N}_chain_horz.data",
                        "lepton_file": "simulation_templates/lepton.inc",
                        "simulation": "Chain_flop",
                        "run": f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}",
                        "extra_vars":{
                            "viscosity": viscosity,
                            "run_steps": run_step,
                            "dt": dt
                        },
                    }, 
                    num_procs=num_procs,
                    num_threads=num_threads,
                    use_kokkos=use_kokkos,
                    use_intel=use_intel)
                    print(f"[{completed}/{total}] Running {config.run}...")
                    runner = SimulationRunner(lammps_executable=self.lammps_executable)
                    runner.run(config, verbose=False, clean_dir=True)
                    eta.update(completed)

    def run_grid_batch_relaxation(self, n_beads: Any = "4,8,12,24,48,100", n_states: int = 50,
                                  output_dir: str = "chain_data/relaxed_grid",
                                  spacing: float = 0.5,
                                  seed: Optional[int] = None,
                                  dt: float = 1e-6,
                                  num_procs: int = 1,
                                  num_threads: int = 1,
                                  use_kokkos: bool = True,
                                  use_intel: bool = True,
                                  lepton_file: str = "simulation_templates/lepton.inc",
                                  dump_file: str = "simulation_templates/quiet_dump.inc",
                                  simulation: str = "Grid_Relaxation_Batch",
                                  motion_steps: int = 500000,
                                  explore_steps: int = 200000,
                                  viscous_relax_steps: int = 300000,
                                  temperature: float = 1e15):
        """
        Runs a massive batch relaxation using a single grid-based simulation.
        This is the most performant way to generate large libraries using GPUs.
        """
        from simulation.grid_manager import GridRelaxManager
        from simulation.runner import SimulationRunner
        from simulation.config import SimulationConfig
        
        if seed is None:
            seed = int(time.time()) % 1000000
            
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        grid_manager = GridRelaxManager(runner)
        
        # 1. Parse N list
        if isinstance(n_beads, str):
            n_beads_list = [int(n.strip()) for n in n_beads.split(",")]
        elif isinstance(n_beads, list):
            n_beads_list = n_beads
        else:
            n_beads_list = [int(n_beads)]
            
        # 2. Generate Grid Data
        print(f"--- Preparing Grid Data for {len(n_beads_list)} N-values x {n_states} states ---")
        grid_data_path, metadata = grid_manager.generate_grid_data(n_beads_list, n_states, spacing)
        
        # 3. Configure Super-Simulation
        # We'll use a specific run name in dumping yard
        run_name = f"Grid_Relax_Batch_{len(n_beads_list)}N_{n_states}S"
        
        config = SimulationConfig(
            template="in.grid_relax",
            data_file=str(grid_data_path),
            lepton_file=lepton_file,
            dump_file=dump_file, 
            simulation=simulation,
            run=run_name,
            extra_vars={
                "seed": seed,
                "motion_steps": motion_steps, 
                "explore_steps": explore_steps, 
                "viscous_relax_steps": viscous_relax_steps, 
                "dt": dt,
                "temperature": temperature
            },
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel
        )
        
        # 4. Run the Super-Simulation
        print(f"--- Starting Super-Simulation: {run_name} ---")
        runner.run(config)
        
        # 5. Split Results
        final_grid_data = Path(config.output_dir) / "relaxed_grid.data"
        if final_grid_data.exists():
            grid_manager.split_grid_results(final_grid_data, metadata, output_dir)
        else:
            print(f"Error: Final grid data not found at {final_grid_data}")

    def run_grid_hopper_filling(self, n_hoppers: int = 4, n_fill: int = 10, N: int = 4,
                                spacing: float = 1.0, relax_steps: int = 500000,
                                seed: Optional[int] = None, dt: float = 1e-6,
                                output_dir: str = "chain_data/grid_filled",
                                source_dir: Optional[str] = None,
                                hopper_template_data: str = "simulation_geometries/2D_hopper.inc",
                                lepton_file: str = "simulation_templates/lepton.inc",
                                dump_file: str = "simulation_templates/quiet_dump.inc",
                                viscosity: float = 0.001,
                                num_procs: int = 1, num_threads: int = 1,
                                use_kokkos: bool = True,
                                mode: str = "2D_stacked",
                                simulation: str = "Grid_Hopper_Filling",
                                template: str = "in.grid_hopper_fill"):
        """
        Main entry point for grid-based batch hopper filling.
        Packs multiple hoppers into one simulation box for faster generation.
        """
        from simulation.grid_hopper_manager import GridHopperManager
        from simulation.runner import SimulationRunner
        
        if seed is None:
            seed = int(time.time()) % 1000000
            
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        grid_manager = GridHopperManager(runner)
        
        return grid_manager.run_grid_filling(
            n_hoppers=n_hoppers,
            n_fill_per_hopper=n_fill,
            N=N,
            spacing=spacing,
            relax_steps=relax_steps,
            seed=seed,
            dt=dt,
            output_dir=output_dir,
            source_dir=source_dir,
            hopper_template_data=hopper_template_data,
            lepton_file=lepton_file,
            dump_file=dump_file,
            viscosity=viscosity,
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            mode=mode,
            simulation=simulation,
            template=template
        )

    def resume_grid_hopper_filling(self, restart_path: str, relax_steps: int = 500000,
                                  dt: float = 1e-6, output_dir: str = "chain_data/grid_filled",
                                  lepton_file: str = "simulation_templates/lepton.inc",
                                  dump_file: str = "simulation_templates/quiet_dump.inc",
                                  viscosity: float = 0.001,
                                  num_procs: int = 1, num_threads: int = 1,
                                  use_kokkos: bool = True,
                                  template: str = "in.grid_hopper_fill_resume"):
        """
        Resumes a grid hopper filling simulation from a restart file.
        """
        from simulation.grid_hopper_manager import GridHopperManager
        from simulation.runner import SimulationRunner
        
        runner = SimulationRunner(lammps_executable=self.lammps_executable)
        grid_manager = GridHopperManager(runner)
        
        grid_manager.resume_grid_filling(
            restart_path=restart_path,
            relax_steps=relax_steps,
            dt=dt,
            output_dir=output_dir,
            lepton_file=lepton_file,
            dump_file=dump_file,
            viscosity=viscosity,
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            template=template
        )

    def run_hopper_flow(self, N: int = 4, run_steps: int = 200000, orifice_width: float = 0.05, dt: float = 1e-6):
        """(Upcoming) Run a hopper discharge/flow simulation."""
        print(f"Hopper Flow Simulation (N={N}) - Backend logic coming soon!")
        # This will be implemented as a separate template logic
