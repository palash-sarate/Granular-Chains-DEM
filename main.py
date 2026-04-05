from analysis.data_manager import SimulationData
from analysis.utilities import get_dt_token, get_viscosity_token, ETAEstimator
import os
from simulation import SimulationConfig, SimulationRunner
from simulation.chain_generator import ChainConfig, write_chain_data
from simulation.library_generator import LibraryGenerator
from simulation.hopper_manager import HopperManager
from pathlib import Path
import argparse
import ast
import sys
import time

def main():
    pass

def run_hopper_fill(source_dir="chain_data/relaxed/N4", n_fill=10, 
                       relax_steps=100000, run_name=None, seed=None, dt = 1e-6,
                       mol_dir = "chain_data/molecules_temp",
                       setup_inc = "simulation_geometries/2D_hopper.inc",
                       dump_inc = "simulation_templates/default_dump.inc"):
    if seed is None:
        seed = int(time.time()) % 1000000

    runner = SimulationRunner(lammps_executable="lmp")
    manager = HopperManager(runner)
    manager.generate_filled_state(source_dir=source_dir, 
                                 n_fill=n_fill, dt = dt, relax_steps=relax_steps,
                                 run_name=run_name,
                                 seed = seed,
                                 mol_dir = mol_dir, 
                                 setup_inc = setup_inc,
                                 dump_inc = dump_inc)

def resume_hopper_fill(source_dir="chain_data/relaxed/N4", n_fill=10, 
                       relax_steps=100000, restart_path=None, 
                       run_name=None, seed=None):
    """Resume filling a hopper from an existing restart file."""
    if not restart_path:
        print("Error: restart_path is required for resuming.")
        return

    if seed is None:
        seed = int(time.time()) % 1000000

    runner = SimulationRunner(lammps_executable="lmp")
    manager = HopperManager(runner)
    manager.resume_filled_state(restart_path=restart_path,
                                source_dir=source_dir,
                                n_fill=n_fill,
                                relax_steps=relax_steps,
                                run_name=run_name,
                                dt=1e-6,
                                seed = seed,
                                mol_dir = "chain_data/molecules_temp",
                                setup_inc="simulation_geometries/2D_hopper.inc",
                                dump_inc = "simulation_templates/default_dump.inc")

def generate_relaxed_chain_states(n_beads=4, n_states=10, forced = False):
    runner = SimulationRunner(lammps_executable="lmp")
    lib_gen = LibraryGenerator(runner, forced)
    lib_gen.generate_library(n_beads=n_beads, n_states=n_states)

def run_flop_simulations(Ns = [6], 
                         run_steps = [50000],
                         viscosities = [i/100000 for i in range(60, 101, 5)], dt = 1e-6):
    # Provide an on-the-fly ETA estimator for the nested loop
    total = len(Ns) * len(run_steps) * len(viscosities)
    eta = ETAEstimator(total=total)
    eta.start()
    completed = 0

    for N in Ns:
        for run_step in run_steps:
            for viscosity in viscosities:
                completed += 1
                try:
                    # Build config here so we can determine output dir and logfile
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
                    })

                    # Ensure log directory exists and create per-simulation logfile
                    log_path = f"{config.output_dir}/lammps.log"
                    # Print concise progress line with overall ETA and current simulation name
                    overall_line = f"[{completed}/{total}] {int(100*eta.progress_fraction() if eta.progress_fraction() is not None else 0):3d}% | {eta}"
                    current_line = f"Running: {config.simulation} => {config.run} | log: {log_path}"
                    # Overwrite two terminal lines
                    sys.stdout.write('\r' + ' ' * 120 + '\r')
                    sys.stdout.write(overall_line + '\n' + current_line + '\n')
                    sys.stdout.flush()

                    # Run simulation with output sent to logfile (no verbose terminal output)
                    runner = SimulationRunner(lammps_executable="lmp")
                    runner.run(config, verbose=False, clean_dir=True)
                except Exception as e:
                    # Log the exception to the simulation log if possible, otherwise print
                    try:
                        with open(log_path, 'a') as fh:
                            fh.write(f"\nERROR: {e}\n")
                    except Exception:
                        print(f"Error processing N={N}: {e}")
                    print("Continuing with the generated files...")
                finally:
                    eta.update(completed)
                    # move cursor up to keep only latest two lines visible
                    sys.stdout.write('\x1b[2A')
                    sys.stdout.flush()
                    visualize_results("Chain_flop", f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}")
                
def run_flop_simulation(N, run_step, viscosity, dt):
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
            # "regions": [
            #     "region hopper_cone cone z 0.0 0.0 0.01 0.05 0.05 0.15 open 1 open 2 move v_x_osc v_zero v_zero",
            #     "region hopper_cyl cylinder z 0.0 0.0 0.01 0.02 0.05 open 1 open 2 move v_x_osc v_zero v_zero",
            #     "region hopper_union union 2 hopper_cone hopper_cyl",
            #     "region floor_plane plane 0 0 -0.45 0 0 1"
            # ],
            # "wall_blocks": [
            #     "fix wall all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region hopper_union",
            #     "fix floor all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region floor_plane"
            # ]
        })
    
    try:
        run_simulation(config)
    except Exception as e:
        print(f"Error running sim N={N}: {e}")
    
    visualize_results("Chain_flop", f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}")

def resume_flop_simulation(N, run_step, viscosity, dt, resume_token):
    viscosity_token = get_viscosity_token(viscosity)
    dt_token = get_dt_token(dt)

    config = SimulationConfig(**{
            "template": "in.chain_flop_resume",
            "resume_file": f"./dumping_yard/Chain_flop/N{N}_Viscosity_{viscosity_token}_dt_{dt_token}/restart/restart.{resume_token}.bin",
            "simulation": "Chain_flop",
            "run": f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}",
            "extra_vars":{
                "viscosity": viscosity,
                "run_steps": run_step,
                "dt": dt
            },
            # regions and walls should be the same as the original template, so no need to redefine them here
        })
    
    try:
        resume_simulation(config)
    except Exception as e:
        print(f"Error Resuming sim N={N}: {e}")
        
    visualize_results("Chain_flop", f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}")

def create_hopper_filled_state(chain_source_dir, n_fill, relax_steps, dt, run_name):
    runner = SimulationRunner(lammps_executable="lmp")
    manager = HopperManager(runner)
    return manager.generate_filled_state(chain_source_dir, n_fill, relax_steps, dt=dt, run_name=run_name)

def resume_simulation(config: SimulationConfig):
    # Initialize runner
    # Ensure 'lmp' is in your PATH or provide absolute path
    runner = SimulationRunner(lammps_executable="lmp")
    
    print(f"Resuming simulation: {config.simulation}=>{config.run}")
    print(f"Output directory: {config.output_dir}")
        
    runner.resume(config)

def run_simulation(config: SimulationConfig):
    # Initialize runner
    # Ensure 'lmp' is in your PATH or provide absolute path
    runner = SimulationRunner(lammps_executable="lmp")

    print(f"Running simulation: {config.simulation}=>{config.run}")
    print(f"Output directory: {config.output_dir}")

    # Run simulation
    # This will create directories: post_chain_flop/Viscosity_03/{bond,angle,restart}
    # and generate a temporary input script 'generated_in.chain_flop'
    runner.run(config)
    
def generate_linear_chains(Ns: list[int], orientation: str, output_dir: str) -> None:
    for N in Ns:
        print(f"Generating linear chain with {N} beads, orientation={orientation}...")
        linear_config = ChainConfig(
            beads=N,
            # center to center spacing 0 to 0.0025 max
            spacing=0.0025,
            mode="linear",
            orientation=orientation,
            output_dir=Path(f"chain_data/{output_dir}"),
        )
        path = write_chain_data(linear_config)
        print(f"Created: {path}")
  
if __name__ == "__main__":
    available_functions = {
        "main": main,
        "run_flop_simulations": run_flop_simulations,
        "run_flop_simulation": run_flop_simulation,
        "resume_flop_simulation": resume_flop_simulation,
        "run_simulation": run_simulation,
        "generate_linear_chains": generate_linear_chains,
        "generate_relaxed_chain_states": generate_relaxed_chain_states,
        "run_hopper_fill": run_hopper_fill,
        "resume_hopper_fill": resume_hopper_fill
    }

    parser = argparse.ArgumentParser(description="Execute functions from main.py")
    parser.add_argument("func", nargs="?", default="main", choices=available_functions.keys())
    parser.add_argument("func_args", nargs="*", help="Positional or key=value arguments")
    parsed = parser.parse_args()

    def _coerce(token: str):
        try:
            return ast.literal_eval(token)
        except (ValueError, SyntaxError):
            return token

    positional, keyword = [], {}
    i = 0
    while i < len(parsed.func_args):
        token = parsed.func_args[i]
        
        # Handle cases like "key=value", "key = value", "key= value", "key =value"
        if "=" in token:
            key, value = token.split("=", 1)
            # If value is empty, the actual value might be the next token
            if not value and i + 1 < len(parsed.func_args):
                value = parsed.func_args[i+1]
                i += 1
            keyword[key.strip()] = _coerce(value.strip())
        elif i + 1 < len(parsed.func_args) and parsed.func_args[i+1] == "=":
            # Handle "key = value"
            key = token
            if i + 2 < len(parsed.func_args):
                value = parsed.func_args[i+2]
                i += 2
            else:
                value = "" # No value provided
                i += 1
            keyword[key.strip()] = _coerce(value.strip())
        else:
            positional.append(_coerce(token))
        i += 1

    try:
        available_functions[parsed.func](*positional, **keyword)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(1)
    except TypeError as e:
        print(f"Error calling '{parsed.func}': {e}")
        print(f"Arguments provided: positional={positional}, keyword={keyword}")
        sys.exit(1)