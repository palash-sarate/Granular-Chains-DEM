from analysis.data_manager import SimulationData
from analysis.utilities import get_angle_series, get_xyz_series, get_distance_series
from analysis.utilities import plot_angle_evolution, plot_xyz_evolution, plot_distance_evolution
from analysis.utilities import get_dt_token, get_viscosity_token, ETAEstimator
# import matplotlib.pyplot as plt
from analysis.animate import Animator
from analysis.vtk_exporter import VTKExporter
import os
# import uuid
from simulation import SimulationConfig, SimulationRunner
from simulation.chain_generator import ChainConfig, write_chain_data
from simulation.library_generator import LibraryGenerator
from simulation.hopper_manager import HopperManager
from pathlib import Path
import argparse
import ast
import sys

def main():
    # session_id = str(uuid.uuid4())[:8]
    # run_flop_simulations()
    run_hopper_simulation()
    # generate_relaxed_chain_states(n_beads=4, n_states=5)
    # generate chain along x for N 4,6,8,10,12,14,16,24,48,100
    # Ns = [4,6,8,10,12,14,16,24,48,100]
    # generate_linear_chains(Ns, orientation="horz", output_dir="chains_linear_x")

def run_hopper_simulation():
    runner = SimulationRunner(lammps_executable="lmp")
    manager = HopperManager(runner)
    manager.generate_filled_state(source_dir="chain_data/relaxed/N4", 
                                n_fill=100, dt = 1e-6, relax_steps=100000,
                                run_name="hopper_fill_N4",
                                seed = 12345,
                                mol_dir = "chain_data/molecules_temp", 
                                setup_inc = "simulation_geometries/2D_hopper_flow_setup.inc")

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
  
def visualize_chain_flop_results(Ns, viscosity, dt):
    viscosity_token = get_viscosity_token(viscosity)
    dt_token = get_dt_token(dt)
    for N in Ns:
        visualize_results("Chain_flop", f"N{N}_Viscosity_{viscosity_token}_dt_{dt_token}")
   
def visualize_results(simulation, run):
    # 1. Initialize Data Manager
    data_dir = f"dumping_yard/{simulation}/{run}"
    save_dir = f"dumping_yard/{simulation}/{run}"
    
    # make save_dir if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    
    sim = SimulationData(data_dir)
    
    # 2. Load Data (Auto-caches)
    df = sim.load_data(force_reload=True)
    
    if df.empty:
        print("No data found!")
        return

    print(f"Loaded simulation with {len(df)} records.")

    # get the number of atoms from the dataframe
    num_atoms = df.index.get_level_values('id').nunique()
    print(f"Number of unique atoms: {num_atoms}")
    
    # 3. Perform Analysis
    # # Example: Calculate angle between atoms 1, 2, and 3
    # print("Calculating angles...")
    if num_atoms < 3:
        print("Need at least three atoms to compute angles.")
        return

    angle_series_list = []
    angle_labels = []

    for start_id in range(1, num_atoms - 1):
        id1, id2, id3 = start_id, start_id + 1, start_id + 2
        angle_series_list.append(get_angle_series(df, id1=id1, id2=id2, id3=id3))
        angle_labels.append(f"{id1}-{id2}-{id3}")

    plot_angle_evolution(angle_series_list, legend_labels=angle_labels, save_path=f"{save_dir}/Angle_evol.png")

    # plot x, y, z over time of an atom
    # atom_id=4
    # xyz_data = get_xyz_series(df, atom_id=atom_id)
    # plot_xyz_evolution(xyz_data, title=f'Atom {atom_id} Position Evolution')
    
    if num_atoms < 2:
        print("Need at least two atoms to compute distances.")
    else:
        distance_series_list = []
        distance_labels = []
        for atom_id in range(1, num_atoms):
            id1, id2 = atom_id, atom_id + 1
            distance_series_list.append(get_distance_series(df, id1=id1, id2=id2))
            distance_labels.append(f"{id1}-{id2}")
        plot_distance_evolution(
            distance_series_list,
            legend_labels=distance_labels,
            save_path=f"{save_dir}/Distance_evol.png",
        )
    # print("Generating animation...")
    
    # Find LAMMPS script in data_dir
    lammps_script = None
    if os.path.exists(data_dir):
        for file in os.listdir(data_dir):
            if file.startswith("in."):
                lammps_script = os.path.join(data_dir, file)
                break
    
    anim = Animator(df, output_file=f"{save_dir}/chain_motion.mp4", lammps_script=lammps_script)
    
    # You can color by 'id', 'vx', 'vy', 'vz', or 'velocity_magnitude' if those columns exist
    # Axis limits set to 15mm (0.015m) as requested
    anim.create_animation(start_frame=1, end_frame=None, fps=24, color_by='id', point_size=100, view='z_left_x_down', axis_limits=None)

def export_vtk(simulation, run):
    """
    Export simulation data to VTK format for ParaView.
    """
    data_dir = f"dumping_yard/{simulation}/{run}"
    save_dir = f"dumping_yard/{simulation}/{run}/vtk_output"
    
    print(f"Exporting VTK for {simulation}/{run}...")
    
    # Load Data
    sim = SimulationData(data_dir)
    df = sim.load_data(force_reload=False)
    
    if df.empty:
        print("No data found!")
        return

    # Find LAMMPS script
    lammps_script = None
    if os.path.exists(data_dir):
        for file in os.listdir(data_dir):
            if file.startswith("in."):
                lammps_script = os.path.join(data_dir, file)
                break
    
    exporter = VTKExporter(df, save_dir, lammps_script)
    exporter.export_series()
    print(f"Export complete. Files saved to {save_dir}")

if __name__ == "__main__":
    available_functions = {
        "main": main,
        "run_flop_simulations": run_flop_simulations,
        "run_flop_simulation": run_flop_simulation,
        "resume_flop_simulation": resume_flop_simulation,
        "run_simulation": run_simulation,
        "generate_linear_chains": generate_linear_chains,
        "visualize_results": visualize_results,
        "visualize_chain_flop_results": visualize_chain_flop_results,
        "generate_relaxed_chain_states": generate_relaxed_chain_states,
        "run_hopper_simulation": run_hopper_simulation,
        "export_vtk": export_vtk,
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
    for token in parsed.func_args:
        if "=" in token:
            key, value = token.split("=", 1)
            keyword[key] = _coerce(value)
        else:
            positional.append(_coerce(token))

    try:
        available_functions[parsed.func](*positional, **keyword)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(1)