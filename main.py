from analysis.data_manager import SimulationData
from analysis.geometry import get_angle_series, get_xyz_series, get_distance_series
from analysis.plotting import plot_angle_evolution, plot_xyz_evolution, plot_distance_evolution
# import matplotlib.pyplot as plt
from analysis.animate import Animator
import os
# import uuid
from simulation import SimulationConfig, SimulationRunner
from simulation.chain_generator import ChainConfig, write_chain_data
from pathlib import Path

def main():
    # session_id = str(uuid.uuid4())[:8]
    run_flop_simulations()
    # visualise_results("Chain_flop", "N4_Viscosity_03_dt_1e6")
    # visualise_results("Chain_flop", "N6_Viscosity_03_dt_1e6")
    
    # generate chain along x for N 4,6,8,10,12,14,16,24,48,100
    # Ns = [4,6,8,10,12,14,16,24,48,100]
    # generate_linear_chains(Ns, orientation="horz", output_dir="chains_linear_x")
    
def run_flop_simulations():
    # Ns = [4,6,8,10,12,14,16,24,48,100]
    Ns = [4,6]
    for N in Ns:
        print(f"Running flop simulation for N={N}...")
        config = SimulationConfig(
            template = "in.chain_flop_template",
            data_file = f"chains_linear_x/N{N}_chain_horz.data",
            simulation = "Chain_flop",
            run = f"N{N}_Viscosity_003_dt_1e6",
            extra_vars={
                "viscosity": 0.003,
                "run_steps": 100000,
                "dt": 1e-6
            }
        )

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
    
def visualise_results(simulation, run):
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
    
    # get distances between atom 1 and atom 2 over time
    distance_series_12 = get_distance_series(df, id1=1, id2=2)
    distance_series_23 = get_distance_series(df, id1=2, id2=3)
    distance_series_34 = get_distance_series(df, id1=3, id2=4)
    plot_distance_evolution([distance_series_12,distance_series_23,distance_series_34], legend_labels=['1-2','2-3','3-4'], save_path=f"{save_dir}/Distance_evol.png")
    
    # print("Generating animation...")
    anim = Animator(df, output_file=f"{save_dir}/chain_motion.mp4")
    
    # You can color by 'id', 'vx', 'vy', 'vz', or 'velocity_magnitude' if those columns exist
    # Axis limits set to 15mm (0.015m) as requested
    anim.create_animation(start_frame=1, end_frame=None, fps=24, color_by='id', point_size=100, view='z_left_x_down', axis_limits=0.015)
    
def run_simulation():
    
    # Create a config with a specific run name
    config = SimulationConfig(
        template = "in.chain_flop_template",
        data_file = "N5_chain.data",
        simulation = "Chain_flop",
        run = "N4_Viscosity_03_dt_1e6",
        extra_vars={
            "viscosity": 0.03,
            "run_steps": 10000,
            "dt": 1e-6
        }
    )

    # Initialize runner
    # Ensure 'lmp' is in your PATH or provide absolute path
    runner = SimulationRunner(lammps_executable="lmp")

    print(f"Running simulation: {config.simulation}=>{config.run}")
    print(f"Output directory: {config.output_dir}")

    # Run simulation
    # This will create directories: post_chain_flop/Viscosity_03/{bond,angle,restart}
    # and generate a temporary input script 'generated_in.chain_flop'
    runner.run(config)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")