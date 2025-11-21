from analysis.data_manager import SimulationData
from analysis.geometry import get_angle_series, get_xyz_series, get_distance_series
from analysis.plotting import plot_angle_evolution, plot_xyz_evolution, plot_distance_evolution

def main():
    # 1. Initialize Data Manager
    data_dir = "post_chain_flop/current"
    sim = SimulationData(data_dir)
    
    # 2. Load Data (Auto-caches)
    df = sim.load_data()
    
    if df.empty:
        print("No data found!")
        return

    print(f"Loaded simulation with {len(df)} records.")
    # print df head
    # print(df.head())
    
    # 3. Perform Analysis
    # # Example: Calculate angle between atoms 1, 2, and 3
    # print("Calculating angles...")
    # angle_series = get_angle_series(df, id1=1, id2=2, id3=3)

    # # 4. Visualize
    # print("Plotting...")
    # plot_angle_evolution(angle_series, title="Angle between Atoms 1-2-3")# filepath: d:\Chains_simulations\main_analysis.py

    # plot x, y, z over time of an atom
    # atom_id=4
    # xyz_data = get_xyz_series(df, atom_id=atom_id)
    # plot_xyz_evolution(xyz_data, title=f'Atom {atom_id} Position Evolution')
    
    # get distances between atom 1 and atom 2 over time
    distance_series_23 = get_distance_series(df, id1=2, id2=3)
    distance_series_34 = get_distance_series(df, id1=3, id2=4)
    plot_distance_evolution(distance_series_23, title='Distance between Atoms 2 and 3')
    plot_distance_evolution(distance_series_34, title='Distance between Atoms 3 and 4')
    
    
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")