from analysis.data_manager import SimulationData
from analysis.geometry import get_angle_series, get_xyz_series, get_distance_series
from analysis.plotting import plot_angle_evolution, plot_xyz_evolution, plot_distance_evolution
import matplotlib.pyplot as plt
from analysis.animate import Animator

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
    
    # # 3. Perform Analysis
    # # # Example: Calculate angle between atoms 1, 2, and 3
    # # print("Calculating angles...")
    # angle_123 = get_angle_series(df, id1=1, id2=2, id3=3)
    # angle_234 = get_angle_series(df, id1=2, id2=3, id3=4)

    # # # 4. Visualize
    # plot_angle_evolution([angle_123, angle_234], legend_labels=["1-2-3","2-3-4"])
    

    # # plot x, y, z over time of an atom
    # # atom_id=4
    # # xyz_data = get_xyz_series(df, atom_id=atom_id)
    # # plot_xyz_evolution(xyz_data, title=f'Atom {atom_id} Position Evolution')
    
    # # get distances between atom 1 and atom 2 over time
    # distance_series_12 = get_distance_series(df, id1=1, id2=2)
    # distance_series_23 = get_distance_series(df, id1=2, id2=3)
    # distance_series_34 = get_distance_series(df, id1=3, id2=4)
    # plot_distance_evolution([distance_series_12,distance_series_23,distance_series_34], legend_labels=['1-2','2-3','3-4'])

    # plt.show()
    
    print("Generating animation...")
    anim = Animator(df, output_file="chain_movement.mp4")
    
    # You can color by 'id', 'vx', 'vy', 'vz', or 'velocity_magnitude' if those columns exist
    anim.create_animation(start_frame=1, end_frame=200, fps=24, color_by='id', point_size=100)

    
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")