from analysis.data_manager import SimulationData
from analysis.geometry import get_angle_series
from analysis.plotting import plot_angle_evolution

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

    # 3. Perform Analysis
    # Example: Calculate angle between atoms 1, 2, and 3
    print("Calculating angles...")
    angle_series = get_angle_series(df, id1=1, id2=2, id3=3)

    # 4. Visualize
    print("Plotting...")
    plot_angle_evolution(angle_series, title="Angle between Atoms 1-2-3")# filepath: d:\Chains_simulations\main_analysis.py

if __name__ == "__main__":
    main()