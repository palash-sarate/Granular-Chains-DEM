# Granular Chains DEM Simulation Framework

This project provides a comprehensive framework for simulating, analyzing, and visualizing granular chains using the Discrete Element Method (DEM) with LAMMPS. It automates the process of generating chain structures, running simulations, and producing detailed analysis and visualizations.

## Features

### 1. Chain Generation
Flexible tools to generate LAMMPS data files for various chain configurations:
- **Linear**: Straight chains aligned along X, Y, or Z axes.
- **Directional**: Chains aligned along any arbitrary 3D vector.
- **Random**: Random walk (polymer-like) chain generation.
- **Loop**: Loop/Helix configurations.

### 2. Simulation Management
- **Automated Execution**: Python wrappers to manage LAMMPS simulations.
- **Templating**: Uses LAMMPS input script templates for flexible simulation setups.
- **Resume Capability**: Easily resume simulations from restart files.
- **Directory Management**: Automatically organizes output files (dumps, restarts, logs).

### 3. Analysis & Visualization
- **Data Parsing**: Efficiently parses LAMMPS dump files into Pandas DataFrames with caching for speed.
- **Geometric Analysis**:
  - Calculate bond angles between triplets of atoms over time.
  - Track inter-atomic distances.
  - Track atom trajectories (XYZ evolution).
- **Plotting**: Generate static plots for angle evolution, distance changes, and chain geometry.
- **3D Animation**: Create high-quality 3D animations of the simulation using Matplotlib and FFmpeg.
  - Supports GPU acceleration (NVIDIA NVENC) for fast rendering.
  - Customizable views (Isometric, Top, Side, etc.).

## Project Structure

```
Chains_simulations/
├── main.py                     # Main entry point for running workflows
├── simulation/                 # Simulation setup and execution tools
│   ├── chain_generator.py      # Generates LAMMPS data files
│   ├── runner.py               # Handles LAMMPS execution
│   └── config.py               # Configuration classes
├── analysis/                   # Analysis and visualization tools
│   ├── data_manager.py         # Loads and caches simulation data
│   ├── geometry.py             # Geometric calculations (angles, distances)
│   ├── plotting.py             # Static plotting functions
│   └── animate.py              # 3D animation generation
├── simulation_templates/       # LAMMPS input script templates
│   └── in.chain_flop_template  # Template for chain flop simulation
├── chain_data/                 # Generated chain structure files
└── dumping_yard/               # Simulation outputs (dumps, logs, plots)
```

## Prerequisites

1.  **Python 3.10+**
2.  **LAMMPS**: The LAMMPS executable (`lmp`) must be installed and accessible in your system's PATH.
3.  **FFmpeg**: Required for generating video animations. Must be in your system's PATH.

### Python Dependencies

Install the required Python packages:

```bash
pip install numpy pandas matplotlib
```

## Usage

### 1. Running the Main Workflow

The `main.py` script is the primary entry point. It can be configured to run simulations, generate chains, or visualize existing results.

```bash
python main.py
```

Modify the `main()` function in `main.py` to select the desired operation:

```python
def main():
    # 1. Run Simulations
    # run_flop_simulations()

    # 2. Visualize Results (Analysis & Animation)
    visualise_results("Chain_flop", "N4_Viscosity_03_dt_1e6")
    
    # 3. Generate New Chains
    # Ns = [4, 6, 8, 10]
    # generate_linear_chains(Ns, orientation="horz", output_dir="chains_linear_x")
```

### 2. Generating Chains Manually

You can use the chain generator as a standalone tool:

```bash
# Generate a linear chain with 50 beads
python simulation/chain_generator.py 50 --spacing 0.0025 --orientation vert

# Generate a random walk chain
python simulation/chain_generator.py 100 --mode random --spacing 0.0025
```

### 3. Visualizing Chain Data

Quickly visualize a generated LAMMPS data file:

```bash
python analysis/plotting.py chain_data/chains_linear_x/N4_chain_horz.data
```

## Outputs

Simulation results are stored in the `dumping_yard/` directory, organized by simulation name and run ID.

Example structure:
```
dumping_yard/Chain_flop/N4_Viscosity_03_dt_1e6/
├── chain/              # LAMMPS dump files
├── restart/            # Restart files
├── Angle_evol.png      # Angle evolution plot
├── Distance_evol.png   # Distance evolution plot
└── chain_motion.mp4    # 3D Animation of the simulation
```

## Customization

-   **Simulation Parameters**: Edit `SimulationConfig` in `main.py` to change viscosity, timestep (`dt`), or run duration.
-   **LAMMPS Script**: Modify `simulation_templates/in.chain_flop_template` to change the physics or boundary conditions of the simulation.
