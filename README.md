# Granular Simulation Framework for LAMMPS: Advanced DEM & Geometry Visualization

A high-performance Python-driven framework for simulating, analyzing, and visualizing complex granular assemblies (Hoppers, Direct Shear Tests, Shear Vanes) using LAMMPS. This project provides a robust workflow for material insertion, state resumption, and **sharp analytical geometry extraction**, serving as a powerful open-source alternative to commercial visualization features.

## 🚀 Key Features

### 1. Advanced Simulation Workflows (Hoppers, DST, Shear Vane)
- **Automated Material Assembly**: Insertion of molecular chains or granular particles into complex boundary conditions (`run_hopper_fill`, etc.).
- **Physical State Resumption**: Seamlessly continue simulations from previous restart files to extend relaxation or change loading conditions.
- **Molecular Library**: Pre-relaxed states for molecular chains (N4, N8, etc.) to ensure stable simulation starts.

### 2. Visualizing LAMMPS Regions (Analytical VTK Mesh Extraction)
- **Sharp Geometry Extraction**: The `analysis/geometry_extractor.py` module converts analytical LAMMPS regions (`block`, `plane`, `intersect`, `union`) into perfectly sharp VTK meshes.
- **High-Fidelity "Pro" Feature**: Achieves the same high-quality visual results as expensive commercial software (like OVITO Pro) using a specialized **Lattice-Based Sampling** method.
- **Customizable Resolution**: Easily adjust mesh resolution for complex `fix wall/gran/region` or `fix wall/region` setups.

### 3. Interactive Analysis & 3D Viewing
- **Interactive 3D Viewer**: The `analysis/ui_viewer.py` tool (powered by `vedo`) offers a lightweight interface for real-time inspection of simulation trajectories and static geometry.
- **Optimized Data Parsing**: Fast loading and caching of LAMMPS dump files, designed for massive datasets.
- **Physics-Aware**: Native support for granular mechanics models like `hertz/history` and region-based friction.

## 📁 Project Structure

```bash
Chains_simulations/
├── main.py                     # Entry point for Granular Workflows (Filling, Resumes)
├── simulation/                 # Simulation Orchestration
│   ├── hopper_manager.py       # Manages material insertion and multi-stage runs
│   ├── runner.py               # Robust LAMMPS execution and variable management
│   ├── chain_generator.py      # Molecule structure design
│   └── config.py               # Centralized simulation settings
├── analysis/                   # Analysis & Visualization Modules
│   ├── geometry_extractor.py   # Visualise LAMMPS regions by extracting sharp meshes
│   ├── ui_viewer.py            # High-fidelity interactive 3D viewer
│   ├── data_manager.py         # Performance-optimized dump file loader
│   └── geometry.py             # Calculations for angles, distances, and orbits
├── simulation_templates/       # Modular LAMMPS input scripts (.filling, .resume)
├── simulation_geometries/      # Complex geometry definitions (.inc files)
├── chain_data/                 # Pre-relaxed chain libraries (N4, N100, etc.)
└── out_mesh/                   # Extracted VTK meshes for visualization
```

## 🛠️ Installation & Getting Started

### Prerequisites
- **Python 3.10+**
- **LAMMPS**: Standard `lmp` executable must be in your system's PATH.
- **FFmpeg**: (Optional) For high-quality video exports using NVENC acceleration.

### Setup
```bash
pip install numpy pandas matplotlib vedo ovito
```

## 📖 Usage Examples

### 1. Visualising LAMMPS Geometry (VTK Export)
To extract a perfectly sharp 3D mesh from your LAMMPS region definitions:
```bash
python analysis/geometry_extractor.py simulation_geometries/2D_hopper.inc \
    --regions funnel_neg_vis hopper_pos_vis \
    --spacing 0.001
```

### 2. Running a Granular Filling Simulation
Insert 100 chains into your geometry using the optimized managers:
```bash
python main.py run_hopper_fill n_fill=100 run_name="DST_initial_test"
```

### 3. Resuming a Simulation (State Persistence)
Continue from a prior restart file to add more material or extend runs:
```bash
python main.py resume_hopper_fill \
    restart_path="path/to/restart.bin" \
    n_fill=50 \
    run_name="DST_resume"
```

## 📐 Why This Framework? (Comparison)

Unlike standard stochastic sampling which creates "fuzzy" surfaces, this framework's **Lattice Sampling** ensures:
- **Zero Noise**: Edges and planes are perfectly flat, reflecting the true analytical geometry of your LAMMPS regions.
- **Memory Efficient**: High-resolution meshes are generated without the need for millions of random points.
- **Open Source**: Full access to advanced visualization without requiring premium commercial licenses.

---
*Created for the LAMMPS and Granular Physics communities.*
