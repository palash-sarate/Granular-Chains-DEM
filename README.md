# Granular Simulation Framework for LAMMPS

A comprehensive framework for the simulation, analysis, and management of granular and molecular chain assemblies using LAMMPS. This repository provides automated workflows for material insertion, high-fidelity geometry extraction, and cluster job management.

## System Architecture

The framework is organized into three primary functional layers:

### 1. Simulation Orchestration
- **Automated Generation**: Facilities for generating pre-relaxed molecular state libraries (`LibraryGenerator`) with support for multi-threaded parallel execution.
- **Workflow Management**: Management of material assembly in complex geometries, including hoppers and direct shear tests.
- **Hardware Portability**: Automatic detection and utilization of KOKKOS and INTEL acceleration packages, as well as GPU offloading.

### 2. Analysis & Visualization
- **Geometry Extraction**: A specialized module for converting analytical LAMMPS regions into sharp VTK meshes using deterministic lattice-based sampling.
- **Interactive Inspection**: A 3D viewer for the visualization of simulation trajectories and extracted geometries.
- **Performance Optimization**: Fast parsing and caching of large-scale LAMMPS datasets.

### 3. Cluster Management (Pulse)
- **Monitoring**: A real-time dashboard for PBS-based HPC clusters.
- **Tracking**: Live monitoring of simulation resource utilization (RAM, CPU, Walltime).
- **History Persistence**: Persistent tracking of historical job performance via a local metadata cache and automated log tracing.

## Project Structure

```text
Granular-Chains-DEM/
├── main.py                     # CLI entry point for simulation workflows
├── simulation/                 # Simulation logic and orchestration
│   ├── library_generator.py    # Parallel generation of relaxed states
│   ├── runner.py               # LAMMPS execution and hardware detection
│   └── hopper_manager.py       # Material insertion management
├── analysis/                   # Post-processing and visualization
│   ├── geometry_extractor.py   # VTK mesh extraction from regions
│   └── ui_viewer.py            # Interactive 3D visualization
├── Pulse/                      # HPC Job Management Suite
│   ├── pulse.py                # Pulse CLI
│   ├── dashboard.py            # Pulse Streamlit Dashboard
│   └── pulse_core.py           # PBS backend and metadata persistence
├── PBS_scripts/                # HPC job submission templates
├── simulation_templates/       # Modular LAMMPS input scripts
└── Documentation/              # Detailed technical documentation
```

## Technical Specifications

### Prerequisites
- **Python**: Version 3.10 or higher.
- **LAMMPS**: Executable with MOLECULE and LEPTON packages (KOKKOS recommended).
- **Environment**: Linux-based HPC cluster (PBS) or local workstation.

### Core Modules Reference
For detailed technical documentation on specific components, please refer to [COMPONENTS.md](Documentation/COMPONENTS.md).

## Usage Methodology

### Parallel Library Generation
Generate a library of relaxed states utilizing multiple CPU cores:
```bash
python main.py generate_relaxed_library --n_beads 4 --n_states 100 --inParallel 16 --num_procs 4
```

### HPC Job Monitoring
Launch the Pulse dashboard for real-time cluster monitoring:
```bash
python Pulse/pulse.py dashboard
```

### Geometry Extraction
Export analytical regions to VTK format for visualization:
```bash
python analysis/geometry_extractor.py <geometry_file.inc> --regions <region_names> --spacing 0.001
```

---
*Technical framework developed for advanced granular mechanics and molecular dynamics research.*
