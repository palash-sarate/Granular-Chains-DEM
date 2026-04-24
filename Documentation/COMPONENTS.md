# Technical Documentation: Framework Components

This document provides a detailed overview of the modules and tools implemented in the Granular Simulation Framework.

## 1. Simulation Orchestration

### `simulation/runner.py` (SimulationRunner)
The core interface for LAMMPS execution.
- **Variable Injection**: Automatically manages LAMMPS variables (lepton, dump, seeds, etc.) by injecting header blocks into input scripts.
- **Hardware Acceleration**: Detects system capabilities (KOKKOS, INTEL packages, NVIDIA GPUs) and automatically applies optimal acceleration flags.
- **Process Management**: Handles MPI execution and OpenMP thread counts for cross-platform portability.

### `simulation/library_generator.py` (LibraryGenerator)
Automates the creation of pre-relaxed state libraries for molecular chains.
- **Parallel Execution**: Utilizes `ThreadPoolExecutor` for concurrent simulation runs, optimizing throughput on high-core-count hardware.
- **State Management**: Handles random seeding, data file preparation, and standardized output naming for large-batch generation.

### `simulation/hopper_manager.py` (HopperManager)
Orchestrates multi-stage material insertion workflows.
- **Filling Logic**: Manages the iterative insertion of particles or chains into defined simulation geometries.
- **Resumption**: Supports continuing simulations from binary restart files while maintaining geometry and boundary conditions.

---

## 2. Analysis & Visualization

### `analysis/geometry_extractor.py`
Extracts high-fidelity meshes from analytical LAMMPS region definitions.
- **Lattice-Based Sampling**: Employs a deterministic grid-sampling method to ensure perfectly sharp edges and smooth planes.
- **VTK Export**: Generates industry-standard meshes compatible with ParaView, OVITO, and the internal UI viewer.

### `analysis/ui_viewer.py`
An interactive 3D inspection tool built on the `vedo` library.
- **Trajectory Visualization**: Supports real-time rendering of LAMMPS dump files and extracted VTK geometries.
- **Interactive UI**: Provides controls for animation playback, color mapping, and cross-sectional analysis.

---

## 3. High-Performance Computing (HPC) Management

### `Pulse/pulse.py` (CLI) & `Pulse/dashboard.py` (Web)
A comprehensive job monitoring and management suite for PBS-based clusters.
- **Active Queue Monitoring**: Real-time status tracking of running, queued, and held jobs.
- **Live Tracker**: Dedicated monitoring for the most recently submitted job, displaying real-time walltime and RAM usage.
- **Historical Analysis**: Persistent metadata storage (`pulse_metadata.json`) allowing for the tracing of completed jobs.
- **Bulk Scanning**: Tool for batch-discovery of historical jobs across large ID ranges, extracting peak performance metrics (Peak RAM, CPU utilization).

### `PBS_scripts/`
Standardized job submission templates for cluster environments.
- **Environment Verification**: Scripts for validating KOKKOS/GPU availability on compute nodes.
- **Batch Generation**: Automated loops for generating state libraries across multiple chain lengths (N=4, 12, 24, 48).
