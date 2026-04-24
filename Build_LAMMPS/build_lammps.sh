#!/bin/bash

# build_lammps.sh - Standardized build script for LAMMPS on master node
# Optimizes for AMD EPYC (Zen 3/4) and NVIDIA GPU (Ampere)

set -e

PROJECT_ROOT=$(pwd)
LAMMPS_DIR="$PROJECT_ROOT/lammps"
BUILD_DIR="$LAMMPS_DIR/build"
CONDA_ENV_PATH="/home/guest/miniconda3/envs/gchain"

# Ensure conda environment is active (now contains compatible GCC 11)
source /home/guest/miniconda3/etc/profile.d/conda.sh
conda activate gchain

echo "=============================================================================="
echo "⚡ Starting LAMMPS Build Process"
echo "Target: AMD EPYC (Zen 3/4) + NVIDIA GPU (Ampere)"
echo "Mode: Using GCC 11 (Conda) for optimal KOKKOS/CUDA compatibility"
echo "=============================================================================="

# 1. Check if lammps directory exists
if [ ! -d "$LAMMPS_DIR" ]; then
    echo "❌ Error: lammps directory not found in $PROJECT_ROOT"
    exit 1
fi

# 2. Clean and Create build directory
echo "Step 0/3: Cleaning old build artifacts..."
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# 3. Configure with CMake
echo "Step 1/3: Configuring with CMake (MOST Preset + Stable Performance Flags)..."
# We use the compilers in the environment (GCC 11)
cmake -C ../cmake/presets/most.cmake ../cmake \
      -D PKG_KOKKOS=ON \
      -D PKG_PYTHON=ON \
      -D PKG_INTEL=OFF \
      -D Kokkos_ENABLE_CUDA=ON \
      -D Kokkos_ENABLE_OPENMP=ON \
      -D Kokkos_ARCH_AMPERE86=ON \
      -D Kokkos_ARCH_ZEN3=ON \
      -D CMAKE_CUDA_FLAGS="-allow-unsupported-compiler"

# 4. Build with multiple cores
CORES=64
echo "Step 2/3: Compiling with $CORES cores (this may take a few minutes)..."
make -j$CORES

# 5. Deploy to Conda environment
echo "Step 3/3: Deploying binary to $CONDA_ENV_PATH/bin/lmp..."
cp lmp "$CONDA_ENV_PATH/bin/lmp"

echo "=============================================================================="
echo "✅ LAMMPS Build and Deployment Successful!"
echo "Verification:"
"$CONDA_ENV_PATH/bin/lmp" -h | grep -E "MOLECULE|RIGID|KOKKOS|LEPTON"
echo "=============================================================================="
