# build_lammps.ps1 - Native Windows Build Script for LAMMPS
# Optimizes for NVIDIA GTX 1060 (Pascal)

$ErrorActionPreference = "Stop"

$PROJECT_ROOT = Get-Location
$LAMMPS_DIR = Join-Path $PROJECT_ROOT "lammps"
$BUILD_DIR = Join-Path $LAMMPS_DIR "build_win"

Write-Host "==============================================================================" -ForegroundColor Cyan
Write-Host "⚡ Starting LAMMPS Build Process (Windows)"
Write-Host "Target: NVIDIA GTX 1060 (Pascal61) + MS-MPI"
Write-Host "==============================================================================" -ForegroundColor Cyan

# 1. Ensure LAMMPS source is present and updated
if (!(Test-Path (Join-Path $LAMMPS_DIR "cmake\CMakeLists.txt"))) {
    Write-Host "Step 0/3: LAMMPS source not found in $LAMMPS_DIR. Initializing..." -ForegroundColor Cyan
    
    if (!(Test-Path $LAMMPS_DIR)) {
        New-Item -ItemType Directory -Path $LAMMPS_DIR -Force | Out-Null
    }
    
    Set-Location $LAMMPS_DIR
    # Initialize as a git repo if not one already
    if (!(Test-Path ".git")) {
        git init
        git remote add origin https://github.com/lammps/lammps.git
    }
    
    Write-Host "Fetching LAMMPS source (stable branch)..." -ForegroundColor Cyan
    git fetch origin stable
    git checkout -f stable
    Set-Location $PROJECT_ROOT
} else {
    Write-Host "Step 0/3: Updating LAMMPS repository..." -ForegroundColor Cyan
    Set-Location $LAMMPS_DIR
    try {
        git pull origin stable
    } catch {
        Write-Host "⚠️ Warning: git pull failed. Proceeding with existing code." -ForegroundColor Yellow
    }
    Set-Location $PROJECT_ROOT
}

# Final validity check
if (!(Test-Path (Join-Path $LAMMPS_DIR "cmake\CMakeLists.txt"))) {
    Write-Host "❌ Error: Could not retrieve LAMMPS source into $LAMMPS_DIR" -ForegroundColor Red
    exit 1
}

# 2. Clean and Create build directory
if (Test-Path $BUILD_DIR) {
    Write-Host "Step 0/3: Cleaning old build artifacts..."
    Remove-Item -Path $BUILD_DIR -Recurse -Force
}
New-Item -ItemType Directory -Path $BUILD_DIR -Force | Out-Null
Set-Location $BUILD_DIR

# 3. Configure with CMake
Write-Host "Step 1/3: Configuring with CMake..."

# Find and load VS environment variables (vcvars64.bat)
$VS_PATH = & "C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe" -latest -property installationPath
$VCVARS_PATH = Join-Path $VS_PATH "VC\Auxiliary\Build\vcvars64.bat"
if (Test-Path $VCVARS_PATH) {
    Write-Host "  -> Loading VS environment variables..."
    $env_vars = cmd /c "`"$VCVARS_PATH`" && set"
    foreach ($line in $env_vars) {
        if ($line -match "^(.*?)=(.*)$") {
            $name = $matches[1]
            $value = $matches[2]
            if ($name -ieq "Path") { $env:Path = $value }
            else { Set-Item -Path "env:$name" -Value $value }
        }
    }
}

# Add Git's usr/bin to path for 'patch' utility (needed for VORO++)
$GIT_PATH = "C:\Program Files\Git\usr\bin"
if (Test-Path $GIT_PATH) { $env:PATH = "$GIT_PATH;$env:PATH" }

# Find x64 clang-cl.exe (avoiding ARM versions)
$CLANG_EXE = Get-ChildItem -Path $VS_PATH -Filter "clang-cl.exe" -Recurse -ErrorAction SilentlyContinue | `
             Where-Object { $_.FullName -like "*\x64\*" -and $_.FullName -notlike "*\ARM*" } | `
             Select-Object -First 1 -ExpandProperty FullName

if ($CLANG_EXE) {
    Write-Host "  -> Found Clang-cl: $CLANG_EXE"
} else {
    # Fallback to just "clang-cl" if we can't find the path but environment is loaded
    $CLANG_EXE = "clang-cl"
}

# Note: Using PASCAL61 for GTX 1060 and Ninja with latest VS 18 ClangCL
cmake -G "Ninja" `
      -D CMAKE_CXX_COMPILER="$CLANG_EXE" `
      -D CMAKE_C_COMPILER="$CLANG_EXE" `
      -C ../cmake/presets/most.cmake ../cmake `
      -D PKG_KOKKOS=ON `
      -D PKG_PYTHON=ON `
      -D PKG_INTEL=OFF `
      -D PKG_VORONOI=OFF `
      -D Kokkos_ENABLE_CUDA=ON `
      -D Kokkos_ENABLE_OPENMP=ON `
      -D Kokkos_ARCH_PASCAL61=ON `
      -D BUILD_SHARED_LIBS=OFF `
      -D CMAKE_BUILD_TYPE=Release `
      -D LAMMPS_SIZES=smallbig

# 4. Build
Write-Host "Step 2/3: Compiling LAMMPS (Release)..." -ForegroundColor Cyan
cmake --build . --config Release --parallel 8

# 5. Verification and Deployment
# Ninja puts the binary directly in the build folder
$LMP_EXE = Join-Path $BUILD_DIR "lmp.exe"
if (Test-Path $LMP_EXE) {
    Write-Host "==============================================================================" -ForegroundColor Green
    Write-Host "✅ LAMMPS Build Successful!"
    Write-Host "Binary Location: $LMP_EXE"
    
    Write-Host "Step 3/3: Deploying binary to project environment..." -ForegroundColor Cyan
    
    # Copy to root
    $ROOT_EXE = Join-Path $PROJECT_ROOT "lmp.exe"
    Copy-Item -Path $LMP_EXE -Destination $ROOT_EXE -Force
    Write-Host "  -> Copied to: $ROOT_EXE"
    
    # Copy to .venv if exists
    $VENV_SCRIPTS = Join-Path $PROJECT_ROOT ".venv\Scripts"
    if (Test-Path $VENV_SCRIPTS) {
        Copy-Item -Path $LMP_EXE -Destination (Join-Path $VENV_SCRIPTS "lmp.exe") -Force
        Write-Host "  -> Copied to: $VENV_SCRIPTS\lmp.exe"
    }

    Write-Host "Verification:"
    & $ROOT_EXE -h | Select-String "MOLECULE|RIGID|KOKKOS|LEPTON"
    Write-Host "==============================================================================" -ForegroundColor Green
} else {
    Write-Host "❌ Build failed. lmp.exe not found." -ForegroundColor Red
}

Set-Location $PROJECT_ROOT
