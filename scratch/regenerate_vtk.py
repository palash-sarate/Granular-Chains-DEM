
import sys
import os
from pathlib import Path

# Fix for Ovito DLL load error on Windows
conda_prefix = Path(r"d:\Chains_simulations\.conda")
library_bin = conda_prefix / "Library" / "bin"
if library_bin.exists():
    os.environ["PATH"] = str(library_bin) + os.pathsep + os.environ["PATH"]
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(library_bin))

# Add project root to sys.path to allow imports
project_root = Path(r"d:\Chains_simulations")
sys.path.append(str(project_root))

from analysis.geometry_extractor import GeometryExtractor

def regenerate_vtk():
    job_dir = project_root / "dumping_yard" / "Grid_Hopper_Filling" / "Grid_Fill_5H_MixedN_S466133"
    if not job_dir.exists():
        # Fallback to simulation type from user state if different
        job_dir = project_root / "dumping_yard" / "Grid_Hopper_Filling_test" / "Grid_Fill_5H_MixedN_S466133"
        
    if not job_dir.exists():
        print(f"Error: Job directory not found at {job_dir}")
        return

    inc_file = job_dir / "replicated_geometry.inc"
    out_dir = job_dir / "Geometry_vtk"
    
    print(f"Regenerating VTK files for {job_dir.name}...")
    
    # We need a lammps command. Assuming 'lmp' or looking it up from a known runner config if possible.
    # Usually it's 'lmp' in this environment or we can try to find it.
    lammps_cmd = "lmp" # Default
    
    extractor = GeometryExtractor(lammps_cmd=lammps_cmd)
    
    # Using 5 hoppers logic to determine bounds if needed
    # But extract() can detect regions and we can pass bounds.
    # From grid_hopper_manager.py:
    n_hoppers = 5
    spacing = 2.0 # Standard spacing? In metadata.json it was 2.0?
    # Let's check metadata.json first
    import json
    with open(job_dir / "metadata.json", 'r') as f:
        meta = json.load(f)
        spacing = meta.get("spacing", 2.0)
    
    import math
    cols = math.ceil(math.sqrt(n_hoppers))
    rows = math.ceil(n_hoppers / cols)
    max_x = cols * spacing
    max_y = rows * spacing
    vtk_bounds = [-spacing, max_x, -spacing, max_y, -0.1, 1.0]

    vtk_files = extractor.extract(
        inc_file=inc_file,
        outdir=out_dir,
        auto_vis=True,
        combined=False,
        bounds=vtk_bounds,
        spacing=0.0005 # The fix
    )
    
    print(f"Done! Generated {len(vtk_files)} VTK files in {out_dir}")

if __name__ == "__main__":
    regenerate_vtk()
