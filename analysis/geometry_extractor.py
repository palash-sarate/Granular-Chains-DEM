#!/usr/bin/env python3

import os
import argparse
import subprocess
from pathlib import Path

from ovito.io import import_file, export_file
from ovito.modifiers import ConstructSurfaceModifier


# --------------------------------------------------
# Helper: create LAMMPS input script
# --------------------------------------------------
def generate_lammps_input(inc_file, dump_file, lattice_spacing,
                         region_bounds, extra_vars, target_regions):

    xlo, xhi, ylo, yhi, zlo, zhi = region_bounds

    var_lines = "\n".join([f"variable {k} equal {v}" for k, v in extra_vars.items()])

    script = f"""
units si
atom_style sphere
boundary f f f
newton off

{var_lines}

# Simulation box (large enough)
region simbox block {xlo} {xhi} {ylo} {yhi} {zlo} {zhi}
create_box 1 simbox

# Include geometry
include {inc_file}

# Fill with structured lattice for perfect edges
lattice sc {lattice_spacing}
"""
    for region in target_regions:
        script += f"create_atoms 1 region {region}\n"

    script += f"\nwrite_dump all custom {dump_file} id type x y z\n"
    return script

    return script


# --------------------------------------------------
# Main pipeline
# --------------------------------------------------
class GeometryExtractor:
    def __init__(self, lammps_cmd: str = "lmp"):
        self.lammps_cmd = lammps_cmd

    def extract(self, inc_file: Path, outdir: Path, spacing: float = 0.001, 
                radius: float = None, regions: list = None, auto_vis: bool = True, 
                combined: bool = True, bounds: list = None,
                num_procs: int = 1, num_threads: int = 1, use_kokkos: bool = False):
        
        inc_file = Path(inc_file).resolve()
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        
        if not bounds:
            # Default wide bounds for grid simulations
            bounds = [-1.0, 1.0, -1.0, 1.0, -0.1, 1.0]

        extra_vars = {}

        # Auto-detect regions if requested or if no regions provided
        if auto_vis or not regions:
            with open(inc_file, "r") as f:
                detected = [line.split()[1] for line in f if line.strip().startswith("region") and "_vis" in line]
            if detected:
                print(f"[INFO] Auto-detected {len(detected)} visualization regions.")
                if regions:
                    regions = list(set(regions + detected))
                else:
                    regions = detected
            else:
                if not regions:
                    regions = ["simbox"]

        vtk_files = []
        batches = [regions] if combined else [[r] for r in regions]

        for batch in batches:
            label = "combined" if combined else batch[0]
            dump_file = outdir / f"dump_{label}.lammpstrj"
            input_file = outdir / f"in.generate_{label}"

            # Generate LAMMPS input
            lmp_script = generate_lammps_input(
                inc_file=inc_file,
                dump_file=dump_file,
                lattice_spacing=spacing,
                region_bounds=bounds,
                extra_vars=extra_vars,
                target_regions=batch
            )

            with open(input_file, "w") as f:
                f.write(lmp_script)

            print(f"[{label}] Running LAMMPS sampling for {len(batch)} regions (nprocs={num_procs}, threads={num_threads})...")
            
            cmd = [self.lammps_cmd, "-in", str(input_file)]
            
            # 2. Prepare environment (consistent with runner.py)
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(num_threads)
            env["OMP_PROC_BIND"] = "spread"
            env["OMP_PLACES"] = "cores"

            # Apply KOKKOS/Acceleration if requested (exact match with runner.py logic)
            if use_kokkos:
                # We assume GPU availability matches main sim preference
                cmd.extend([
                    "-k", "on", "t", str(num_threads), 
                    "-sf", "kk", 
                    "-pk", "kokkos", "newton", "on", "neigh", "half"
                ])
            elif num_threads > 1:
                cmd.extend(["-sf", "omp", "-pk", "omp", str(num_threads)])

            # Wrap in mpiexec
            if num_procs > 1:
                import sys
                if sys.platform == "win32":
                    cmd = ["mpiexec", "-n", str(num_procs)] + cmd
                else:
                    # Linux/macOS bind-to logic
                    if num_threads > 1:
                        cmd = ["mpiexec", "-n", str(num_procs), "--map-by", f"socket:PE={num_threads}", "--bind-to", "core"] + cmd
                    else:
                        cmd = ["mpiexec", "-n", str(num_procs), "--bind-to", "core"] + cmd

            try:
                subprocess.run(cmd, check=True, capture_output=True, env=env)
            except subprocess.CalledProcessError as e:
                print(f"Error running LAMMPS: {e.stderr.decode()}")
                return []

            print(f"[{label}] Constructing surface with OVITO...")
            pipeline = import_file(str(dump_file))
            recon_radius = radius if radius else spacing * 1.2
            pipeline.modifiers.append(
                ConstructSurfaceModifier(
                    radius=recon_radius,
                    smoothing_level=0
                )
            )

            vtk_file = outdir / f"{label}_mesh.vtk"
            export_file(pipeline, str(vtk_file), format="vtk/trimesh", key="surface")
            vtk_files.append(vtk_file)
            
            # Cleanup temporary files
            if dump_file.exists(): dump_file.unlink()
            if input_file.exists(): input_file.unlink()

        return vtk_files

def main():
    parser = argparse.ArgumentParser(description="Convert LAMMPS .inc geometry to mesh using OVITO")
    parser.add_argument("--inc", type=str, required=True, help="Input .inc file")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory")
    parser.add_argument("--spacing", type=float, default=0.001, help="Lattice sampling spacing")
    parser.add_argument("--radius", type=float, default=None, help="Surface reconstruction radius")
    parser.add_argument("--regions", type=str, nargs="+", default=None, help="Specific regions to extract")
    parser.add_argument("--auto_vis", action="store_true", help="Auto-detect visualization regions")
    parser.add_argument("--combined", action="store_true", help="Combine all regions into one mesh")
    parser.add_argument("--bounds", type=float, nargs=6, default=None, help="Bounding box (xlo xhi ylo yhi zlo zhi)")
    parser.add_argument("--lammps_cmd", type=str, default="lmp", help="LAMMPS executable command")

    args = parser.parse_args()
    
    extractor = GeometryExtractor(lammps_cmd=args.lammps_cmd)
    vtk_files = extractor.extract(
        inc_file=args.inc,
        outdir=args.outdir,
        spacing=args.spacing,
        radius=args.radius,
        regions=args.regions,
        auto_vis=args.auto_vis,
        combined=args.combined,
        bounds=args.bounds
    )

    if not os.environ.get("DISPLAY"):
        print("[INFO] No DISPLAY detected. Skipping interactive visualization.")
        return

    try:
        import vedo
        print("[INFO] Launching vedo viewer... (Close the viewer window to exit)")
        
        colors = ["lightblue", "salmon", "lightgreen", "gold", "plum", "tomato", "cyan", "orchid", "khaki"]
        mesh_objects = []
        for i, v_file in enumerate(vtk_files):
            obj = vedo.load(str(v_file))
            mesh = obj.tomesh() if hasattr(obj, "tomesh") else obj
            color = colors[i % len(colors)]
            mesh.c(color).alpha(0.8)
            mesh_objects.append(mesh)
            
        title = f"{inc_file.name} - Regions: {', '.join(args.regions)}"
        vedo.show(mesh_objects, axes=1, bg="white", title=title)
    except ImportError:
        print("[INFO] 'vedo' module not found. Skipping UI visualization.")


if __name__ == "__main__":
    main()