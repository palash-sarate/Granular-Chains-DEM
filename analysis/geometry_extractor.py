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
                         region_bounds, extra_vars, target_region):

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
create_atoms 1 region {target_region}

write_dump all custom {dump_file} id type x y z
"""

    return script


# --------------------------------------------------
# Main pipeline
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Convert LAMMPS .inc geometry to mesh using OVITO")

    parser.add_argument("inc", help="Path to .inc file")
    parser.add_argument("--outdir", default="out_mesh", help="Output directory")
    parser.add_argument("--spacing", type=float, default=0.001, help="Lattice spacing for sampling")
    parser.add_argument("--radius", type=float, default=None, help="Surface reconstruction radius (defaults to 1.2 * spacing)")
    parser.add_argument("--lammps_cmd", default="lmp", help="LAMMPS executable")

    parser.add_argument("--regions", nargs="+", default=["simbox"], help="List of LAMMPS region names to fill")
    parser.add_argument("--bounds", nargs=6, type=float,
                        default=[-0.2, 0.2, -0.3, 0.3, -0.1, 0.6],
                        help="Sampling box bounds: xlo xhi ylo yhi zlo zhi")

    parser.add_argument("--var", action="append",
                        help="Extra LAMMPS variables (format: name=value)")

    args = parser.parse_args()

    inc_file = Path(args.inc).resolve()
    outdir = Path(args.outdir) / inc_file.stem
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse extra variables
    extra_vars = {}
    if args.var:
        for v in args.var:
            k, val = v.split("=")
            extra_vars[k] = val

    vtk_files = []

    for region in args.regions:
        dump_file = outdir / f"dump_{region}.lammpstrj"
        input_file = outdir / f"in.generate_{region}"

        # Generate LAMMPS input
        lmp_script = generate_lammps_input(
            inc_file=inc_file,
            dump_file=dump_file,
            lattice_spacing=args.spacing,
            region_bounds=args.bounds,
            extra_vars=extra_vars,
            target_region=region
        )

        with open(input_file, "w") as f:
            f.write(lmp_script)

        print(f"[{region}] Running LAMMPS (Lattice fill)...")
        subprocess.run([args.lammps_cmd, "-in", str(input_file)], check=True)

        print(f"[{region}] Loading dump into OVITO...")
        pipeline = import_file(str(dump_file))

        print(f"[{region}] Constructing surface (Zero smoothing)...")
        recon_radius = args.radius if args.radius else args.spacing * 1.2
        pipeline.modifiers.append(
            ConstructSurfaceModifier(
                radius=recon_radius,
                smoothing_level=0
            )
        )

        vtk_file = outdir / f"{region}_mesh.vtk"

        print(f"[{region}] Exporting mesh to {vtk_file}...")
        export_file(pipeline, str(vtk_file), format="vtk/trimesh", key="surface")
        vtk_files.append(vtk_file)

    print(f"[DONE] Generated {len(vtk_files)} meshes.")

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