"""Utilities for generating LAMMPS `chain.data` files.

This module builds linear bead chains compatible with the existing
`in.chain_flop.lmp` workflow.  You can import and re‑use the helper
functions or execute the module directly from the command line.

Example
-------

Generate a vertical chain with 12 beads spaced 2.5 mm apart::

	python chain_generator.py 12 --orientation vert --spacing 2.5e-3

The output file will be written as ``chain_vert_12.data`` in the current
working directory unless ``--output-dir`` is provided.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


ORIENTATION_AXES = {
	"vert": 2,  # z-axis
	"horz": 0,  # x-axis
}


@dataclass
class ChainConfig:
	"""Configuration parameters for a bead chain."""

	beads: int
	orientation: str
	spacing: float  # centre-to-centre spacing in metres
	diameter: float = 2.0e-3  # metres
	density: float = 8.5e3  # kg/m^3
	output_dir: Path = Path(".")

	def validate(self) -> None:
		if self.beads < 1:
			raise ValueError("The chain must contain at least one bead.")
		if self.spacing <= 0:
			raise ValueError("Spacing must be positive and > 2bd.")
		if self.diameter <= 0:
			raise ValueError("Diameter must be positive.")
		if self.density <= 0:
			raise ValueError("Density must be positive.")
		if self.orientation not in ORIENTATION_AXES:
			raise ValueError(
				f"Unsupported orientation '{self.orientation}'. Use one of {tuple(ORIENTATION_AXES)}."
			)

	@property
	def output_path(self) -> Path:
		name = f"chain_{self.orientation}_{self.beads}.data"
		return self.output_dir / name

	@property
	def radius(self) -> float:
		return 0.5 * self.diameter

	@property
	def bead_mass(self) -> float:
		volume = (4.0 / 3.0) * math.pi * self.radius**3
		return volume * self.density


def generate_positions(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	"""Return the bead centre coordinates for the requested chain."""

	axis = ORIENTATION_AXES[cfg.orientation]
	span = cfg.spacing * (cfg.beads - 1)
	start = -0.5 * span

	positions: List[List[float]] = []
	for bead in range(cfg.beads):
		coord = [0.0, 0.0, 0.0]
		coord[axis] = start + bead * cfg.spacing
		positions.append(coord)

	return [tuple(values) for values in positions]


def estimate_box(cfg: ChainConfig, positions: Sequence[Tuple[float, float, float]]) -> Tuple[float, float, float, float, float, float]:
	"""Compute a bounding box that comfortably encloses the chain."""

	margin = max(cfg.spacing, cfg.diameter)
	radius = cfg.radius

	xs, ys, zs = zip(*positions)

	xlo = min(xs) - radius - margin
	xhi = max(xs) + radius + margin
	ylo = min(ys) - radius - margin
	yhi = max(ys) + radius + margin
	zlo = min(zs) - radius - margin
	zhi = max(zs) + radius + margin

	# Guarantee non-zero extents even for single-bead chains aligned to zero.
	if math.isclose(xlo, xhi):
		xlo -= margin
		xhi += margin
	if math.isclose(ylo, yhi):
		ylo -= margin
		yhi += margin
	if math.isclose(zlo, zhi):
		zlo -= margin
		zhi += margin

	return xlo, xhi, ylo, yhi, zlo, zhi


def format_float(value: float) -> str:
	return f"{value:.16e}"


def write_chain_data(cfg: ChainConfig) -> Path:
	"""Build and write the LAMMPS data file, returning the output path."""

	cfg.validate()
	positions = generate_positions(cfg)
	box = estimate_box(cfg, positions)
	mass = cfg.bead_mass

	lines: List[str] = []
	lines.append(f"LAMMPS data file for {cfg.orientation} chain of {cfg.beads} beads")
	lines.append("")
	lines.append(f"{cfg.beads} atoms")
	lines.append("1 atom types")

	bond_count = max(cfg.beads - 1, 0)
	angle_count = max(cfg.beads - 2, 0)
	lines.append(f"{bond_count} bonds")
	lines.append("1 bond types")
	lines.append(f"{angle_count} angles")
	lines.append("1 angle types")
	lines.append("")

	xlo, xhi, ylo, yhi, zlo, zhi = box
	lines.append(f"{format_float(xlo)} {format_float(xhi)} xlo xhi")
	lines.append(f"{format_float(ylo)} {format_float(yhi)} ylo yhi")
	lines.append(f"{format_float(zlo)} {format_float(zhi)} zlo zhi")
	lines.append("")

	lines.append("Masses")
	lines.append("")
	lines.append(f"1 {format_float(mass)}")
	lines.append("")

	lines.append("Atoms # hybrid sphere molecular")
	lines.append("")

	for atom_id, (x, y, z) in enumerate(positions, start=1):
		lines.append(
			" ".join(
				[
					str(atom_id),
					"1",  # molecule-ID (single chain)
					"1",  # atom type
					format_float(x),
					format_float(y),
					format_float(z),
					format_float(cfg.diameter),
					format_float(cfg.density),
					"1.0",  # quaternion w
					"0.0",
					"0.0",
					"0.0",
				]
			)
		)

	lines.append("")
	lines.append("Velocities")
	lines.append("")
	for atom_id in range(1, cfg.beads + 1):
		lines.append(f"{atom_id} 0 0 0 0 0 0")

	lines.append("")
	lines.append("Bonds")
	lines.append("")
	for bond_id in range(1, bond_count + 1):
		lines.append(f"{bond_id} 1 {bond_id} {bond_id + 1}")

	lines.append("")
	lines.append("Angles")
	lines.append("")
	for angle_id in range(1, angle_count + 1):
		lines.append(f"{angle_id} 1 {angle_id} {angle_id + 1} {angle_id + 2}")

	output_path = cfg.output_path
	cfg.output_dir.mkdir(parents=True, exist_ok=True)
	output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

	return output_path


def parse_args(argv: Iterable[str] | None = None) -> ChainConfig:
	parser = argparse.ArgumentParser(description="Generate linear bead-chain LAMMPS data files.")
	parser.add_argument("beads", type=int, help="Number of beads in the chain.")
	parser.add_argument(
		"--orientation",
		choices=ORIENTATION_AXES.keys(),
		default="vert",
		help="Chain orientation: 'vert' aligns beads along z, 'horz' along x.",
	)
	parser.add_argument(
		"--spacing",
		type=float,
		required=True,
		help="Centre-to-centre spacing between beads (metres).",
	)
	parser.add_argument("--diameter", type=float, default=2.0e-3, help="Bead diameter (metres).")
	parser.add_argument("--density", type=float, default=8.5e3, help="Bead material density (kg/m^3).")
	parser.add_argument(
		"--output-dir",
		type=Path,
		default=Path("."),
		help="Directory where the data file will be written.",
	)

	args = parser.parse_args(argv)

	return ChainConfig(
		beads=args.beads,
		orientation=args.orientation,
		spacing=args.spacing,
		diameter=args.diameter,
		density=args.density,
		output_dir=args.output_dir,
	)


def main(argv: Iterable[str] | None = None) -> None:
	cfg = parse_args(argv)
	output_path = write_chain_data(cfg)
	print(f"Wrote {output_path}")
	print(f"Bead mass: {cfg.bead_mass:.6e} kg")
	print(
		"Spacing: {0:.6e} m, diameter: {1:.6e} m, density: {2:.4e} kg/m^3".format(
			cfg.spacing, cfg.diameter, cfg.density
		)
	)
	print("For stiff bonds, consider Δt ≲ 0.1 * sqrt(m/k) using the mass above.")


if __name__ == "__main__":
	main()

