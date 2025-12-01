"""Utilities for generating LAMMPS `chain.data` files.

This module builds bead chains in various configurations (linear, directional,
random, loop) compatible with the existing `in.chain_flop.lmp` workflow.

You can import and re‑use the helper functions or execute the module directly
from the command line.

Run with --help for detailed usage and examples.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


ORIENTATION_AXES = {
	"vert": 2,  # z-axis
	"horz": 0,  # x-axis
}

MODES = ["linear", "direction", "random", "loop"]


@dataclass
class ChainConfig:
	"""Configuration parameters for a bead chain."""

	beads: int
	orientation: str = "vert"
	spacing: float = 0.0 # centre-to-centre spacing in metres
	diameter: float = 2.0e-3  # metres
	density: float = 7084.7  # kg/m^3
	output_dir: Path = Path("chain_data/default")
	
	# New configuration fields
	mode: str = "linear"
	direction: Tuple[float, float, float] = (0.0, 0.0, 1.0)
	max_spacing: float = 0.0
	loop_radius: float = 0.0
	box: Tuple[float, float, float, float, float, float] | None = None

	def validate(self) -> None:
		if self.beads < 1:
			raise ValueError("The chain must contain at least one bead.")
		if self.spacing <= 0:
			raise ValueError("Spacing must be positive.")
		if self.diameter <= 0:
			raise ValueError("Diameter must be positive.")
		if self.density <= 0:
			raise ValueError("Density must be positive.")
		
		if self.mode not in MODES:
			raise ValueError(f"Unsupported mode '{self.mode}'. Use one of {MODES}.")
			
		if self.mode == "linear" and self.orientation not in ORIENTATION_AXES:
			raise ValueError(
				f"Unsupported orientation '{self.orientation}'. Use one of {tuple(ORIENTATION_AXES)}."
			)
			
		if self.mode == "random" and self.max_spacing > 0 and self.max_spacing < self.diameter:
			raise ValueError("Max spacing must be >= diameter.")

	@property
	def output_path(self) -> Path:
		if self.mode == "linear":
			name = f"N{self.beads}_chain_{self.orientation}.data"
		else:
			name = f"N{self.beads}_chain_{self.mode}.data"
		return Path(f"{self.output_dir}/{name}")

	@property
	def radius(self) -> float:
		return 0.5 * self.diameter

	@property
	def bead_mass(self) -> float:
		volume = (4.0 / 3.0) * math.pi * self.radius**3
		return volume * self.density


def normalize(v: Tuple[float, float, float]) -> Tuple[float, float, float]:
	norm = math.sqrt(sum(x*x for x in v))
	if norm == 0: return (0.0, 0.0, 1.0)
	return (v[0]/norm, v[1]/norm, v[2]/norm)


def generate_positions(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	"""Return the bead centre coordinates for the requested chain."""
	if cfg.mode == "linear":
		return _generate_linear(cfg)
	elif cfg.mode == "direction":
		return _generate_direction(cfg)
	elif cfg.mode == "random":
		return _generate_random(cfg)
	elif cfg.mode == "loop":
		return _generate_loop(cfg)
	else:
		raise ValueError(f"Unknown mode: {cfg.mode}")


def _generate_linear(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	axis = ORIENTATION_AXES[cfg.orientation]
	span = cfg.spacing * (cfg.beads - 1)
	start = -0.5 * span

	positions: List[List[float]] = []
	for bead in range(cfg.beads):
		coord = [0.0, 0.0, 0.0]
		coord[axis] = start + bead * cfg.spacing
		positions.append(coord)

	return [tuple(values) for values in positions]


def _generate_direction(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	ux, uy, uz = normalize(cfg.direction)
	span = cfg.spacing * (cfg.beads - 1)
	start_dist = -0.5 * span
	
	positions = []
	for i in range(cfg.beads):
		d = start_dist + i * cfg.spacing
		positions.append((d * ux, d * uy, d * uz))
	return positions


def _generate_random(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	positions = [(0.0, 0.0, 0.0)]
	for _ in range(cfg.beads - 1):
		prev = positions[-1]
		
		# Uniform random direction on sphere
		# z is uniform in [-1, 1]
		# phi is uniform in [0, 2pi]
		z = random.uniform(-1.0, 1.0)
		phi = random.uniform(0.0, 2.0 * math.pi)
		
		r_xy = math.sqrt(1.0 - z*z)
		dx = r_xy * math.cos(phi)
		dy = r_xy * math.sin(phi)
		dz = z
		
		# Spacing
		if cfg.max_spacing > 0:
			dist = random.uniform(cfg.diameter, cfg.max_spacing)
		else:
			dist = cfg.spacing
			
		positions.append((prev[0] + dist * dx, prev[1] + dist * dy, prev[2] + dist * dz))
	return positions


def _generate_loop(cfg: ChainConfig) -> List[Tuple[float, float, float]]:
	# Parameters
	# Default radius is 2 * diameter if not specified
	R = cfg.loop_radius if cfg.loop_radius > 0 else 2.0 * cfg.diameter
	pitch = 1.5 * cfg.diameter # Small pitch to avoid intersection
	
	# Calculate beads for one full loop (2pi)
	loop_length = math.sqrt((2 * math.pi * R)**2 + pitch**2)
	beads_per_loop = int(math.ceil(loop_length / cfg.spacing))
	
	n_loop = beads_per_loop
	remaining = cfg.beads - n_loop
	
	if remaining < 0:
		# If requested beads are fewer than a full loop, just make the loop part
		n_loop = cfg.beads
		n_start = 0
		n_end = 0
	else:
		# Distribute remaining beads to start and end straight sections
		n_start = remaining // 2
		n_end = remaining - n_start
		
	positions = []
	
	# 1. Straight part (along Y axis, ending at origin)
	# We align so that the straight part flows into the helix at (R, 0, 0)
	# Helix starts at t=0 => (R, 0, 0). Tangent is (0, R, h/2pi) ~ +Y
	# So straight part comes from -Y direction.
	
	y_start = - (n_start * cfg.spacing)
	for i in range(n_start):
		# x=R, z=0, y increasing
		positions.append((R, y_start + i * cfg.spacing, 0.0))
		
	# 2. Loop part (Helix)
	dt = cfg.spacing / math.sqrt(R**2 + (pitch / (2 * math.pi))**2)
	current_t = 0.0
	
	# If we had straight part, the last point was (R, -spacing, 0).
	# The first loop point at t=0 is (R, 0, 0). Distance is spacing. Perfect.
	# If n_start=0, we start at (R, 0, 0).
	
	for i in range(n_loop):
		x = R * math.cos(current_t)
		y = R * math.sin(current_t)
		z = (pitch / (2 * math.pi)) * current_t
		positions.append((x, y, z))
		current_t += dt
		
	# 3. Straight part end
	# Continue along the tangent of the helix at the last point
	last_t = current_t - dt
	last_pos = positions[-1]
	
	# Tangent vector at last_t
	tx = -R * math.sin(last_t)
	ty = R * math.cos(last_t)
	tz = pitch / (2 * math.pi)
	norm = math.sqrt(tx*tx + ty*ty + tz*tz)
	ux, uy, uz = tx/norm, ty/norm, tz/norm
	
	for i in range(n_end):
		px = last_pos[0] + (i + 1) * cfg.spacing * ux
		py = last_pos[1] + (i + 1) * cfg.spacing * uy
		pz = last_pos[2] + (i + 1) * cfg.spacing * uz
		positions.append((px, py, pz))
		
	return positions


def estimate_box(cfg: ChainConfig, positions: Sequence[Tuple[float, float, float]]) -> Tuple[float, float, float, float, float, float]:
	"""Compute a bounding box that comfortably encloses the chain."""

	# Calculate geometric center
	xs, ys, zs = zip(*positions)
	center_x = (min(xs) + max(xs)) / 2.0
	center_y = (min(ys) + max(ys)) / 2.0
	center_z = (min(zs) + max(zs)) / 2.0

	# Create a generous cubic box: size = 2 * number of beads * 0.0025
	box_size = 2.0 * cfg.beads * 0.003
	half_size = box_size / 2.0

	xlo = center_x - half_size
	xhi = center_x + half_size
	ylo = center_y - half_size
	yhi = center_y + half_size
	zlo = center_z - half_size
	zhi = center_z + half_size

	return xlo, xhi, ylo, yhi, zlo, zhi


def format_float(value: float) -> str:
	return f"{value:.16e}"


def write_chain_data(cfg: ChainConfig) -> Path:
	"""Build and write the LAMMPS data file, returning the output path."""

	cfg.validate()
	positions = generate_positions(cfg)
	
	if cfg.box:
		box = cfg.box
	else:
		box = estimate_box(cfg, positions)
		
	mass = cfg.bead_mass

	lines: List[str] = []
	lines.append(f"LAMMPS data file for {cfg.mode} chain of {cfg.beads} beads")
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
					"1",  # atom type
					format_float(x),
					format_float(y),
					format_float(z),
					format_float(cfg.diameter),
					format_float(cfg.density),
					"1",  # molecule-ID (single chain)
					"0",
					"0",
					"0",
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
	description = """Generate bead-chain LAMMPS data files with various geometries.

Modes:
  linear    - Straight chain along a primary axis (vert/z or horz/x).
  direction - Straight chain along an arbitrary 3D vector.
  random    - 3D random walk (squiggle).
  loop      - Straight section flowing into a helical loop and back out.
"""
	epilog = """
Examples:
  1. Linear vertical chain (default):
     python chain_generator.py 50 --spacing 0.0025

  2. Linear horizontal chain:
     python chain_generator.py 50 --mode linear --orientation horz --spacing 0.0025

  3. Chain along a specific direction vector (1, 1, 1):
     python chain_generator.py 50 --mode direction --direction 1 1 1 --spacing 0.0025

  4. Random walk chain with variable spacing (2mm to 4mm):
     python chain_generator.py 100 --mode random --spacing 0.002 --max-spacing 0.004

  5. Loop configuration:
     python chain_generator.py 100 --mode loop --spacing 0.0025
"""
	parser = argparse.ArgumentParser(
		description=description,
		epilog=epilog,
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument("beads", type=int, help="Number of beads in the chain.")
	
	# Mode selection
	parser.add_argument(
		"--mode",
		choices=MODES,
		default="linear",
		help="Generation mode: 'linear', 'direction', 'random', or 'loop'.",
	)

	# Linear mode options
	parser.add_argument(
		"--orientation",
		choices=ORIENTATION_AXES.keys(),
		default="vert",
		help="[Linear] Chain orientation: 'vert' aligns beads along z, 'horz' along x.",
	)
	
	# Direction mode options
	parser.add_argument(
		"--direction",
		type=float,
		nargs=3,
		default=[0.0, 0.0, 1.0],
		help="[Direction] Vector defining the chain direction (x y z).",
	)
	
	# Random mode options
	parser.add_argument(
		"--max-spacing",
		type=float,
		default=0.0,
		help="[Random] Maximum spacing for random walk. If 0, uses fixed spacing.",
	)
	
	# Loop mode options
	parser.add_argument(
		"--loop-radius",
		type=float,
		default=0.0,
		help="[Loop] Radius of the loop section. Defaults to 2 * diameter.",
	)

	# Common options
	parser.add_argument(
		"--spacing",
		type=float,
		required=True,
		help="Centre-to-centre spacing between beads (metres).",
	)
	parser.add_argument("--diameter", type=float, default=ChainConfig.diameter, help="Bead diameter (metres).")
	parser.add_argument("--density", type=float, default=ChainConfig.density, help="Bead material density (kg/m^3).")
	parser.add_argument(
		"--output-dir",
		type=Path,
		default=Path("chain_datas"),
		help="Directory where the data file will be written.",
	)
	
	parser.add_argument(
		"--box",
		type=float,
		nargs=6,
		metavar=("XLO", "XHI", "YLO", "YHI", "ZLO", "ZHI"),
		help="Override simulation box bounds (xlo xhi ylo yhi zlo zhi).",
	)

	args = parser.parse_args(argv)

	return ChainConfig(
		beads=args.beads,
		mode=args.mode,
		orientation=args.orientation,
		direction=tuple(args.direction),
		spacing=args.spacing,
		max_spacing=args.max_spacing,
		loop_radius=args.loop_radius,
		diameter=args.diameter,
		density=args.density,
		output_dir=args.output_dir,
		box=tuple(args.box) if args.box else None,
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

