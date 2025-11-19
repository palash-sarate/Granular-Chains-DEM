"""Generate tabulated angle potentials for LAMMPS.

This utility constructs an "infinite" angular well by providing a
zero-torque plateau inside a user-specified window and quadratic walls
outside the window.  The resulting file can be consumed with
``angle_style table spline`` using the generated section name
``ANGLE_WELL``.

Run ``python Angle_table_generator.py --help`` for CLI usage details.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np


@dataclass
class AngleWellConfig:
    """Configuration for the angular well generator."""

    theta0: float = 38.0  # centre of the window (deg)
    window: float = 38.0  # total width with zero torque (deg)
    stiffness: float = 1.0  # quadratic stiffness outside the window
    theta_min: float = 0.0  # lower limit of the tabulated range (deg)
    theta_max: float = 180.0  # upper limit of the tabulated range (deg)
    n_points: int = 2001  # resolution of the table
    output: Path = Path("bond_angle.table")
    plot: bool = False

    def validate(self) -> None:
        if self.window < 0:
            raise ValueError("window must be non-negative")
        if self.n_points < 2:
            raise ValueError("n_points must be >= 2")
        if self.theta_min >= self.theta_max:
            raise ValueError("theta_min must be smaller than theta_max")
        if self.stiffness < 0:
            raise ValueError("stiffness must be non-negative")


def generate_angle_profile(cfg: AngleWellConfig) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return angle grid, energy, and torque arrays."""

    theta = np.linspace(cfg.theta_min, cfg.theta_max, cfg.n_points)
    half_window = 0.5 * cfg.window
    lo = cfg.theta0 - half_window
    hi = cfg.theta0 + half_window

    # Distances outside the allowed window.  Any negative value means the
    # point lies inside the window so we clip it to zero.
    delta_low = np.clip(lo - theta, 0.0, None)
    delta_high = np.clip(theta - hi, 0.0, None)

    # Quadratic energy contributions on each side.
    energy = 0.5 * cfg.stiffness * (delta_low ** 2 + delta_high ** 2)

    # Torque is -dU/dtheta.  Differentiating the quadratic gives the
    # expressions below, which are continuous at the window boundaries.
    torque = cfg.stiffness * delta_low - cfg.stiffness * delta_high

    return theta, energy, torque


def write_table(cfg: AngleWellConfig, theta: np.ndarray, energy: np.ndarray, torque: np.ndarray) -> Path:
    section_name = "ANGLE_WELL"
    header = (
        "# LAMMPS angle table generated with Angle_table_generator.py\n"
        f"# Zero-torque window centred at {cfg.theta0:.6f} deg with width {cfg.window:.6f} deg\n"
        "# Columns: index  theta(deg)  energy  torque\n"
    )

    lines = [header, section_name + "\n", f"N {len(theta)}\n", "\n"]
    for idx, (th, en, tq) in enumerate(zip(theta, energy, torque), start=1):
        lines.append(f"{idx} {th:.6f} {en:.12e} {tq:.12e}\n")

    cfg.output.write_text("".join(lines), encoding="utf-8")
    return cfg.output


def plot_profile(theta: np.ndarray, energy: np.ndarray, torque: np.ndarray) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - plotting is optional
        raise RuntimeError("matplotlib is required for plotting") from exc

    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Angle (deg)")
    ax1.set_ylabel("Energy", color="tab:blue")
    ax1.plot(theta, energy, color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Torque", color="tab:red")
    ax2.plot(theta, torque, color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")

    plt.title("Angular well energy and torque")
    plt.grid(True)
    fig.tight_layout()
    plt.show()


def parse_args(argv: Iterable[str] | None = None) -> AngleWellConfig:
    parser = argparse.ArgumentParser(description="Generate an angular well table for LAMMPS.")
    parser.add_argument("--theta0", type=float, default=38.0, help="Centre angle of the zero-torque window (deg).")
    parser.add_argument(
        "--window",
        type=float,
        default=38.0,
        help="Full angular width of the zero-torque window (deg).",
    )
    parser.add_argument(
        "--stiffness",
        type=float,
        default=1.0,
        help="Quadratic stiffness applied outside the window (energy units per deg^2).",
    )
    parser.add_argument("--theta-min", type=float, default=0.0, help="Minimum tabulated angle (deg).")
    parser.add_argument("--theta-max", type=float, default=180.0, help="Maximum tabulated angle (deg).")
    parser.add_argument("--points", type=int, default=2001, help="Number of sample points between the limits.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("bond_angle.table"),
        help="Path to the angle table that will be written.",
    )
    parser.add_argument("--plot", action="store_true", help="Display a plot of the generated table.")

    args = parser.parse_args(argv)

    cfg = AngleWellConfig(
        theta0=args.theta0,
        window=args.window,
        stiffness=args.stiffness,
        theta_min=args.theta_min,
        theta_max=args.theta_max,
        n_points=args.points,
        output=args.output,
        plot=args.plot,
    )
    cfg.validate()
    return cfg


def main(argv: Iterable[str] | None = None) -> None:
    cfg = parse_args(argv)
    theta, energy, torque = generate_angle_profile(cfg)
    output_path = write_table(cfg, theta, energy, torque)
    print(f"Wrote {output_path}")

    if cfg.plot:
        plot_profile(theta, energy, torque)


if __name__ == "__main__":
    main()
# Generating a LAMMPS angle table file that implements an "infinite" angular well:
# - Assumption: user wants a free (zero-torque) plateau of width 38 degrees centered at theta0=38 deg
#   (i.e., theta in [19, 57] degrees is free). Adjust `window_width` if you meant differently.
# - Outside that window we make a smooth quadratic 'wall' with stiffness K_wall so motion is practically forbidden.
# - The table format produced is: index  angle(deg)  energy  torque
#   (LAMMPS angle_style table expects angles in degrees when you use degrees, and torque = -dU/dtheta)
# - File written to /mnt/data/angle_well.table
#
# You can change: theta0, window_width, delta (transition width in deg), K_wall (stiffness of walls),
#      N_points (resolution across 0..180 deg).
#
# After this I will print the first 40 lines so you can inspect, and provide a download link.

import numpy as np
from pathlib import Path

# Parameters (edit here if you want different)
theta0 = 38.0           # center (deg)
window_width = 38.0     # free region: [theta0 - width/2, theta0 + width/2]
delta = 3.0             # transition half-width (deg) for smooth connection to the wall
K_wall = 1.0            # stiffness of the quadratic wall (energy units per deg^2)
N_points = 2001         # number of table points (cover 0..180 deg)
outfile = Path("./bond_angle.table")

# Derived
theta_min = 0.0
theta_max = 180.0
theta = np.linspace(theta_min, theta_max, N_points)

# allowed window
w_half = window_width / 2.0
allowed_lo = theta0 - w_half
allowed_hi = theta0 + w_half

# energy function: smooth quadratic walls that go to 0 inside [allowed_lo+delta, allowed_hi-delta]
U = np.zeros_like(theta)
torque = np.zeros_like(theta)  # torque = -dU/dtheta (deg units)

# Build left wall for theta < allowed_lo
left_transition_start = allowed_lo - delta
left_transition_end = allowed_lo + delta  # within this we'll smoothly go to zero energy at allowed_lo+delta
# similarly right wall
right_transition_start = allowed_hi - delta
right_transition_end = allowed_hi + delta

def plot_table(cfg, rows):
    import matplotlib.pyplot as plt

    angles = [r[1] for r in rows]
    energies = [r[2] for r in rows]
    torques = [r[3] for r in rows]

    fig, ax1 = plt.subplots()

    color = "tab:blue"
    ax1.set_xlabel("Angle (deg)")
    ax1.set_ylabel("Energy", color=color)
    ax1.plot(angles, energies, color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = "tab:red"
    ax2.set_ylabel("Torque", color=color)  # we already handled the x-label with ax1
    ax2.plot(angles, torques, color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title("Angle well energy and torque")
    plt.grid()
    plt.show()

# Helper: smoothstep for blending
def smoothstep(x):
    # cubic smoothstep from 0->1 for x in [0,1]
    return 3*x**2 - 2*x**3

for i, th in enumerate(theta):
    if th <= left_transition_start:
        # fully in left hard wall: quadratic anchored at left_transition_start (so U large when far)
        dx = left_transition_start - th
        U[i] = 0.5 * K_wall * (dx)**2
        torque[i] = -(-K_wall * dx)  # -dU/dtheta ; dU/dtheta = -K_wall*dx -> torque = -dU/dtheta = K_wall*dx
    elif th < left_transition_end:
        # in smoothing region: blend quadratic (at left_transition_start) to zero at left_transition_end (=allowed_lo+delta)
        x = (th - left_transition_start) / (left_transition_end - left_transition_start)  # 0->1
        blend = smoothstep(x)
        # quadratic value at left_transition_start:
        dx0 = 0.0  # at left_transition_start itself dx=0
        # but energy should be 0 at left_transition_end
        # simpler approach: define a quadratic wall anchored at allowed_lo (so U=0 at allowed_lo)
        dx = allowed_lo - th
        Uq = 0.5 * K_wall * (dx)**2
        U[i] = blend * 0.0 + (1-blend) * Uq
        # torque -> derivative of Uq = -K_wall*dx ; include blend derivative approx
        torque[i] = -(-K_wall * dx) * (1-blend)  # approximate (ok since delta small)
    elif th <= right_transition_start:
        # fully inside plateau -> zero energy and torque
        U[i] = 0.0
        torque[i] = 0.0
    elif th < right_transition_end:
        # smoothing region at right side
        x = (th - right_transition_start) / (right_transition_end - right_transition_start)
        blend = smoothstep(x)
        dx = th - allowed_hi
        Uq = 0.5 * K_wall * (dx)**2
        U[i] = (1-blend) * 0.0 + blend * Uq
        torque[i] = -(K_wall * dx) * blend  # -dU/dtheta ; dU/dtheta = K_wall*dx
    else:
        # fully in right hard wall
        dx = th - right_transition_end
        # anchor quadratic at right_transition_end
        Uq = 0.5 * K_wall * (dx + delta)**2  # shift so energy increases away from allowed_hi
        # simpler: use quadratic about allowed_hi
        dx2 = th - allowed_hi
        U[i] = 0.5 * K_wall * (dx2)**2
        torque[i] = -(K_wall * dx2)

# Ensure torque = -dU/dtheta numerically (impose consistency) by numerical derivative of U
dU = np.gradient(U, theta)  # dU/dtheta (per degree)
torque = -dU

# Write file in a simple LAMMPS-friendly table format:
# Section name header then lines: index  theta(deg)  energy  torque
section_name = "ANGLE_WELL"
with outfile.open("w", newline="\n") as fh:
    fh.write("# LAMMPS angle table generated: allowed window [{:.3f},{:.3f}] deg centered at {:.3f} deg\n".format(allowed_lo, allowed_hi, theta0))
    fh.write("# Columns: index  theta(deg)  energy  torque\n")
    fh.write(section_name + "\n")
    fh.write(f"N {N_points}\n")
    # fh.write(f"theta_min {theta_min:.6f}\n")
    # fh.write(f"theta_max {theta_max:.6f}\n")
    fh.write("\n")
    for idx, (th, u, t) in enumerate(zip(theta, U, torque), start=1):
        fh.write(f"{idx} {th:.6f} {u:.12e} {t:.12e}\n")

# plot the table for visual inspection
plot_table(None, [(i+1, th, u, t) for i, (th, u, t) in enumerate(zip(theta, U, torque))])
