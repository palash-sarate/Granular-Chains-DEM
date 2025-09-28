"""Generate the tabulated bond potential used by LAMMPS.

The goal is to keep grains at a nominal centre spacing of 2.5 mm with a very
stiff restoring force once the separation exits a narrow slack window. LAMMPS
only uses the energy column for `bond_style table spline`, so we construct the
energy such that its slope yields the desired force profile. The force column
is emitted for easier inspection but is not consumed by LAMMPS in this style.

Feel free to tweak the parameters below to explore different well widths or
stiffness values. Re-run this script whenever the table needs to be updated.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class BondTableConfig:
    section: str = "WIN"
    r0: float = 2.5e-3             # nominal centre spacing (m)
    slack_half_width: float = 1.5e-4  # half-width of the near-zero-force window (m)
    stiffness: float = 1.0e8       # N/m applied outside the slack window
    r_min: float = 2.0e-3          # lower bound of the tabulated range (m)
    r_max: float = 5.2e-3          # upper bound of the tabulated range (m)
    n_points: int = 401            # samples across the interval
    output_path: Path = Path("bond_window.table")


def compute_energy_force(cfg: BondTableConfig, r: float) -> tuple[float, float]:
    """Return energy (J) and force (N) for distance *r*.

    The potential is flat within the slack window [r0 - w, r0 + w]. Outside that
    region the energy follows a quadratic well with curvature chosen so that the
    resulting force reaches roughly 2.5e5 N after a 0.5 mm excursion.
    """

    lower = cfg.r0 - cfg.slack_half_width
    upper = cfg.r0 + cfg.slack_half_width

    if lower <= r <= upper:
        return 0.0, 0.0

    if r < lower:
        delta = lower - r
        energy = 0.5 * cfg.stiffness * delta * delta
        force = cfg.stiffness * delta  # pulls beads apart (positive force)
        return energy, force

    delta = r - upper
    energy = 0.5 * cfg.stiffness * delta * delta
    force = -cfg.stiffness * delta  # pulls beads back together (negative force)
    return energy, force


def generate_table(cfg: BondTableConfig) -> list[tuple[int, float, float, float]]:
    dr = (cfg.r_max - cfg.r_min) / (cfg.n_points - 1)
    rows: list[tuple[int, float, float, float]] = []

    for idx in range(cfg.n_points):
        r = cfg.r_min + idx * dr
        energy, force = compute_energy_force(cfg, r)
        rows.append((idx + 1, r, energy, force))

    return rows


def write_table(cfg: BondTableConfig) -> None:
    rows = generate_table(cfg)
    cfg.output_path.write_text(_format_table(cfg, rows), encoding="utf-8")


def _format_table(
    cfg: BondTableConfig, rows: list[tuple[int, float, float, float]]
) -> str:
    header = [
        cfg.section,
        f"N {len(rows)}",
        "",
    ]

    body = [
        f"{idx} {r:0.6f} {energy: .12e} {force: .12e}"
        for idx, r, energy, force in rows
    ]

    return "\n".join(header + body) + "\n"

# plot the table for visual inspection
def plot_table(cfg: BondTableConfig, rows: list[tuple[int, float, float, float]]) -> None:
    import matplotlib.pyplot as plt

    rs = [r * 1e3 for _, r, _, _ in rows]  # convert to mm
    energies = [energy for _, _, energy, _ in rows]
    forces = [force for _, _, _, force in rows]

    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel('Distance (mm)')
    ax1.set_ylabel('Energy (J)', color=color)
    ax1.plot(rs, energies, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel('Force (N)', color=color)  # we already handled the x-label with ax1
    ax2.plot(rs, forces, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title('Bond Potential and Force Profile')
    plt.show()

if __name__ == "__main__":
    write_table(BondTableConfig())
    # suggest time step as del t = 0.1 * sqrt(m/k)
    del_t = 0.1 * (1.0e-6 / BondTableConfig().stiffness) ** 0.5
    print (f"Suggested time step: {del_t:.3e} s")
    plot_table(BondTableConfig(), generate_table(BondTableConfig()))