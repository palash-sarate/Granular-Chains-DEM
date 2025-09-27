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

# # Print first 40 lines for inspection
# print("Wrote:", outfile)
# print("--- first 40 lines ---")
# with outfile.open() as fh:
#     for i, line in enumerate(fh):
#         if i<40:
#             print(line.rstrip())
#         else:
#             break
