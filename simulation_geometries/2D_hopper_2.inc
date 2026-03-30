# --- Geometric Variables ---
variable x_half      equal 0.00125      # 1.25 mm (Half of 2.5 mm gap)
variable y_half      equal 0.155        # 15.5 cm (Half of 31 cm width)
variable orifice_half equal 0.05        # 5 cm (Half of 10 cm)
variable z_bin_top   equal 0.4          # Height where entry funnel starts
variable funnel_ang  equal 45.0         # Entry funnel angle to horizontal
variable pi equal 3.14159265
variable sin_f       equal sin(${funnel_ang}*${pi}/180.0)
variable cos_f       equal cos(${funnel_ang}*${pi}/180.0)
variable hopper_ang  equal 60.0         # Entry funnel angle to horizontal
variable sin_c       equal sin(${hopper_ang}*${pi}/180.0)
variable cos_c       equal cos(${hopper_ang}*${pi}/180.0)
variable eps equal 0.01
variable funnel_thickness equal ${x_half}+${eps}
variable funnel_ht equal ${z_bin_top}+0.1

# 1. Entry Funnel (Top) - Funneling in the X-direction
# These planes converge from a wide X-opening down to the 2.5mm gap at z_bin_top

region entry_f_neg prism ${x_half} ${funnel_thickness} -${y_half} ${y_half} ${z_bin_top} ${funnel_ht} 0 0.1 0 side out
region entry_f_pos prism -${funnel_thickness} -${x_half} -${y_half} ${y_half} ${z_bin_top} ${funnel_ht} 0 -0.1 0 side out

fix wall_f_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region entry_f_neg
fix wall_f_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region entry_f_pos
# fix wall_f_neg all wall/gran/region granular hertz/material 1e5 1e3 0.3 tangential mindlin NULL 1.0 0.5 region entry_f_neg
# fix wall_f_pos all wall/gran/region granular hertz/material 1e5 1e3 0.3 tangential mindlin NULL 1.0 0.5 region entry_f_pos

# 2. Hopper
# Inclined in the Y-direction as per your previous setup
region exit_v_neg plane 0 -${orifice_half} 0 0  ${sin_c} ${cos_c}
region exit_v_pos plane 0  ${orifice_half} 0 0 -${sin_c} ${cos_c}
region planes_h intersect 2 exit_v_neg exit_v_pos
# Combine the exit planes with front/back limits (X-direction)
region hopper_bounds block -${x_half} ${x_half} -${y_half} ${y_half} 0.0 ${z_bin_top}
region hopper intersect 2 planes_h hopper_bounds

# 6. Apply Granular Wall Fix
fix wall_all all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region hopper

# 5 Orifice cover
# This stops the particles from leaving hopper when filling
region orifice_cover plane 0 0 0 0 0 1
fix lid all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region orifice_cover