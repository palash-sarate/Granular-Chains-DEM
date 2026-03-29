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
# variable converging_sec_ht equal tan(${hopper_ang}*${pi}/180.0)*(${y_half} - ${orifice_half})
variable x_cover equal 4*${x_half}

# 1. Entry Funnel (Top) - Funneling in the X-direction
# These planes converge from a wide X-opening down to the 2.5mm gap at z_bin_top
region entry_f_neg plane -${x_half} 0 ${z_bin_top} -${sin_f} 0 ${cos_f}
region entry_f_pos plane  ${x_half} 0 ${z_bin_top}  ${sin_f} 0 ${cos_f}
# This ensures particles don't escape the funnel
region funnel_bounds block -0.5 0.5 -${y_half} ${y_half} ${z_bin_top} 0.5 open 5 open 6
# Combine the entry planes with side limits
region funnel intersect 3 entry_f_neg entry_f_pos funnel_bounds

# 2. Hopper
# Inclined in the Y-direction as per your previous setup
region exit_v_neg plane 0 -${orifice_half} 0 0  ${sin_c} ${cos_c}
region exit_v_pos plane 0  ${orifice_half} 0 0 -${sin_c} ${cos_c}

# Combine the exit planes with front/back limits (X-direction)
region hopper_bounds block -${x_half} ${x_half} -${y_half} ${y_half} 0.0 ${z_bin_top} open 5 open 6
region hopper intersect 3 exit_v_neg exit_v_pos hopper_bounds


# 5 Orifice cover
# This stops the particles from leaving hopper when filling
region orifice_cover block -${x_cover} ${x_cover} -${y_half} ${y_half} -0.005 0.0

fix lid all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region orifice_cover

# Final Union of all parts
region complete_hopper union 2 funnel hopper

# 6. Apply Granular Wall Fix
fix wall_all all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region complete_hopper