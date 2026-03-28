# --- Geometric Variables ---
variable x_half      equal 0.00125      # 1.25 mm (Half of 2.5 mm gap)
variable y_half      equal 0.155        # 15.5 cm (Half of 31 cm width)
variable orifice_half equal 0.05        # 5 cm (Half of 10 cm)
variable z_bin_top   equal 0.4          # Height where entry funnel starts
variable funnel_ang  equal 45.0         # Entry funnel angle to horizontal
variable pi equal 3.14159265
variable sin_f       equal sin(${funnel_ang}*${pi}/180.0)
variable cos_f       equal cos(${funnel_ang}*${pi}/180.0)
variable converging_sec_ang  equal 60.0         # Entry funnel angle to horizontal
variable sin_c       equal sin(${converging_sec_ang}*${pi}/180.0)
variable cos_c       equal cos(${converging_sec_ang}*${pi}/180.0)

# 1. Entry Funnel (Top) - Funneling in the X-direction
# These planes converge from a wide X-opening down to the 2.5mm gap at z_bin_top
region entry_f_neg plane -${x_half} 0 ${z_bin_top} -${sin_f} 0 ${cos_f}
region entry_f_pos plane  ${x_half} 0 ${z_bin_top}  ${sin_f} 0 ${cos_f}

# 2. Main Hopper Bin (Middle) - The 2.5mm thin section
region bin_middle block -${x_half} ${x_half} -${y_half} ${y_half} 0.0 ${z_bin_top}

# 3. Converging Section (Bottom) - The 60-degree exit
# Inclined in the Y-direction as per your previous setup
region exit_v_neg plane 0 -${orifice_half} 0 0 -${sin_c} ${cos_c}
region exit_v_pos plane 0  ${orifice_half} 0 0  ${sin_c} ${cos_c}

# 4. Side Walls (Y-limits) for the Top and Middle sections
# This ensures particles don't escape the 31 cm width
region entry_funnel_bounds block -1 1 -${y_half} ${y_half} 0.0 1.0

# 5. Logical Assembly
# Combine the entry planes with side limits
region funnel_top intersect 3 entry_f_neg entry_f_pos entry_funnel_bounds

# Combine the exit planes with front/back limits (X-direction)
region xyz_limits block -${x_half} ${x_half} -${y_half} ${y_half} 0.0 ${z_bin_top}
region funnel_bottom intersect 3 exit_v_neg exit_v_pos xyz_limits

# 5 Orifice cover
# This stops the particles from leaving hopper when filling
# Orifice cover plane (blocks particles from leaving at z = 0)
region orifice_plane plane 0.0 0.0 0.0   0.0 0.0 1.0

# Limit it to hopper cross-section
region orifice_cover intersect 2 orifice_plane xyz_limits

fix lid all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region orifice_cover

# Final Union of all parts
region complete_hopper union 3 funnel_top bin_middle funnel_bottom

# 6. Apply Granular Wall Fix
fix wall_all all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region complete_hopper