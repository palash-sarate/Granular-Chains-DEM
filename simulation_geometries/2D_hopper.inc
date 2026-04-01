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
variable converging_sec_ht equal tan(${hopper_ang}*${pi}/180.0)*(${y_half} - ${orifice_half})
variable eps equal 0.01
variable funnel_top equal 0.5 # funnel extent in z direction
variable f_x_half equal 0.1 # funnel extent in x direction
variable h_y_half equal ${y_half}+${eps} # hopper extent in y direction
variable h_x_half equal ${x_half}+${eps} # hopper extent in y direction

# use _vis regions for visualising geometry

# Entry funnel (Top) - Funneling in the X-direction
region entry_f_neg plane -${x_half} 0 ${z_bin_top} -${sin_f} 0 -${cos_f}
region entry_f_pos plane  ${x_half} 0 ${z_bin_top}  ${sin_f} 0 -${cos_f}
region funnel_bounds block -${f_x_half} ${f_x_half} -${y_half} ${y_half} ${z_bin_top} ${funnel_top}
region funnel_neg_vis intersect 2 entry_f_neg funnel_bounds
region funnel_pos_vis intersect 2 entry_f_pos funnel_bounds
region funnel_neg intersect 2 entry_f_neg funnel_bounds side out
region funnel_pos intersect 2 entry_f_pos funnel_bounds side out
fix wall_funnel_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region funnel_neg
fix wall_funnel_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region funnel_pos
# Funnel side walls
region funnel_left_vis block -${f_x_half} ${f_x_half} ${y_half} ${h_y_half} ${z_bin_top} ${funnel_top}
region funnel_right_vis block -${f_x_half} ${f_x_half} -${h_y_half} -${y_half} ${z_bin_top} ${funnel_top}
region funnel_left block -${f_x_half} ${f_x_half} ${y_half} ${h_y_half} ${z_bin_top} ${funnel_top} side out
region funnel_right block -${f_x_half} ${f_x_half} -${h_y_half} -${y_half} ${z_bin_top} ${funnel_top} side out
fix wall_funnel_left all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region funnel_left
fix wall_funnel_right all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region funnel_right
# hopper side walls
region hopper_pos_vis block -${x_half} ${x_half} ${y_half} ${h_y_half} ${converging_sec_ht} ${z_bin_top}
region hopper_neg_vis block -${x_half} ${x_half} -${h_y_half} -${y_half} ${converging_sec_ht} ${z_bin_top}
region hopper_pos block -${x_half} ${x_half} ${y_half} ${h_y_half} ${converging_sec_ht} ${z_bin_top} side out
region hopper_neg block -${x_half} ${x_half} -${h_y_half} -${y_half} ${converging_sec_ht} ${z_bin_top} side out
fix wall_hopper_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region hopper_pos
fix wall_hopper_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region hopper_neg
# hopper front and back covers
region hopper_front_vis block ${x_half} ${h_x_half} -${y_half} ${y_half} 0 ${z_bin_top}
region hopper_back_vis block -${h_x_half} -${x_half} -${y_half} ${y_half} 0 ${z_bin_top}
region hopper_front block ${x_half} ${h_x_half} -${y_half} ${y_half} 0 ${z_bin_top} side out
region hopper_back block -${h_x_half} -${x_half} -${y_half} ${y_half} 0 ${z_bin_top} side out
fix wall_hopper_front all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region hopper_front
fix wall_hopper_back all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region hopper_back
# Converging section
region exit_v_neg plane 0 -${orifice_half} 0 0 -${sin_c} -${cos_c}
region exit_v_pos plane 0  ${orifice_half} 0 0  ${sin_c} -${cos_c}
region exit_bounds block -${x_half} ${x_half} -${y_half} ${y_half} 0.0 ${converging_sec_ht}
region exit_neg_vis intersect 2 exit_v_neg exit_bounds
region exit_pos_vis intersect 2 exit_v_pos exit_bounds
region exit_neg intersect 2 exit_v_neg exit_bounds side out
region exit_pos intersect 2 exit_v_pos exit_bounds side out
fix wall_exit_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region exit_neg
fix wall_exit_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region exit_pos

# Orifice cover
region orifice_cover_vis block -${x_half} ${x_half} -${orifice_half} ${orifice_half} -${eps} 0
region orifice_cover block -${x_half} ${x_half} -${orifice_half} ${orifice_half} -${eps} 0 side out
fix lid all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 0.0 0.3 1 region orifice_cover