# use _vis regions for visualising geometry

# Entry funnel (Top) - Funneling in the X-direction
region entry_f_neg plane -${x_half} 0 ${z_bin_top} -${sin_f} 0 -${cos_f}
region entry_f_pos plane  ${x_half} 0 ${z_bin_top}  ${sin_f} 0 -${cos_f}
region funnel_bounds block -${f_x_half} ${f_x_half} -${y_half} ${y_half} ${z_bin_top} ${funnel_top}
region funnel_neg_vis intersect 2 entry_f_neg funnel_bounds
region funnel_pos_vis intersect 2 entry_f_pos funnel_bounds
region funnel_neg intersect 2 entry_f_neg funnel_bounds side out
region funnel_pos intersect 2 entry_f_pos funnel_bounds side out
fix wall_funnel_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region funnel_neg
fix wall_funnel_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region funnel_pos
# Funnel side walls
region funnel_left_vis block -${f_x_half} ${f_x_half} ${y_half} ${h_y_half} ${z_bin_top} ${funnel_top}
region funnel_right_vis block -${f_x_half} ${f_x_half} -${h_y_half} -${y_half} ${z_bin_top} ${funnel_top}
region funnel_left block -${f_x_half} ${f_x_half} ${y_half} ${h_y_half} ${z_bin_top} ${funnel_top} side out
region funnel_right block -${f_x_half} ${f_x_half} -${h_y_half} -${y_half} ${z_bin_top} ${funnel_top} side out
fix wall_funnel_left all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region funnel_left
fix wall_funnel_right all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region funnel_right
