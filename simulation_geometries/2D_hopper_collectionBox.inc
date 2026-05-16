# --- Collection Box Geometry Definitions ---
# Dimensions: 2y_half * 2y_half * 2y_half
# Located right below the hopper (z from -2*y_half to 0)

# Bottom floor of the collection box
region cb_floor_vis block -${x_half} ${x_half} -${y_half} ${y_half} -2*${y_half}-${eps} -2*${y_half}
region cb_floor block -${x_half} ${x_half} -${y_half} ${y_half} -2*${y_half}-${eps} -2*${y_half} side out
fix wall_cb_floor all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region cb_floor

# Left and Right side walls (y-direction)
region cb_wall_pos_vis block -${x_half} ${x_half} ${y_half} ${y_half}+${eps} -2*${y_half} 0
region cb_wall_pos block -${x_half} ${x_half} ${y_half} ${y_half}+${eps} -2*${y_half} 0 side out
fix wall_cb_wall_pos all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region cb_wall_pos

region cb_wall_neg_vis block -${x_half} ${x_half} -${y_half}-${eps} -${y_half} -2*${y_half} 0
region cb_wall_neg block -${x_half} ${x_half} -${y_half}-${eps} -${y_half} -2*${y_half} 0 side out
fix wall_cb_wall_neg all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region cb_wall_neg

# Front and Back covers (x-direction)
region cb_front_vis block ${x_half} ${x_half}+${eps} -${y_half} ${y_half} -2*${y_half} 0
region cb_front block ${x_half} ${x_half}+${eps} -${y_half} ${y_half} -2*${y_half} 0 side out
fix wall_cb_front all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region cb_front

region cb_back_vis block -${x_half}-${eps} -${x_half} -${y_half} ${y_half} -2*${y_half} 0
region cb_back block -${x_half}-${eps} -${x_half} -${y_half} ${y_half} -2*${y_half} 0 side out
fix wall_cb_back all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region cb_back
