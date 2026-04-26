# --- Global Variables (Resolved in Python) ---
variable x_half equal 0.00125
variable y_half equal 0.155
variable orifice_half equal 0.05
variable z_bin_top equal 0.4
variable funnel_ang equal 45.0
variable sin_f equal 0.7071067805519557
variable cos_f equal 0.7071067818211394
variable hopper_ang equal 60.0
variable sin_c equal 0.8660254031861399
variable cos_c equal 0.500000001036284
variable converging_sec_ht equal 0.18186533429216112
variable eps equal 0.01
variable funnel_top equal 0.5
variable f_x_half equal 0.1
variable h_y_half equal 0.165
variable h_x_half equal 0.01125

# --- Replicated Hopper 0 at [0. 0. 0.] ---
region hopper_pos_vis_0 block -0.00125 0.00125 0.155 0.165 0.18186533429216112 0.4 
region hopper_neg_vis_0 block -0.00125 0.00125 -0.165 -0.155 0.18186533429216112 0.4 
region hopper_pos_0 block -0.00125 0.00125 0.155 0.165 0.18186533429216112 0.4 side out
region hopper_neg_0 block -0.00125 0.00125 -0.165 -0.155 0.18186533429216112 0.4 side out
fix wall_hopper_pos_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_pos_0
fix wall_hopper_neg_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_neg_0
region hopper_front_vis_0 block 0.00125 0.01125 -0.155 0.155 0.0 0.4 
region hopper_back_vis_0 block -0.01125 -0.00125 -0.155 0.155 0.0 0.4 
region hopper_front_0 block 0.00125 0.01125 -0.155 0.155 0.0 0.4 side out
region hopper_back_0 block -0.01125 -0.00125 -0.155 0.155 0.0 0.4 side out
fix wall_hopper_front_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_front_0
fix wall_hopper_back_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_back_0
region exit_v_neg_0 plane 0.0 -0.05 0.0 0 -${sin_c} -${cos_c}
region exit_v_pos_0 plane 0.0 0.05 0.0 0 ${sin_c} -${cos_c}
region exit_bounds_0 block -0.00125 0.00125 -0.155 0.155 0.0 0.18186533429216112 
region exit_neg_vis_0 intersect 2 exit_v_neg_0 exit_bounds_0
region exit_pos_vis_0 intersect 2 exit_v_pos_0 exit_bounds_0
region exit_neg_0 intersect 2 exit_v_neg_0 exit_bounds_0 side out
region exit_pos_0 intersect 2 exit_v_pos_0 exit_bounds_0 side out
fix wall_exit_neg_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_neg_0
fix wall_exit_pos_0 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_pos_0

# --- Replicated Hopper 1 at [2. 0. 0.] ---
region hopper_pos_vis_1 block 1.99875 2.00125 0.155 0.165 0.18186533429216112 0.4 
region hopper_neg_vis_1 block 1.99875 2.00125 -0.165 -0.155 0.18186533429216112 0.4 
region hopper_pos_1 block 1.99875 2.00125 0.155 0.165 0.18186533429216112 0.4 side out
region hopper_neg_1 block 1.99875 2.00125 -0.165 -0.155 0.18186533429216112 0.4 side out
fix wall_hopper_pos_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_pos_1
fix wall_hopper_neg_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_neg_1
region hopper_front_vis_1 block 2.00125 2.01125 -0.155 0.155 0.0 0.4 
region hopper_back_vis_1 block 1.98875 1.99875 -0.155 0.155 0.0 0.4 
region hopper_front_1 block 2.00125 2.01125 -0.155 0.155 0.0 0.4 side out
region hopper_back_1 block 1.98875 1.99875 -0.155 0.155 0.0 0.4 side out
fix wall_hopper_front_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_front_1
fix wall_hopper_back_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_back_1
region exit_v_neg_1 plane 2.0 -0.05 0.0 0 -${sin_c} -${cos_c}
region exit_v_pos_1 plane 2.0 0.05 0.0 0 ${sin_c} -${cos_c}
region exit_bounds_1 block 1.99875 2.00125 -0.155 0.155 0.0 0.18186533429216112 
region exit_neg_vis_1 intersect 2 exit_v_neg_1 exit_bounds_1
region exit_pos_vis_1 intersect 2 exit_v_pos_1 exit_bounds_1
region exit_neg_1 intersect 2 exit_v_neg_1 exit_bounds_1 side out
region exit_pos_1 intersect 2 exit_v_pos_1 exit_bounds_1 side out
fix wall_exit_neg_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_neg_1
fix wall_exit_pos_1 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_pos_1

# --- Replicated Hopper 2 at [0. 2. 0.] ---
region hopper_pos_vis_2 block -0.00125 0.00125 2.155 2.165 0.18186533429216112 0.4 
region hopper_neg_vis_2 block -0.00125 0.00125 1.835 1.845 0.18186533429216112 0.4 
region hopper_pos_2 block -0.00125 0.00125 2.155 2.165 0.18186533429216112 0.4 side out
region hopper_neg_2 block -0.00125 0.00125 1.835 1.845 0.18186533429216112 0.4 side out
fix wall_hopper_pos_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_pos_2
fix wall_hopper_neg_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_neg_2
region hopper_front_vis_2 block 0.00125 0.01125 1.845 2.155 0.0 0.4 
region hopper_back_vis_2 block -0.01125 -0.00125 1.845 2.155 0.0 0.4 
region hopper_front_2 block 0.00125 0.01125 1.845 2.155 0.0 0.4 side out
region hopper_back_2 block -0.01125 -0.00125 1.845 2.155 0.0 0.4 side out
fix wall_hopper_front_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_front_2
fix wall_hopper_back_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_back_2
region exit_v_neg_2 plane 0.0 1.95 0.0 0 -${sin_c} -${cos_c}
region exit_v_pos_2 plane 0.0 2.05 0.0 0 ${sin_c} -${cos_c}
region exit_bounds_2 block -0.00125 0.00125 1.845 2.155 0.0 0.18186533429216112 
region exit_neg_vis_2 intersect 2 exit_v_neg_2 exit_bounds_2
region exit_pos_vis_2 intersect 2 exit_v_pos_2 exit_bounds_2
region exit_neg_2 intersect 2 exit_v_neg_2 exit_bounds_2 side out
region exit_pos_2 intersect 2 exit_v_pos_2 exit_bounds_2 side out
fix wall_exit_neg_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_neg_2
fix wall_exit_pos_2 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_pos_2

# --- Replicated Hopper 3 at [2. 2. 0.] ---
region hopper_pos_vis_3 block 1.99875 2.00125 2.155 2.165 0.18186533429216112 0.4 
region hopper_neg_vis_3 block 1.99875 2.00125 1.835 1.845 0.18186533429216112 0.4 
region hopper_pos_3 block 1.99875 2.00125 2.155 2.165 0.18186533429216112 0.4 side out
region hopper_neg_3 block 1.99875 2.00125 1.835 1.845 0.18186533429216112 0.4 side out
fix wall_hopper_pos_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_pos_3
fix wall_hopper_neg_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_neg_3
region hopper_front_vis_3 block 2.00125 2.01125 1.845 2.155 0.0 0.4 
region hopper_back_vis_3 block 1.98875 1.99875 1.845 2.155 0.0 0.4 
region hopper_front_3 block 2.00125 2.01125 1.845 2.155 0.0 0.4 side out
region hopper_back_3 block 1.98875 1.99875 1.845 2.155 0.0 0.4 side out
fix wall_hopper_front_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_front_3
fix wall_hopper_back_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region hopper_back_3
region exit_v_neg_3 plane 2.0 1.95 0.0 0 -${sin_c} -${cos_c}
region exit_v_pos_3 plane 2.0 2.05 0.0 0 ${sin_c} -${cos_c}
region exit_bounds_3 block 1.99875 2.00125 1.845 2.155 0.0 0.18186533429216112 
region exit_neg_vis_3 intersect 2 exit_v_neg_3 exit_bounds_3
region exit_pos_vis_3 intersect 2 exit_v_pos_3 exit_bounds_3
region exit_neg_3 intersect 2 exit_v_neg_3 exit_bounds_3 side out
region exit_pos_3 intersect 2 exit_v_pos_3 exit_bounds_3 side out
fix wall_exit_neg_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_neg_3
fix wall_exit_pos_3 all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region exit_pos_3
