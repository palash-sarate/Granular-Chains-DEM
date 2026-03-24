# Static hopper bodies (no move for filling)
region hopper_cone cone z 0.0 0.0 0.01 0.05 0.05 0.15 open 1 open 2
region hopper_cyl cylinder z 0.0 0.0 0.01 0.02 0.05 open 1 open 2
region hopper_union union 2 hopper_cone hopper_cyl
fix wall all wall/gran/region hertz/history 1.0e8 1.0e8 5e-4 0.0 0.3 1 region hopper_union