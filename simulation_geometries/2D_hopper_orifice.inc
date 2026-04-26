# Orifice cover
region orifice_cover_vis block -${x_half} ${x_half} -${orifice_half} ${orifice_half} -${eps} 0
region orifice_cover block -${x_half} ${x_half} -${orifice_half} ${orifice_half} -${eps} 0 side out
fix lid all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region orifice_cover
