# --- 2D Hopper Stubs Array Geometry ---
# Defines a lattice of cylindrical stubs parallel to the x-axis.

# Parameters
variable stub_radius  equal 0.002    # 2 mm
variable stub_spacing equal 0.028    # 28 mm
variable stub_top_offset equal 0.01   # Vertical offset from funnel_top

# --- Calculations ---
# Y direction (Horizontal): ranges from -y_half to y_half
variable ny equal floor(2.0*${y_half}/${stub_spacing})+1
variable y_span equal (${ny}-1)*${stub_spacing}
variable y_start equal -${y_half}+(2.0*${y_half}-${y_span})/2.0

# Z direction (Vertical): concentrated in the top half (closer to funnel_top)
variable z_target_min equal ${z_bin_top} + (${funnel_top} - ${z_bin_top}) / 2.0
variable z_target_max equal ${funnel_top} - ${stub_top_offset}
variable nz equal floor((${z_target_max} - ${z_target_min}) / ${stub_spacing}) + 1
variable z_span equal (${nz} - 1) * ${stub_spacing}
variable z_start equal ${z_target_max} - ${z_span}

# --- Loop Generation ---
# Note: This creates many regions and fixes. 
# Ensure your LAMMPS version/configuration allows this many fixes (~180).

label loop_z
variable iz loop ${nz}
    variable cz equal ${z_start}+(${iz}-1)*${stub_spacing}
    
    label loop_y
    variable iy loop ${ny}
        variable cy equal ${y_start}+(${iy}-1)*${stub_spacing}
        
        # Region name: s_Z_Y
        # Cylinder parallel to X-axis
        region s_${iz}_${iy} cylinder x ${cy} ${cz} ${stub_radius} -${f_x_half} ${f_x_half} side out
        region s_${iz}_${iy}_vis cylinder x ${cy} ${cz} ${stub_radius} -${f_x_half} ${f_x_half}
        print "Region ${iz} ${iy} created"
        # Physical wall fix
        fix f_s_${iz}_${iy} all wall/gran/region hertz/history 1.0e6 1.0e6 5e-4 5e-4 0.3 1 region s_${iz}_${iy}
        
        next iy
        jump SELF loop_y
    variable iy delete # Clean up inner loop variable
    
    next iz
    jump SELF loop_z
variable iz delete # Clean up outer loop variable

# --- Cleanup temporary calculation variables ---
print "Stubs Array Generated: ${ny} columns x ${nz} rows"
variable y_start delete
variable y_span delete
variable z_start delete
variable z_span delete


# vis command
# 1. Generate the updated list for the top 2 rows with correct suffixing
# $regions = (1..2 | % { $iz=$_; 1..12 | % { "s_${iz}_$($_)_vis" } })

# # 2. Run the extractor
# python analysis/geometry_extractor.py simulation_geometries/2D_hopper_stubs_array.inc `
#     --regions $regions `
#     --var y_half=0.155 `
#     --var z_bin_top=0.4 `
#     --var funnel_top=0.5 `
#     --var f_x_half=0.1 `
#     --spacing 0.0005
