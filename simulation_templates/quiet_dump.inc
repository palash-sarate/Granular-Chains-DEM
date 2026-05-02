# --- Optimized Quiet Dump for Relaxation Library Generation ---
# Reduces I/O overhead by using single files and larger intervals.

# Re-define computes
compute 1 chain_particles bond/local dist engpot force 
compute 2 chain_particles angle/local theta eng

thermo 50000
thermo_style custom step time cpu pe ke ebond eangle etotal

# Dump to a single trajectory file instead of many small files
# Using .lammpstrj suffix for easy identification
dump 1 all custom 50000 ${outdir}/chain/trajectory.lammpstrj id mol type x y z vx vy vz fx fy fz diameter mass
dump_modify 1 time yes
# dump_modify 1 sort id

# Single files for local data
# dump 2 all local 50000 ${outdir}/bond/bonds.dump index c_1[1] c_1[2] c_1[3]
# dump 3 all local 50000 ${outdir}/angle/angles.dump index c_2[1] c_2[2]
