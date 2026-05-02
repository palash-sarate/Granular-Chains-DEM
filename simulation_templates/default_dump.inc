# Re-define computes
compute 1 chain_particles bond/local dist engpot force 
compute 2 chain_particles angle/local theta eng

thermo 1000
thermo_style custom step time cpu pe ke ebond eangle etotal ecoul emol edihed epair evdwl

# Dumps (Note: 'append yes' is not strictly needed if filenames have step numbers, 
# but good if you use single filenames)
dump 1 all custom 1000 ${outdir}/chain/chain_*.dump id mol type x y z vx vy vz fx fy fz diameter mass
dump_modify 1 sort id time yes
dump 2 all local 1000 ${outdir}/bond/bond_*.dump index c_1[1] c_1[2] c_1[3]
dump_modify 2 colname 2 dist
dump_modify 2 colname 3 energy
dump_modify 2 colname 4 force

dump 3 all local 1000 ${outdir}/angle/angle_*.dump index c_2[1] c_2[2]
dump_modify 3 colname 2 theta 
dump_modify 3 colname 3 energy