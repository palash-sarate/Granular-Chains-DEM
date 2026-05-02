# Benchmark specialized dump/thermo for Grid Hopper Filling
# Optimized for performance tracking with minimal I/O overhead

thermo 1000
thermo_style custom step time cpu pe ke etotal

# Disable heavy dumps for benchmarking purposes
# If you need to see the result, use a standard dump instead
