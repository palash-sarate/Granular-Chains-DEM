from simulation.runner import SimulationRunner
runner = SimulationRunner()
features = runner._get_lammps_features()
print(f"Features: {features}")
print(f"OPENMP in features: {'OPENMP' in features}")
