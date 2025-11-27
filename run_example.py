from simulation import SimulationConfig, SimulationRunner

# Create a config with a specific run name
config = SimulationConfig(
    run_name="Viscosity_03",
    base_output_dir="post_chain_flop",
    run_steps=10000  # Short run for testing
)

# Initialize runner
# Ensure 'lmp' is in your PATH or provide absolute path
runner = SimulationRunner(lammps_executable="lmp")

print(f"Running simulation: {config.run_name}")
print(f"Output directory: {config.output_dir}")

# Run simulation
# This will create directories: post_chain_flop/Viscosity_03/{bond,angle,restart}
# and generate a temporary input script 'generated_in.chain_flop'
runner.run(config, template_path="in.chain_flop")
