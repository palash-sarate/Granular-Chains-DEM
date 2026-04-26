import subprocess
from pathlib import Path

# 1. Path settings
inc_file = "temp/replicated_geometry.inc"
extractor = "analysis/geometry_extractor.py"

if not Path(inc_file).exists():
    print(f"Error: {inc_file} not found. Run the replication test first.")
    exit(1)

# 2. Automatically find all regions with '_vis' in their name
with open(inc_file, 'r') as f:
    regions = [line.split()[1] for line in f if line.strip().startswith("region") and "_vis" in line]

print(f"Found {len(regions)} visualization regions.")

# 3. Define the sampling bounds to cover the 2x2 grid (spacing 2.0 + padding)
# Format: xlo xhi ylo yhi zlo zhi
bounds = ["-0.5", "5.0", "-0.5", "5.0", "-0.1", "0.6"]

# 4. Construct and run the command
cmd = [
    "python", extractor, inc_file,
    "--regions"
] + regions + [
    "--bounds"
] + bounds + [
    "--combined",       # Combine all hoppers into a single mesh for easier viewing
    "--spacing", "0.005" # 5mm resolution for good detail
]

print("Launching Geometry Extractor...")
subprocess.run(cmd)
