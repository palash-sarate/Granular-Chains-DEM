from pathlib import Path
from simulation.chain_generator import ChainConfig, write_chain_data

def generate_chains():
    # Example 1: Basic linear chain
    print("Generating linear chain...")
    linear_config = ChainConfig(
        beads=50,
        spacing=0.0025,
        mode="linear",
        orientation="vert",
        output_dir=Path("generated_chains")
    )
    path1 = write_chain_data(linear_config)
    print(f"Created: {path1}")

    # Example 2: Chain in a specific direction
    print("\nGenerating directional chain...")
    direction_config = ChainConfig(
        beads=30,
        spacing=0.0025,
        mode="direction",
        direction=(1.0, 1.0, 0.0),  # Diagonal in XY plane
        output_dir=Path("generated_chains")
    )
    path2 = write_chain_data(direction_config)
    print(f"Created: {path2}")

    # Example 3: Random walk chain
    print("\nGenerating random walk chain...")
    random_config = ChainConfig(
        beads=100,
        spacing=0.0025,
        mode="random",
        max_spacing=0.003,  # Optional: variable spacing
        output_dir=Path("generated_chains")
    )
    path3 = write_chain_data(random_config)
    print(f"Created: {path3}")

    # Example 4: Loop chain
    print("\nGenerating loop chain...")
    loop_config = ChainConfig(
        beads=80,
        spacing=0.0025,
        mode="loop",
        loop_radius=0.01,
        output_dir=Path("generated_chains")
    )
    path4 = write_chain_data(loop_config)
    print(f"Created: {path4}")

if __name__ == "__main__":
    generate_chains()
