import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def validate_chain_spacing(chain_data, min_spacing=0.0, max_spacing=0.0025):
    """
    Validates spacing between consecutive atoms in the chain.
    Returns (is_valid, messages)
    """
    if not chain_data or len(chain_data) < 2:
        return True, []

    messages = []
    is_valid = True
    
    # Extract just x,y,z
    coords = [np.array(p[:3]) for p in chain_data]
    
    for i in range(len(coords) - 1):
        dist = np.linalg.norm(coords[i+1] - coords[i])
        if dist < min_spacing or dist > max_spacing:
            is_valid = False
            messages.append(f"Gap {i+1}-{i+2}: {dist:.6f} m (Expected {min_spacing}-{max_spacing})")
            
    return is_valid, messages

# function to visualize chain.data
def plot_chain_data(chain_data_source, title=None, save_path=None):
    """
    Plots the 3D coordinates of the chain data.
    chain_data_source: list of tuples/lists with (x, y, z) coordinates OR path to LAMMPS data file
    """
    chain_data = []
    default_title = "Chain Data Visualization"
    
    if isinstance(chain_data_source, (str, Path)):
        filepath = Path(chain_data_source)
        if not filepath.exists():
            print(f"File not found: {filepath}")
            return

        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            if lines:
                # Use the first line as title if it's not empty
                file_title = lines[0].strip().lstrip('#').strip()
                if file_title:
                    default_title = file_title

            in_atoms_section = False
            for line in lines:
                stripped = line.strip()
                if not stripped:
                    continue
                
                if stripped.startswith("Atoms"):
                    in_atoms_section = True
                    continue
                
                if in_atoms_section:
                    # Check if we hit the next section (e.g., "Bonds", "Velocities", "Masses")
                    if stripped[0].isalpha(): 
                        in_atoms_section = False
                        break
                    
                    parts = stripped.split()
                    # Assuming format: id type x y z diameter ...
                    if len(parts) >= 6:
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            d = float(parts[5])
                            chain_data.append((x, y, z, d))
                        except ValueError:
                            continue
                    elif len(parts) >= 5:
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            chain_data.append((x, y, z, 0.002)) # Default diameter
                        except ValueError:
                            continue
        except Exception as e:
            print(f"Error reading data file: {e}")
            return
    else:
        chain_data = chain_data_source

    if not chain_data:
        print("No chain data to plot.")
        return

    if title is None:
        title = default_title

    # Validate spacing
    valid, msgs = validate_chain_spacing(chain_data)
    if not valid:
        print("WARNING: Spacing validation failed!")
        for msg in msgs:
            print(f"  - {msg}")
    else:
        print("Spacing validation passed.")

    # Unpack data
    if len(chain_data[0]) == 4:
        xs, ys, zs, ds = zip(*chain_data)
    else:
        xs, ys, zs = zip(*chain_data)
        ds = [0.002] * len(xs)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot spheres
    for x, y, z, d in zip(xs, ys, zs, ds):
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 20)
        r = d / 2.0
        
        x_sphere = r * np.outer(np.cos(u), np.sin(v)) + x
        y_sphere = r * np.outer(np.sin(u), np.sin(v)) + y
        z_sphere = r * np.outer(np.ones(np.size(u)), np.cos(v)) + z
        
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color='b', alpha=0.6)

    # Also plot the line connecting centers for clarity
    ax.plot(xs, ys, zs, color='k', linewidth=1, alpha=0.5)

    # Set equal aspect ratio
    # Matplotlib 3D doesn't have 'equal' aspect ratio built-in easily, 
    # so we manually set limits to be cubic
    range_x = max(xs) - min(xs)
    range_y = max(ys) - min(ys)
    range_z = max(zs) - min(zs)
    
    max_range = np.array([range_x, range_y, range_z]).max() / 2.0
    
    # Handle case where all points are identical or collinear along axis
    if max_range == 0:
        max_range = 0.001 # Default small range

    mid_x = (max(xs)+min(xs)) * 0.5
    mid_y = (max(ys)+min(ys)) * 0.5
    mid_z = (max(zs)+min(zs)) * 0.5
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Force the box aspect ratio to be equal (cubic)
    try:
        ax.set_box_aspect([1, 1, 1])
    except AttributeError:
        # Older matplotlib versions might not support this
        pass

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title(title)

    if save_path:
        plt.savefig(save_path)
    # plt.show()

def plot_angle_evolution(angle_data_list, legend_labels=None, title="Angle Evolution", save_path=None):
    """
    Plots multiple series of angle vs time.
    angle_data_list: list of lists, where each inner list contains tuples (timestep, angle)
    legend_labels: optional list of strings for the legend
    """
    if not angle_data_list:
        print("No data to plot.")
        return

    plt.figure(figsize=(10, 6))
    
    # If legend_labels is not provided or doesn't match length, generate default labels
    if legend_labels is None or len(legend_labels) != len(angle_data_list):
        legend_labels = [f'Series {i+1}' for i in range(len(angle_data_list))]

    colors = plt.cm.viridis(np.linspace(0, 1, len(angle_data_list)))

    for i, data in enumerate(angle_data_list):
        if not data:
            continue
        steps, values = zip(*data)
        plt.plot(steps, values, label=legend_labels[i], color=colors[i])

    plt.xlabel('Timestep')
    plt.ylabel('Angle (degrees)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    if save_path:
        plt.savefig(save_path)
    # plt.show()
    
def plot_distance_evolution(distance_data_list, legend_labels=None, title="Distance Evolution", save_path=None):
    """
    Plots multiple series of distance vs time.
    distance_data_list: list of lists, where each inner list contains tuples (timestep, distance)
    legend_labels: optional list of strings for the legend
    """
    if not distance_data_list:
        print("No data to plot.")
        return

    plt.figure(figsize=(10, 6))
    
    # If legend_labels is not provided or doesn't match length, generate default labels
    if legend_labels is None or len(legend_labels) != len(distance_data_list):
        legend_labels = [f'Series {i+1}' for i in range(len(distance_data_list))]

    colors = plt.cm.viridis(np.linspace(0, 1, len(distance_data_list)))

    for i, data in enumerate(distance_data_list):
        if not data:
            continue
        steps, values = zip(*data)
        plt.plot(steps, values, label=legend_labels[i], color=colors[i])

    plt.xlabel('Timestep')
    plt.ylabel('Distance (m)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    if save_path:
        plt.savefig(save_path)
    
def plot_xyz_evolution(xyz_data_list, legend_labels=None, title="X, Y, Z Evolution", save_path=None):
    """
    Plots multiple series of x, y, z coordinates vs time.
    xyz_data_list: list of lists, where each inner list contains tuples (timestep, x, y, z)
    legend_labels: optional list of strings for the legend prefix of each series
    """
    if not xyz_data_list:
        print("No data to plot.")
        return

    plt.figure(figsize=(10, 6))

    # If legend_labels is not provided or doesn't match length, generate default labels
    if legend_labels is None or len(legend_labels) != len(xyz_data_list):
        legend_labels = [f'Series {i+1}' for i in range(len(xyz_data_list))]

    # Use different line styles or markers to distinguish series if needed, 
    # but here we will rely on color/label combinations or just standard colors.
    # To keep it readable, we might cycle colors or line styles.
    linestyles = ['-', '--', '-.', ':']
    
    for i, data in enumerate(xyz_data_list):
        if not data:
            continue
        
        steps, x_vals, y_vals, z_vals = zip(*data)
        ls = linestyles[i % len(linestyles)]
        prefix = legend_labels[i]
        
        plt.plot(steps, x_vals, label=f'{prefix} - X', color='r', linestyle=ls)
        plt.plot(steps, y_vals, label=f'{prefix} - Y', color='g', linestyle=ls)
        plt.plot(steps, z_vals, label=f'{prefix} - Z', color='b', linestyle=ls)

    plt.xlabel('Timestep')
    plt.ylabel('Coordinates (units)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    if save_path:
        plt.savefig(save_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize LAMMPS chain data file.")
    parser.add_argument("file", type=str, help="Path to the LAMMPS data file")
    parser.add_argument("--title", type=str, help="Title of the plot", default=None)
    parser.add_argument("--save", type=str, help="Path to save the plot image", default=None)
    
    args = parser.parse_args()
    
    plot_chain_data(args.file, title=args.title, save_path=args.save)
    if not args.save:
        plt.show()