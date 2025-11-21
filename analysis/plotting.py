import matplotlib.pyplot as plt

def plot_angle_evolution(angle_data, title="Angle Evolution"):
    """
    Plots angle vs time.
    angle_data: list of tuples (timestep, angle)
    """
    if not angle_data:
        print("No data to plot.")
        return

    steps, values = zip(*angle_data)
    
    plt.figure(figsize=(10, 6))
    plt.plot(steps, values, label='Angle', color='b')
    plt.xlabel('Timestep')
    plt.ylabel('Angle (degrees)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()
    
def plot_distance_evolution(distance_data, title="Distance Evolution"):
    """
    Plots distance vs time.
    distance_data: list of tuples (timestep, distance)
    """
    if not distance_data:
        print("No data to plot.")
        return

    steps, values = zip(*distance_data)
    
    plt.figure(figsize=(10, 6))
    plt.plot(steps, values, label='Distance', color='r')
    plt.xlabel('Timestep')
    plt.ylabel('Distance (units)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()
    
def plot_xyz_evolution(xyz_data, title="X, Y, Z Evolution"):
    """
    Plots x, y, z coordinates vs time.
    xyz_data: list of tuples (timestep, x, y, z)
    """
    if not xyz_data:
        print("No data to plot.")
        return

    steps, x_vals, y_vals, z_vals = zip(*xyz_data)
    
    plt.figure(figsize=(10, 6))
    plt.plot(steps, x_vals, label='X', color='r')
    plt.plot(steps, y_vals, label='Y', color='g')
    plt.plot(steps, z_vals, label='Z', color='b')
    plt.xlabel('Timestep')
    plt.ylabel('Coordinates (units)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()