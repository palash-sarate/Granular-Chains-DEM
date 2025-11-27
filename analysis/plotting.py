import numpy as np

import matplotlib.pyplot as plt

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