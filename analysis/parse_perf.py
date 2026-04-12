import argparse
import re
import matplotlib.pyplot as plt
import os
import sys

def parse_lammps_log(log_file_path):
    """
    Parses a LAMMPS log file for loop time and atom counts.
    Returns a list of (atoms, loop_time) tuples.
    """
    data = []
    # Pattern: Loop time of 1064.48 on 8 procs for 69949 steps with 5256 atoms
    pattern = re.compile(r"Loop time of ([\d.]+) on \d+ procs for \d+ steps with (\d+) atoms")
    
    try:
        with open(log_file_path, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    loop_time = float(match.group(1))
                    atoms = int(match.group(2))
                    data.append((atoms, loop_time))
    except FileNotFoundError:
        print(f"Error: File '{log_file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    return data

def plot_performance(data, output_path):
    """
    Plots loop time vs number of atoms.
    """
    if not data:
        print("No matching data found in the log file.")
        return

    # Sort data by number of atoms for a clean line plot
    data.sort(key=lambda x: x[0])
    atoms, times = zip(*data)

    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

    # Use a vibrant color for the line and points
    ax.plot(atoms, times, marker='o', linestyle='-', color='#00f2ff', label='Scan Data', linewidth=2, markersize=8, alpha=0.8)
    
    # Add labels and title with premium styling
    ax.set_xlabel('Number of Atoms', fontsize=12, fontweight='bold', color='#eeeeee')
    ax.set_ylabel('Loop Time (s)', fontsize=12, fontweight='bold', color='#eeeeee')
    ax.set_title('LAMMPS Performance: Loop Time vs. Number of Atoms', fontsize=14, fontweight='bold', pad=20, color='#ffffff')
    
    # Grid styling
    ax.grid(True, linestyle='--', alpha=0.3, color='#555555')
    
    # Customize ticks
    ax.tick_params(colors='#cccccc', labelsize=10)
    
    # Add a legend
    ax.legend(facecolor='#222222', edgecolor='#444444', labelcolor='#ffffff')

    # Tight layout for better spacing
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path)
    print(f"Plot saved to: {output_path}")
    
    # Optionally show the plot (might not work in all environments, but good to have)
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description="Parse LAMMPS log for loop times and plot vs atom count.")
    parser.add_argument("log_file", help="Path to the LAMMPS log file")
    parser.add_argument("--output", "-o", default="performance_plot.png", help="Path to save the output plot (default: performance_plot.png)")
    
    args = parser.parse_args()
    
    print(f"Parsing log file: {args.log_file}")
    data = parse_lammps_log(args.log_file)
    
    if data:
        print(f"Found {len(data)} data points.")
        # If output path is just a filename, save it in the same dir as the log file
        if not os.path.isabs(args.output) and os.path.dirname(args.output) == '':
            log_dir = os.path.dirname(os.path.abspath(args.log_file))
            output_path = os.path.join(log_dir, args.output)
        else:
            output_path = args.output
            
        plot_performance(data, output_path)
    else:
        print("No performance data found in the provided log file.")

if __name__ == "__main__":
    main()
