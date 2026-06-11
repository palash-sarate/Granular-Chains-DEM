import os
import time
import numpy as np
import pandas as pd

def generate_large_dump(filename, num_atoms=500000):
    print(f"Generating synthetic dump file with {num_atoms} atoms...")
    # Generate random coordinates and atom properties
    np.random.seed(42)
    ids = np.arange(1, num_atoms + 1)
    mols = (ids - 1) // 4 + 1
    types = np.ones(num_atoms, dtype=int)
    coords = np.random.rand(num_atoms, 3) * 10
    vels = np.random.rand(num_atoms, 3) - 0.5
    forces = np.random.rand(num_atoms, 3) - 0.5
    diam = np.ones(num_atoms) * 0.002
    mass = np.ones(num_atoms) * 2.9e-5
    
    with open(filename, 'w') as f:
        f.write("ITEM: TIMESTEP\n100000\n")
        f.write(f"ITEM: NUMBER OF ATOMS\n{num_atoms}\n")
        f.write("ITEM: BOX BOUNDS ff ff ff\n-2 2.5\n-2 2.5\n-2 3.5\n")
        f.write("ITEM: ATOMS id mol type x y z vx vy vz fx fy fz diameter mass\n")
        
        # Write data rows
        for i in range(num_atoms):
            f.write(f"{ids[i]} {mols[i]} {types[i]} {coords[i,0]:.6f} {coords[i,1]:.6f} {coords[i,2]:.6f} "
                    f"{vels[i,0]:.6f} {vels[i,1]:.6f} {vels[i,2]:.6f} {forces[i,0]:.6f} {forces[i,1]:.6f} {forces[i,2]:.6f} "
                    f"{diam[i]:.3f} {mass[i]:.5e}\n")
    print("Generation complete.")

def parse_pandas_current(filepath, count, cols, data_start_line):
    df = pd.read_csv(filepath, 
                     skiprows=data_start_line,
                     names=cols, 
                     nrows=count if count > 0 else None,
                     sep=r'\s+', 
                     engine='c',
                     memory_map=True)
    return df

def parse_numpy_fromstring(filepath, cols):
    with open(filepath, 'r') as f:
        content = f.read()
    
    idx = content.find("ITEM: ATOMS")
    if idx == -1:
        idx = content.find("ITEM: ENTRIES")
    if idx == -1:
        return None
        
    end_line = content.find("\n", idx)
    data_str = content[end_line+1:]
    
    arr = np.fromstring(data_str, dtype=float, sep=' ')
    arr = arr.reshape((-1, len(cols)))
    
    df = pd.DataFrame(arr, columns=cols)
    return df

def run_bench():
    large_file = 'scratch/large_synthetic.dump'
    if not os.path.exists(large_file):
        generate_large_dump(large_file, 200000) # use 200,000 for standard test
        
    # Read headers
    with open(large_file, 'r') as f:
        header_lines = [f.readline() for _ in range(15)]
    
    count = 0
    data_start_line = 0
    cols = []
    for i, line in enumerate(header_lines):
        if "ITEM: NUMBER OF" in line:
            count = int(header_lines[i+1].strip())
        elif "ITEM: ATOMS" in line or "ITEM: ENTRIES" in line:
            cols = line.split()[2:]
            data_start_line = i + 1
            break
            
    print(f"\nBenchmarking on large file ({count} atoms)...")
    
    # Method 1: Current pandas
    t0 = time.time()
    for _ in range(3):
        df1 = parse_pandas_current(large_file, count, cols, data_start_line)
    t1 = time.time()
    print(f"Current Pandas parser: {(t1 - t0)/3*1000:.2f} ms")
    
    # Method 2: NumPy fromstring
    t0 = time.time()
    for _ in range(3):
        df3 = parse_numpy_fromstring(large_file, cols)
    t1 = time.time()
    print(f"NumPy fromstring: {(t1 - t0)/3*1000:.2f} ms")

if __name__ == '__main__':
    run_bench()
