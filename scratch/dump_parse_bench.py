import os
import time
import numpy as np
import pandas as pd

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
    
    # Fast parse floats using numpy.fromstring
    arr = np.fromstring(data_str, dtype=float, sep=' ')
    arr = arr.reshape((-1, len(cols)))
    
    # Create DataFrame from array
    df = pd.DataFrame(arr, columns=cols)
    return df

def run_bench():
    d = 'dumping_yard/Hopper_Fill/Grid_Fill_1H_S808086'
    chain_file = os.path.join(d, 'chain/chain_100000.dump')
    
    # Read headers to get count and cols
    with open(chain_file, 'r') as f:
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
            
    print(f"Benchmarking parser on {chain_file} ({count} atoms, {len(cols)} columns)...")
    
    # Warmup
    _ = parse_pandas_current(chain_file, count, cols, data_start_line)
    
    # Method 1: Current pandas
    t0 = time.time()
    for _ in range(10):
        df1 = parse_pandas_current(chain_file, count, cols, data_start_line)
    t1 = time.time()
    time_curr = (t1 - t0) / 10
    print(f"Current Pandas parser: {time_curr*1000:.2f} ms")
    
    # Method 3: NumPy fromstring
    t0 = time.time()
    for _ in range(10):
        df3 = parse_numpy_fromstring(chain_file, cols)
    t1 = time.time()
    time_np = (t1 - t0) / 10
    print(f"NumPy fromstring: {time_np*1000:.2f} ms")
    
    # Verify dataframes are equal
    df1_comp = df1.astype(float)
    df3_comp = df3.astype(float)
    max_diff = np.max(np.abs(df1_comp.values - df3_comp.values))
    print(f"Verification: max difference between current and NumPy parse: {max_diff:.2e}")

if __name__ == '__main__':
    run_bench()
