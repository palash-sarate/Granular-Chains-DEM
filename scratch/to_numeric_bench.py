import os
import time
import numpy as np
import pandas as pd

def parse_without_dtype(filepath, count, cols, data_start_line):
    df = pd.read_csv(filepath, 
                     skiprows=data_start_line,
                     names=cols, 
                     nrows=count if count > 0 else None,
                     sep=r'\s+', 
                     engine='c',
                     memory_map=True)
    return df

def parse_with_dtype(filepath, count, cols, data_start_line):
    # Specify dtypes
    dtype_dict = {}
    for col in cols:
        if col in ['id', 'mol', 'type', 'index']:
            dtype_dict[col] = np.int32
        else:
            dtype_dict[col] = np.float32 # use float32 to save memory and parse faster!
            
    df = pd.read_csv(filepath, 
                     skiprows=data_start_line,
                     names=cols, 
                     nrows=count if count > 0 else None,
                     sep=r'\s+', 
                     engine='c',
                     dtype=dtype_dict,
                     memory_map=True)
    return df

def run_bench():
    large_file = 'scratch/large_synthetic.dump'
    
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
            
    print(f"Benchmarking with 200,000 rows...")
    
    # Without dtype
    t0 = time.time()
    for _ in range(5):
        df1 = parse_without_dtype(large_file, count, cols, data_start_line)
    t1 = time.time()
    time_without = (t1 - t0) / 5
    print(f"Without dtype hints: {time_without*1000:.2f} ms")
    
    # With dtype
    t0 = time.time()
    for _ in range(5):
        df2 = parse_with_dtype(large_file, count, cols, data_start_line)
    t1 = time.time()
    time_with = (t1 - t0) / 5
    print(f"With dtype hints (float32): {time_with*1000:.2f} ms")
    
    print(f"Speedup from specifying dtypes: {time_without / time_with:.2f}x")

if __name__ == '__main__':
    run_bench()
