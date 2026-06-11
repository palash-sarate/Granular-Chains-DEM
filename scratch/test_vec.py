import os
import pandas as pd
import numpy as np
from analysis.data_manager import _parse_single_dump_fast

def calculate_angle(p1, p2, p3):
    v21 = p1 - p2
    v23 = p3 - p2
    cos_t = np.dot(v21, v23) / (np.linalg.norm(v21) * np.linalg.norm(v23))
    return np.arccos(np.clip(cos_t, -1.0, 1.0)) * 180.0 / np.pi

def run_test():
    d = 'dumping_yard/Hopper_Fill/Grid_Fill_1H_S808086'
    
    df_atoms = _parse_single_dump_fast(os.path.join(d, 'chain/chain_100000.dump'))
    df_angles_ref = _parse_single_dump_fast(os.path.join(d, 'angle/angle_100000.dump'))
    
    # Let's map atom ID to coordinates
    coords = {}
    for idx, row in df_atoms.iterrows():
        coords[int(row['id'])] = row[['x', 'y', 'z']].values.astype(float)
        
    # Angle definitions from final_grid.data
    defs = {
        6535: (13193, 13194, 13195),
        6536: (13125, 13126, 13127),
        6537: (13150, 13151, 13152),
        6538: (13250, 13251, 13252),
        6539: (13249, 13250, 13251),
        6540: (13190, 13191, 13192),
        6541: (13149, 13150, 13151),
        6542: (13189, 13190, 13191),
        6543: (13186, 13187, 13188),
        6544: (13185, 13186, 13187),
        6545: (13182, 13183, 13184)
    }
    
    print("Angle ID | Atoms | Calculated | Dumped")
    print("---------------------------------------")
    for aid, (a1, a2, a3) in defs.items():
        p1 = coords[a1]
        p2 = coords[a2]
        p3 = coords[a3]
        calc = calculate_angle(p1, p2, p3)
        dump_row = df_angles_ref[df_angles_ref['index'] == aid]
        dump = dump_row['theta'].iloc[0] if not dump_row.empty else 'N/A'
        print(f"{aid:8d} | {a1},{a2},{a3} | {calc:10.3f} | {dump}")

if __name__ == '__main__':
    run_test()
