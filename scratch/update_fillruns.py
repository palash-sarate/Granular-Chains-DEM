import os
import json
from pathlib import Path

runs_to_update = [
    'Grid_Fill_1H_S762224', 'Grid_Fill_1H_S762251', 'Grid_Fill_1H_S762269',
    'Grid_Fill_1H_S762071', 'Grid_Fill_1H_S762056', 'Grid_Fill_1H_S761743',
    'Grid_Fill_1H_S761319', 'Grid_Fill_1H_S782876', 'Grid_Fill_1H_S782733',
    'Grid_Fill_1H_S782860', 'Grid_Fill_1H_S782642'
]

# Set to your cluster's dumping_yard path
dumping_yard = Path('/Data/palash_data/dumping_yard')
count_updated = 0
found_runs = set()

# Scan all json files in dumping_yard
for json_file in dumping_yard.rglob('*.json'):
    if json_file.name not in ['metadata.json', 'grid_metadata.json']:
        continue
        
    is_target = False
    for r in runs_to_update:
        if r in json_file.parent.name:
            is_target = True
            found_runs.add(r)
            break
            
    if not is_target:
        continue
        
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        if 'hopper_template_data' in data:
            old_val = data['hopper_template_data']
            if old_val == 'simulation_geometries/2D_hopper.inc':
                data['hopper_template_data'] = 'simulation_geometries/2D_hopper_with_orifice_cover.inc'
                with open(json_file, 'w') as f:
                    json.dump(data, f, indent=4)
                count_updated += 1
                print(f'✅ Updated: {json_file.parent.name}')
            elif old_val == 'simulation_geometries/2D_hopper_with_orifice_cover.inc':
                print(f'⏭️ Already updated: {json_file.parent.name}')
                
    except Exception as e:
        print(f'❌ Error reading {json_file}: {e}')

print(f'\n--- Summary ---')
print(f'Total target runs identified: {len(found_runs)}/{len(runs_to_update)}')
missing = set(runs_to_update) - found_runs
if missing:
    print(f'Missing runs (not found in dumping yard): {missing}')
print(f'Total JSON files successfully updated: {count_updated}')
