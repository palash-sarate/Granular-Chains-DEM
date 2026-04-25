import os

files = [
    "chain_0.dump",
    "chain_100000.dump",
    "chain_200000.dump",
    "chain_300000.dump",
    "chain_357000.dump"
]
base_path = r"d:\Chains_simulations\temp"
particle_id = "5500"
dt = 1e-6

results = []

for filename in files:
    filepath = os.path.join(base_path, filename)
    step = int(filename.split("_")[1].split(".")[0])
    time = step * dt
    
    with open(filepath, 'r') as f:
        atom_data_started = False
        for line in f:
            if line.startswith("ITEM: ATOMS"):
                atom_data_started = True
                continue
            if atom_data_started:
                parts = line.split()
                if parts and parts[0] == particle_id:
                    # id mol type x y z vx vy vz ...
                    x, y, z = map(float, parts[3:6])
                    results.append({
                        "step": step,
                        "time": time,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                    break

print("Step | Time | X | Y | Z")
for r in results:
    print(f"{r['step']} | {r['time']:.6f} | {r['x']} | {r['y']} | {r['z']}")

print("\nVelocities between consecutive pairs:")
for i in range(len(results) - 1):
    r1 = results[i]
    r2 = results[i+1]
    dt_pair = r2['time'] - r1['time']
    dx = r2['x'] - r1['x']
    dy = r2['y'] - r1['y']
    dz = r2['z'] - r1['z']
    v_mag = (dx**2 + dy**2 + dz**2)**0.5 / dt_pair
    vx = dx / dt_pair
    vy = dy / dt_pair
    vz = dz / dt_pair
    print(f"Steps {r1['step']} to {r2['step']}:")
    print(f"  dt: {dt_pair:.6f} s")
    print(f"  Vx: {vx:.4f}, Vy: {vy:.4f}, Vz: {vz:.4f}")
    print(f"  V_mag: {v_mag:.4f}")
