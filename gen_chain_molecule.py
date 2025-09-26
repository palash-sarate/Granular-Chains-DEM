#!/usr/bin/env python3
"""
Generate a LAMMPS molecule file for a flexible chain of spheres.
- Chain length N is adjustable
- Bead spacing: 1 mm (0.001 m)
- Bonds: between consecutive beads
- Angles: between consecutive triplets, can be used for angle constraints

Usage:
    python gen_chain_molecule.py N > chain.mol
"""
import sys
import math

# Parameters
radius = 0.0005  # 0.5 mm bead radius
spacing = 0.001  # 1 mm bond length
atom_type = 1
bond_type = 1
angle_type = 1

if len(sys.argv) < 2:
    print("Usage: python gen_chain_molecule.py N > chain.mol")
    sys.exit(1)

N = int(sys.argv[1])

print(f"LAMMPS molecule file for chain of {N} beads\n")
print(f"{N} atoms")
print(f"{N-1} bonds")
print(f"{N-2} angles\n")

print("Coords")
for i in range(N):
    x = i * spacing
    print(f"{i+1} {atom_type} {x:.6f} 0.0 0.0")

print("\nBonds")
for i in range(N-1):
    print(f"{i+1} {bond_type} {i+1} {i+2}")

print("\nAngles")
for i in range(N-2):
    print(f"{i+1} {angle_type} {i+1} {i+2} {i+3}")
