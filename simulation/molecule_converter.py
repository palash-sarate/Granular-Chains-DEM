import math
import sys
import os

def parse_lammps_data(data_file):
    """
    Parses a LAMMPS data file (atom_style hybrid sphere angle/molecular).
    Returns a dictionary with atoms, bonds, angles, and box info.
    """
    with open(data_file, 'r') as f:
        lines = f.readlines()

    data = {
        'atoms': [],
        'bonds': [],
        'angles': [],
        'masses': {}
    }

    section = None
    
    # Mapping for atom_style sphere: id type diameter density x y z
    # But write_data might output differently depending on flags.
    # We will assume standard "id type diameter density x y z" for sphere
    # or "id type x y z ..." if it detects molecular.
    # Let's try to be robust.
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        if "Atoms" in line:
            section = "atoms"
            continue
        elif "Bonds" in line:
            section = "bonds"
            continue
        elif "Angles" in line:
            section = "angles"
            continue
        elif "Velocities" in line:
            section = "velocities"
            continue
        elif "Masses" in line:
            section = "masses"
            continue
            
        # Parse content based on section
        parts = line.split()
        
        if section == "atoms":
            # Observed format for hybrid sphere angle: id type x y z diameter density ...
            if len(parts) >= 7:
                try:
                    atom_id = int(parts[0])
                    atom_type = int(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    diameter = float(parts[5])
                    density = float(parts[6])
                    
                    data['atoms'].append({
                        'id': atom_id,
                        'type': atom_type,
                        'diameter': diameter,
                        'density': density,
                        'x': x, 'y': y, 'z': z
                    })
                except ValueError:
                    pass # Header line or malformed
                    
        elif section == "bonds":
            if len(parts) >= 4:
                # id type atom1 atom2
                data['bonds'].append({
                    'id': int(parts[0]),
                    'type': int(parts[1]),
                    'atom1': int(parts[2]),
                    'atom2': int(parts[3])
                })
                
        elif section == "angles":
            if len(parts) >= 5:
                # id type atom1 atom2 atom3
                data['angles'].append({
                    'id': int(parts[0]),
                    'type': int(parts[1]),
                    'atom1': int(parts[2]),
                    'atom2': int(parts[3]),
                    'atom3': int(parts[4])
                })

    # Sort atoms by ID to ensure order
    data['atoms'].sort(key=lambda x: x['id'])
    return data

def write_molecule_file(data, output_file):
    """
    Writes a LAMMPS molecule file from the parsed data.
    Centers the molecule at (0,0,0).
    """
    atoms = data['atoms']
    if not atoms:
        print("No atoms found to write.")
        return

    # Calculate Center of Mass (or geometric center)
    # Using geometric center for simplicity
    xs = [a['x'] for a in atoms]
    ys = [a['y'] for a in atoms]
    zs = [a['z'] for a in atoms]
    
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    
    with open(output_file, 'w') as f:
        f.write(f"# Molecule file generated from relaxed chain\n\n")
        
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(data['bonds'])} bonds\n")
        f.write(f"{len(data['angles'])} angles\n\n")
        
        f.write("Coords\n\n")
        for i, a in enumerate(atoms):
            # Molecule IDs usually start at 1
            f.write(f"{i+1} {a['x'] - cx:.6f} {a['y'] - cy:.6f} {a['z'] - cz:.6f}\n")
        f.write("\n")
        
        f.write("Types\n\n")
        for i, a in enumerate(atoms):
            f.write(f"{i+1} {a['type']}\n")
        f.write("\n")
        
        f.write("Diameters\n\n")
        for i, a in enumerate(atoms):
            f.write(f"{i+1} {a['diameter']:.6f}\n")
        f.write("\n")
        
        f.write("Masses\n\n")
        for i, a in enumerate(atoms):
            # Calculate mass from density and diameter
            radius = a['diameter'] / 2.0
            volume = (4.0/3.0) * math.pi * (radius**3)
            mass = a['density'] * volume
            f.write(f"{i+1} {mass:.6e}\n")
        f.write("\n")
        
        if data['bonds']:
            f.write("Bonds\n\n")
            for i, b in enumerate(data['bonds']):
                # Remap atom IDs to 1..N
                # Assuming input IDs are contiguous or we map them
                # Here we assume the list is sorted and corresponds to 1..N
                # But we should be careful if IDs are not 1-based in input
                # For a single chain data file, they usually are.
                f.write(f"{i+1} {b['type']} {b['atom1']} {b['atom2']}\n")
            f.write("\n")
            
        if data['angles']:
            f.write("Angles\n\n")
            for i, ang in enumerate(data['angles']):
                f.write(f"{i+1} {ang['type']} {ang['atom1']} {ang['atom2']} {ang['atom3']}\n")
            f.write("\n")

    print(f"Converted to molecule: {output_file}")

def convert_data_to_molecule(data_file, output_mol_file):
    data = parse_lammps_data(data_file)
    write_molecule_file(data, output_mol_file)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        convert_data_to_molecule(sys.argv[1], sys.argv[2])
    else:
        print("Usage: python molecule_converter.py <input_data_file> <output_mol_file>")
