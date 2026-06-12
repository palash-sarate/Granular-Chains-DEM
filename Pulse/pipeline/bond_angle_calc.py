import os
import pandas as pd
from analysis.data_manager import calculate_bonds_angles_vectorized

def run_bond_angle_calc(run_path, atoms_path, bonds_path, angles_path):
    print(f"Calculating bond lengths and angles from: {atoms_path}")
    if not os.path.exists(atoms_path):
        raise FileNotFoundError(f"Input atomic coordinate cache not found: {atoms_path}")
        
    df_atoms = pd.read_parquet(atoms_path)
    df_bonds, df_angles = calculate_bonds_angles_vectorized(df_atoms)
    
    if df_bonds.empty or df_angles.empty:
        raise ValueError("Calculated bond or angle datasets are empty. Verify raw coordinate structure.")
        
    os.makedirs(os.path.dirname(bonds_path), exist_ok=True)
    os.makedirs(os.path.dirname(angles_path), exist_ok=True)
    
    df_bonds.to_parquet(bonds_path)
    df_angles.to_parquet(angles_path)
    print(f"Successfully saved bonds to: {bonds_path} and angles to: {angles_path}")
