import numpy as np

def calculate_angle_3d(pos_a, pos_b, pos_c):
    """
    Calculates angle ABC (at vertex B) in degrees.
    Inputs: numpy arrays or lists [x, y, z]
    """
    a = np.array(pos_a)
    b = np.array(pos_b)
    c = np.array(pos_c)

    ba = a - b
    bc = c - b

    norm_ba = np.linalg.norm(ba)
    norm_bc = np.linalg.norm(bc)

    if norm_ba == 0 or norm_bc == 0:
        return 0.0

    cosine_angle = np.dot(ba, bc) / (norm_ba * norm_bc)
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    
    return np.degrees(np.arccos(cosine_angle))

def get_angle_series(df, id1, id2, id3):
    """
    Returns a list of (timestep, angle) tuples for the whole simulation.
    """
    timesteps = df.index.get_level_values('timestep').unique()
    results = []

    for step in timesteps:
        try:
            # .loc lookup is very fast on MultiIndex
            p1 = df.loc[(step, id1), ['x', 'y', 'z']].values
            p2 = df.loc[(step, id2), ['x', 'y', 'z']].values
            p3 = df.loc[(step, id3), ['x', 'y', 'z']].values
            
            angle = calculate_angle_3d(p1, p2, p3)
            results.append((step, angle))
        except KeyError:
            # Handle cases where an atom might be missing (e.g., lost atoms)
            continue
            
    return results

def get_distance_series(df, id1, id2):
    """
    Returns a list of (timestep, distance) tuples for the whole simulation.
    """
    timesteps = df.index.get_level_values('timestep').unique()
    results = []

    for step in timesteps:
        try:
            p1 = df.loc[(step, id1), ['x', 'y', 'z']].values
            p2 = df.loc[(step, id2), ['x', 'y', 'z']].values
            
            distance = np.linalg.norm(p1 - p2)
            results.append((step, distance))
        except KeyError:
            continue
            
    return results