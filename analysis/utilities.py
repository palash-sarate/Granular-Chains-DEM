import sys
import argparse
import ast
import time
from typing import Optional
import numpy as np
    
def get_dt_token(dt: float) -> str:
    """
    Converts a float dt into a string token suitable for filenames.
    E.g., 1e-6 -> "1e6", 2.5e-7 -> "25e8"
    """
    mantissa, exponent = f"{dt:e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".") or "0"
    exponent = exponent.lstrip("+-")
    exponent = exponent.lstrip("0") or "0"
    dt_token = f"{mantissa}e{exponent}"
    # print(f"DT token: {dt_token}")
    return dt_token

def get_viscosity_token(viscosity: float) -> str:
    """
    Converts a float viscosity into a string token suitable for filenames.
    E.g., 0.003 -> "003", 0.025 -> "025"
    """
    token = int(round(viscosity * 1e9))
    token_str = f"{token:09d}"
    # remove trailing zeros
    token_str = token_str.rstrip("0")
    # print(f"Viscosity token: {token_str}")
    return token_str

def validate_chain_spacing(chain_data, min_spacing=0.0, max_spacing=0.0025):
    """
    Validates spacing between consecutive atoms in the chain.
    Returns (is_valid, messages)
    """
    if not chain_data or len(chain_data) < 2:
        return True, []

    messages = []
    is_valid = True
    
    # Extract just x,y,z
    coords = [np.array(p[:3]) for p in chain_data]
    
    for i in range(len(coords) - 1):
        dist = np.linalg.norm(coords[i+1] - coords[i])
        if dist < min_spacing or dist > max_spacing:
            is_valid = False
            messages.append(f"Gap {i+1}-{i+2}: {dist:.6f} m (Expected {min_spacing}-{max_spacing})")
            
    return is_valid, messages

def get_xyz_series(df, atom_id):
    """
    Returns a list of (timestep, x, y, z) tuples for the whole simulation.
    """
    timesteps = df.index.get_level_values('timestep').unique()
    results = []

    for step in timesteps:
        try:
            pos = df.loc[(step, atom_id), ['x', 'y', 'z']].values
            results.append((step, pos[0], pos[1], pos[2]))
        except KeyError:
            continue
            
    return results

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
class ETAEstimator:
    """Estimate remaining time for a looped job using exponential smoothing.

    Usage:
        eta = ETAEstimator(total=100)
        eta.start()
        for i in range(100):
            # do work
            eta.update(i+1)
            print(eta)
    """

    def __init__(self, total: Optional[int] = None, smoothing: float = 0.2):
        self.total = total
        self.smoothing = float(smoothing)
        self.start_time: Optional[float] = None
        self.last_time: Optional[float] = None
        self.last_count = 0
        self.ema_per_item: Optional[float] = None

    def start(self) -> None:
        self.start_time = time.monotonic()
        self.last_time = self.start_time
        self.last_count = 0
        self.ema_per_item = None

    def update(self, completed: int) -> None:
        """Call after completing `completed` items (1-based count).

        Estimates the per-item time using the delta since last update and
        updates an exponential moving average. Must call `start()` first.
        """
        now = time.monotonic()
        if self.start_time is None:
            self.start()
            now = time.monotonic()

        delta_count = completed - self.last_count
        delta_time = max(1e-9, now - (self.last_time or self.start_time))
        if delta_count > 0:
            per_item = delta_time / delta_count
            if self.ema_per_item is None:
                self.ema_per_item = per_item
            else:
                alpha = self.smoothing
                self.ema_per_item = alpha * per_item + (1 - alpha) * self.ema_per_item

        self.last_time = now
        self.last_count = completed

    def elapsed(self) -> float:
        if self.start_time is None:
            return 0.0
        return time.monotonic() - self.start_time

    def eta_seconds(self) -> Optional[float]:
        if self.ema_per_item is None or self.total is None:
            return None
        remaining = max(0, int(self.total) - int(self.last_count))
        return remaining * self.ema_per_item

    def progress_fraction(self) -> Optional[float]:
        if self.total is None:
            return None
        return min(1.0, float(self.last_count) / float(self.total))

    def format_seconds(self, s: Optional[float]) -> str:
        if s is None:
            return "--:--:--"
        s = int(round(s))
        h = s // 3600
        m = (s % 3600) // 60
        sec = s % 60
        return f"{h:02d}:{m:02d}:{sec:02d}"

    def __str__(self) -> str:
        elapsed = self.elapsed()
        eta = self.eta_seconds()
        frac = self.progress_fraction()
        parts = [f"elapsed={self.format_seconds(elapsed)}"]
        if frac is not None:
            parts.append(f"{int(100*frac):3d}%")
        parts.append(f"ETA={self.format_seconds(eta)}")
        return " ".join(parts)

def main():
    available_functions = {
        "main": main,
        "get_viscosity_token": get_viscosity_token,
        "get_dt_token": get_dt_token,
    }

    parser = argparse.ArgumentParser(description="Execute functions from main.py")
    parser.add_argument("func", nargs="?", default="main", choices=available_functions.keys())
    parser.add_argument("func_args", nargs="*", help="Positional or key=value arguments")
    parsed = parser.parse_args()

    def _coerce(token: str):
        try:
            return ast.literal_eval(token)
        except (ValueError, SyntaxError):
            return token

    positional, keyword = [], {}
    for token in parsed.func_args:
        if "=" in token:
            key, value = token.split("=", 1)
            keyword[key] = _coerce(value)
        else:
            positional.append(_coerce(token))

    try:
        available_functions[parsed.func](*positional, **keyword)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(1)

if __name__ == "__main__":
    main()