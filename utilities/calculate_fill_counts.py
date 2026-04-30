import math
import argparse

def calculate_chains(orifice_w, hopper_w, angle_deg, fill_height_pct, total_height, bead_diameter=2.0):
    """
    Calculates number of chains required to fill a 2D hopper.
    
    Args:
        orifice_w: Width of the orifice (cm)
        hopper_w: Total width of the hopper (cm)
        angle_deg: Angle of converging section to horizontal (degrees)
        fill_height_pct: Percentage of total height to fill (%)
        total_height: Total height of the hopper from orifice to top (cm)
        bead_diameter: Diameter of a single bead (mm)
    """
    
    # Constants: Area fractions extracted from provided data
    phi_values = {
        4: 0.755,
        12: 0.58,
        24: 0.56,
        48: 0.50
    }
    
    # Convert units to cm
    bead_radius_cm = (bead_diameter / 10.0) / 2.0
    bead_area = math.pi * (bead_radius_cm**2)
    
    # Calculate height of converging section
    # H = tan(angle) * (W_hopper - W_orifice) / 2
    angle_rad = math.radians(angle_deg)
    h_converging = math.tan(angle_rad) * (hopper_w - orifice_w) / 2.0
    
    target_height = (fill_height_pct / 100.0) * total_height
    
    # Calculate target area
    if target_height <= h_converging:
        # Fill is within the converging trapezoid
        width_at_height = orifice_w + (2.0 * target_height / math.tan(angle_rad))
        area = (orifice_w + width_at_height) / 2.0 * target_height
    else:
        # Fill covers entire converging section and part of straight section
        area_converging = (orifice_w + hopper_w) / 2.0 * h_converging
        area_straight = hopper_w * (target_height - h_converging)
        area = area_converging + area_straight
        
    print(f"\n--- Hopper Geometry ---")
    print(f"Hopper Width: {hopper_w} cm")
    print(f"Orifice Width: {orifice_w} cm")
    print(f"Angle: {angle_deg}°")
    print(f"Total Height: {total_height} cm")
    print(f"Converging Section Height: {h_converging:.2f} cm")
    print(f"Target Fill Height: {target_height:.2f} cm ({fill_height_pct}%)")
    print(f"Target Fill Area: {area:.2f} cm^2")
    
    print(f"\n--- Chain Counts (per hopper) ---")
    print(f"{'N':<5} | {'Phi':<10} | {'Chains Required':<15} | {'Total Beads':<12}")
    print("-" * 50)
    
    results = {}
    for N, phi in phi_values.items():
        total_beads = (area * phi) / bead_area
        chains = total_beads / N
        results[N] = round(chains)
        print(f"{N:<5} | {phi:<10.3f} | {round(chains):<15} | {round(total_beads):<12}")
        
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate chain counts for hopper filling.")
    parser.add_argument("--orifice", type=float, default=10.0, help="Orifice width (cm)")
    parser.add_argument("--width", type=float, default=31.0, help="Hopper width (cm)")
    parser.add_argument("--angle", type=float, default=60.0, help="Hopper angle (deg)")
    parser.add_argument("--height", type=float, default=40.0, help="Total hopper height (cm)")
    parser.add_argument("--fill", type=float, default=80.0, help="Fill height percentage (0-100)")
    parser.add_argument("--bead-d", type=float, default=2.0, help="Bead diameter (mm)")
    parser.add_argument("--n-hoppers", type=int, default=1, help="Multiply counts for N hoppers")

    args = parser.parse_args()
    
    counts = calculate_chains(args.orifice, args.width, args.angle, args.fill, args.height, args.bead_d)
    
    if args.n_hoppers > 1:
        print(f"\n--- Total Counts for {args.n_hoppers} Hoppers ---")
        total_str = ",".join([str(counts[N] * args.n_hoppers) for N in sorted(counts.keys())])
        print(f"--n_fill {total_str}")
        print(f"--N {','.join([str(N) for N in sorted(counts.keys())])}")
