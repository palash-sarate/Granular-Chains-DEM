# Centralized Lepton bond/angle definitions (lepton.inc)

variable r_max equal 0.0005 #pot center to center
# Theta limit in Radians (140 degrees * PI / 180)
variable theta_lim equal 2.4435
# Stiffnesses
variable knt2   equal 1.0e6  # High stiffness for bond limit
variable ktheta equal 1.0e1   # Stiffness for angle limit
# smoothing widths (tune as needed)
variable w   equal 1.0e-4    # smoothing length (m) 0.1 mm
variable delta equal 1.0e-4    # smoothing angle (rad) ~ 1.15 deg

# Granular pair style parameters
# normal contact
# variable k_n equal 1e6 # Normal stiffness
# variable nu_n0 equal 50.0 # damping prefactor or a coefficient of restitution
# rolling
# variable k_roll equal 1.0e6   # Rolling Friction Stiffness
# variable gamma_roll equal 5000 # Rolling Friction Damping
# variable mu_roll equal 0.3 # Rolling Friction Coefficient
# tangential
# variable k_t equal 1.0e6   # Tangential stiffness coefficient.
# variable x_gamma_t equal 1.0 # dimensionless multiplier for the normal damping
# variable mu_s equal 0.3 # is the tangential (or sliding) friction coefficient

pair_style  granular
pair_coeff  * * hooke 1e6 50.0 & 
            tangential linear_history 1.0e6 1.0 0.3 & 
            rolling sds 1.0e6 5000 0.3 &
            twisting marshall

# BOND: Lepton
bond_style lepton no_offset
bond_coeff 1 0.003 "0.5 * knt2 * (r - r_max)^2 * 0.5 * (1 + tanh((r - r_max)/w)); knt2=v_knt2; r_max=v_r_max; w=v_w"

# ANGLE: Lepton
angle_style lepton no_offset
angle_coeff 1 0 "0.5 * ktheta * (theta - t_lim)^2 * 0.5 * (1 + tanh((t_lim - theta)/w_a)); ktheta=v_ktheta; t_lim=v_theta_lim; w_a=v_delta"