# Centralized Lepton bond/angle definitions (lepton.inc)

if "! $(is_defined(variable,r_max))" then "variable r_max equal 0.0005" #pot center to center
# Theta limit in Radians (140 degrees * PI / 180)
if "! $(is_defined(variable,theta_lim))" then "variable theta_lim equal 2.4435"
# Stiffnesses
if "! $(is_defined(variable,knt2))" then "variable knt2   equal 1.0e6"  # High stiffness for bond limit
if "! $(is_defined(variable,ktheta))" then "variable ktheta equal 1.0e1"   # Stiffness for angle limit
# smoothing widths (tune as needed)
if "! $(is_defined(variable,w))" then "variable w   equal 1.0e-4"    # smoothing length (m) 0.1 mm
if "! $(is_defined(variable,delta))" then "variable delta equal 1.0e-4"    # smoothing angle (rad) ~ 1.15 deg

# Granular pair style parameters
# normal contact
if "! $(is_defined(variable,k_n))" then "variable k_n equal 1e6" # Normal stiffness
if "! $(is_defined(variable,nu_n0))" then "variable nu_n0 equal 50.0" # damping prefactor or a coefficient of restitution
# tangential
if "! $(is_defined(variable,k_t))" then "variable k_t equal 1.0e6"   # Tangential stiffness coefficient.
if "! $(is_defined(variable,x_gamma_t))" then "variable x_gamma_t equal 1.0" # dimensionless multiplier for the normal damping
if "! $(is_defined(variable,mu_s))" then "variable mu_s equal 0.3" # is the tangential (or sliding) friction coefficient
# rolling
if "! $(is_defined(variable,k_roll))" then "variable k_roll equal 1.0e6"   # Rolling Friction Stiffness
if "! $(is_defined(variable,gamma_roll))" then "variable gamma_roll equal 5000" # Rolling Friction Damping
if "! $(is_defined(variable,mu_roll))" then "variable mu_roll equal 0.3" # Rolling Friction Coefficient

pair_style  granular
pair_coeff  * * hooke ${k_n} ${nu_n0} & 
            tangential linear_history ${k_t} ${x_gamma_t} ${mu_s} & 
            rolling sds ${k_roll} ${gamma_roll} ${mu_roll} &
            twisting marshall

# BOND: Lepton
bond_style lepton no_offset
bond_coeff 1 0.003 "0.5 * knt2 * (r - r_max)^2 * 0.5 * (1 + tanh((r - r_max)/w)); knt2=v_knt2; r_max=v_r_max; w=v_w"

# ANGLE: Lepton
angle_style lepton no_offset
angle_coeff 1 0 "0.5 * ktheta * (theta - t_lim)^2 * 0.5 * (1 + tanh((t_lim - theta)/w_a)); ktheta=v_ktheta; t_lim=v_theta_lim; w_a=v_delta"