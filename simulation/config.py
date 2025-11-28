from dataclasses import dataclass, field
from typing import Dict, Any, Optional

@dataclass
class SimulationConfig:
    # Simulation parameters
    timestep: float = 1.0e-6
    run_steps: int = 1000000
    thermo_interval: int = 10000
    dump_interval: int = 1000
    restart_interval: int = 100000
    
    # Physical parameters
    r_max: float = 0.0005
    theta_lim: float = 2.4435
    
    # Stiffness and pot. smoothing
    knt2: float = 1.0e6
    ktheta: float = 1.0e6
    w: float = 1.0e-4
    delta: float = 1.0e-4
    viscosity: float = 0.027
    
    # Granular pair style parameters
    k_n: float = 1e6
    nu_n0: float = 50.0
    k_roll: float = 1.0e6
    gamma_roll: float = 5000.0
    mu_roll: float = 0.3
    k_t: float = 1.0e6
    x_gamma_t: float = 1.0
    mu_s: float = 0.3
    
    # Input/Output
    input_script: Optional[str] = None
    template: str = "in.default_chain_template"
    data_file: str = "default_chain.data"
    run: str = "default_run"
    simulation: str = "default_simulation"
    
    # Additional LAMMPS variables
    extra_vars: Dict[str, Any] = field(default_factory=dict)

    @property
    def output_dir(self) -> str:
        """Constructs the output directory path from simulation and run."""
        # Use forward slashes for LAMMPS compatibility even on Windows
        return f"dumping_yard/{self.simulation}/{self.run}"

    def to_lammps_vars(self) -> Dict[str, str]:
        """Converts config to a dictionary of LAMMPS variable arguments."""
        vars_dict = {
            "knt2": str(self.knt2),
            "ktheta": str(self.ktheta),
            "w": str(self.w),
            "delta": str(self.delta),
            "k_n": str(self.k_n),
            "nu_n0": str(self.nu_n0),
            "k_roll": str(self.k_roll),
            "gamma_roll": str(self.gamma_roll),
            "mu_roll": str(self.mu_roll),
            "k_t": str(self.k_t),
            "x_gamma_t": str(self.x_gamma_t),
            "mu_s": str(self.mu_s),
            "r_max": str(self.r_max),
            "theta_lim": str(self.theta_lim),
            "outdir": str(self.output_dir),
            "template": str(self.template),
            "data_file":str(self.data_file),
            "run":str(self.run),
            "simulation":str(self.simulation)
        }
        # Add any extra variables
        for k, v in self.extra_vars.items():
            vars_dict[k] = str(v)
            
        return vars_dict
