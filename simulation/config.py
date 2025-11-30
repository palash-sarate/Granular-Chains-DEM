from dataclasses import dataclass, field
from typing import Dict, Any, Optional

@dataclass
class SimulationConfig:    
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
