from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List
import json

@dataclass
class SimulationConfig:    
    # Input/Output
    input_script: Optional[str] = None
    template: str = "in.default_chain_template"
    resume_file: Optional[str] = None  # New field for resume file path
    data_file: str = "default_chain.data"
    run: str = "default_run"
    simulation: str = "default_simulation"
    
    # Additional LAMMPS variables
    extra_vars: Dict[str, Any] = field(default_factory=dict)

    # New optional fields for raw LAMMPS lines from JSON
    regions: Optional[List[str]] = None        # each entry: a 'region ...' line or multiline chunk
    wall_blocks: Optional[List[str]] = None    # each entry: a 'fix ...' wall block (string / multiline)
    walls: Optional[Dict[str, Any]] = None     # keep-compatible: structured walls (converted by runner if present)

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

    @classmethod
    def from_json(cls, path: str) -> "SimulationConfig":
        with open(path, "r") as f:
            data = json.load(f)
        return cls(**data)