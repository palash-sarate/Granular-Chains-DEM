import os
from pathlib import Path

def load_env(env_path: str = None) -> dict:
    """
    Loads variables from a .env file into os.environ.
    If env_path is not specified, it searches for a .env file
    in the project root directory relative to this file.
    Returns a dictionary of loaded environment variables.
    """
    if env_path is None:
        # Search upward from this file's directory to find the project root containing .env
        current_dir = Path(__file__).resolve().parent
        while current_dir != current_dir.parent:
            candidate = current_dir / ".env"
            if candidate.exists():
                env_path = str(candidate)
                break
            current_dir = current_dir.parent
        
        # Fallback to local .env if not found
        if not env_path:
            env_path = ".env"

    env_vars = {}
    if os.path.exists(env_path):
        with open(env_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#"):
                    continue
                
                # Split at first '=' sign
                if "=" in line:
                    key, val = line.split("=", 1)
                    key = key.strip()
                    val = val.strip()
                    
                    # Remove inline comments if present
                    if " #" in val or "\t#" in val:
                        val = val.split(" #", 1)[0].strip()
                    elif val.startswith("#"):
                        continue
                    else:
                        # Handle trailing comment with #
                        parts = val.split("#", 1)
                        if len(parts) > 1:
                            val = parts[0].strip()
                    
                    # Strip wrapping quotes
                    if val.startswith('"') and val.endswith('"'):
                        val = val[1:-1]
                    elif val.startswith("'") and val.endswith("'"):
                        val = val[1:-1]
                        
                    os.environ[key] = val
                    env_vars[key] = val
                    
    # Also set default fallback variables if they don't exist
    defaults = {
        "LAMMPS_EXECUTABLE": "lmp",
        "DUMPING_YARD": "dumping_yard",
        "QSTAT_PATH": "/opt/pbs/bin/qstat",
        "QSUB_PATH": "/opt/pbs/bin/qsub"
    }
    for k, v in defaults.items():
        if k not in os.environ:
            os.environ[k] = v
            env_vars[k] = v
            
    return env_vars
