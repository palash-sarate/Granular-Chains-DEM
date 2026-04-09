import os
import sys
import argparse
import ast
import time
from pathlib import Path

from simulation.orchestrator import SimulationOrchestrator

def main():
    """Main entry point for CLI. Instantiate orchestrator and handle help/defaults."""
    print("LAMMPS Chain Simulation CLI")
    print("Use --help to see available functions and arguments.")

if __name__ == "__main__":
    orchestrator = SimulationOrchestrator(lammps_executable="lmp")
    
    # Map CLI command names to Orchestrator methods
    available_functions = {
        "main": main,
        "run_flop_simulations": orchestrator.run_flop_batch,
        "run_flop_simulation": orchestrator.run_flop_simulation,
        "resume_flop_simulation": orchestrator.resume_flop_simulation,
        "generate_linear_chains": orchestrator.generate_chains,
        "generate_relaxed_chain_states": orchestrator.generate_relaxed_library,
        "run_hopper_fill": orchestrator.run_hopper_fill,
        "resume_hopper_fill": orchestrator.resume_hopper_fill,
    }

    parser = argparse.ArgumentParser(description="Execute functions from the Simulation Orchestrator")
    parser.add_argument("func", nargs="?", default="main", choices=available_functions.keys())
    parser.add_argument("func_args", nargs="*", help="Positional or key=value arguments")
    parsed = parser.parse_args()

    def _coerce(token: str):
        try:
            return ast.literal_eval(token)
        except (ValueError, SyntaxError):
            return token

    positional, keyword = [], {}
    i = 0
    while i < len(parsed.func_args):
        token = parsed.func_args[i]
        
        if "=" in token:
            key, value = token.split("=", 1)
            if not value and i + 1 < len(parsed.func_args):
                value = parsed.func_args[i+1]
                i += 1
            keyword[key.strip()] = _coerce(value.strip())
        elif i + 1 < len(parsed.func_args) and parsed.func_args[i+1] == "=":
            key = token
            if i + 2 < len(parsed.func_args):
                value = parsed.func_args[i+2]
                i += 2
            else:
                value = ""
                i += 1
            keyword[key.strip()] = _coerce(value.strip())
        else:
            positional.append(_coerce(token))
        i += 1

    try:
        available_functions[parsed.func](*positional, **keyword)
    except KeyboardInterrupt:
        print("\nProgram interrupted by user. Exiting gracefully.")
        sys.exit(1)
    except TypeError as e:
        print(f"Error calling '{parsed.func}': {e}")
        print(f"Arguments provided: positional={positional}, keyword={keyword}")
        sys.exit(1)
    except Exception as e:
        print(f"Simulation Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)