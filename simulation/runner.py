import subprocess
import os
import sys
import re
from typing import Optional, List
from .config import SimulationConfig

class SimulationRunner:
    def __init__(self, lammps_executable: str = "lmp"):
        self.lammps_exe = lammps_executable

    def generate_input_script(self, config: SimulationConfig, template_path: str, output_path: str):
        """
        Generates a LAMMPS input script by replacing variables in a template 
        with values from the configuration.
        """
        with open(template_path, 'r') as f:
            content = f.read()
            
        vars_dict = config.to_lammps_vars()
        
        for key, value in vars_dict.items():
            # Regex to find "variable key equal/index/string value"
            # We replace the whole line with the new definition
            # We use 'equal' for numbers to be safe, or string for strings.
            
            # Check if value is a number
            try:
                float(value)
                var_type = "equal"
            except ValueError:
                var_type = "string"
                
            # Special case: outdir is a string
            if key == "outdir":
                var_type = "string"

            pattern = re.compile(r'^\s*variable\s+' + re.escape(key) + r'\s+(equal|index|string)\s+.*$', re.MULTILINE)
            replacement = f"variable {key} {var_type} {value} # Generated from config"
            
            if pattern.search(content):
                content = pattern.sub(replacement, content)
            else:
                # If not found, append it to the beginning of the file (after comments?) or just top
                # But LAMMPS reads sequentially. Variables must be defined before use.
                # If it's not in the template, we assume it's needed.
                # Let's prepend it after the first few lines or just at the top.
                # However, if we prepend, it might be overridden later if we missed the regex.
                # But if regex didn't match, it's not there.
                # Safest is to prepend.
                content = f"variable {key} {var_type} {value} # Injected variable\n" + content
                
        with open(output_path, 'w') as f:
            f.write(content)
        print(f"Generated input script: {output_path}")

    def _prepare_directories(self, config: SimulationConfig):
        """Creates the necessary output directories."""
        outdir = config.output_dir
        # Create main output dir and subdirectories
        subdirs = ["bond", "angle", "restart", "chain"]
        
        # Create base dir
        os.makedirs(outdir, exist_ok=True)
        
        for subdir in subdirs:
            os.makedirs(os.path.join(outdir, subdir), exist_ok=True)
            
        print(f"Prepared output directories in: {outdir}")

    def run(self, config: SimulationConfig, verbose: bool = True):
        """
        Runs the simulation using the provided configuration.
        If template_path is provided, generates a new input script.
        Otherwise, runs config.input_script directly (assuming it's ready).
        """
        # Prepare directories first
        self._prepare_directories(config)
        
        if config.input_script:
            script_to_run = config.input_script        
        elif config.template:
            template_path = f"simulation_templates/{config.template}"
            # Generate a temporary or specific input script
            script_to_run = f"{config.output_dir}/in.{config.simulation}"
            self.generate_input_script(config, template_path, script_to_run)            
        else:
            raise ValueError("Either input_script or template must be provided in the config.")
        
        cmd = [self.lammps_exe, "-in", script_to_run]
        
        # We can still pass variables via command line as a backup or for variables not in the script
        # vars_dict = config.to_lammps_vars()
        # for key, value in vars_dict.items():
        #     cmd.extend(["-var", key, value])
            
        self._execute(cmd, verbose)

    def resume(self, config: SimulationConfig, restart_file: str, resume_template: str = "in.chain_flop_resume", verbose: bool = True):
        """
        Resumes the simulation from a restart file.
        """
        # Generate resume script
        script_to_run = f"generated_{os.path.basename(resume_template)}"
        
        # We need to inject the restart file path into the script
        # The current resume script has 'read_restart post_chain_flop/current/restart.final.bin'
        # We should replace that line.
        
        with open(resume_template, 'r') as f:
            content = f.read()
            
        # Replace read_restart line
        # Look for 'read_restart ...'
        pattern = re.compile(r'^\s*read_restart\s+.*$', re.MULTILINE)
        replacement = f"read_restart {restart_file}"
        content = pattern.sub(replacement, content)
        
        # Also apply config variables
        vars_dict = config.to_lammps_vars()
        for key, value in vars_dict.items():
            try:
                float(value)
                var_type = "equal"
            except ValueError:
                var_type = "string"
            if key == "outdir": var_type = "string"

            var_pattern = re.compile(r'^\s*variable\s+' + re.escape(key) + r'\s+(equal|index|string)\s+.*$', re.MULTILINE)
            var_replacement = f"variable {key} {var_type} {value} # Generated from config"
            if var_pattern.search(content):
                content = var_pattern.sub(var_replacement, content)
                
        with open(script_to_run, 'w') as f:
            f.write(content)
            
        cmd = [self.lammps_exe, "-in", script_to_run]
        self._execute(cmd, verbose)

    def _execute(self, cmd: List[str], verbose: bool):
        print(f"Executing: {' '.join(cmd)}")
        
        if verbose:
            # Stream output to stdout
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if process.stdout:
                for line in process.stdout:
                    print(line, end='')
            process.wait()
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, cmd)
        else:
            subprocess.run(cmd, check=True)
            print("Simulation completed.")
