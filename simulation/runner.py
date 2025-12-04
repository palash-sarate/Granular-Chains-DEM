import subprocess
import os
import sys
import re
from typing import Optional, List
from .config import SimulationConfig
# from lammps import lammps

import shutil
import json
import datetime

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
        
        # Remove existing outdir variable definition if present to avoid conflicts
        pattern = re.compile(r'^\s*variable\s+outdir\s+string\s+.+$', re.MULTILINE)
        content = pattern.sub("", content)
        
        # Inject ALL extra_vars as LAMMPS variables at the top of the script
        header_vars = []
        
        # Always inject outdir
        outdir_normalized = config.output_dir.replace("\\", "/")
        header_vars.append(f"variable outdir string {outdir_normalized}")

        for key, value in config.extra_vars.items():
            header_vars.append(f"variable {key} string {value}")
        
        # Also inject data_file as a variable if present
        if config.data_file:
            # Check if it's already a full path or needs prefix
            if config.data_file.startswith("chain_data/") or config.data_file.startswith("chain_data\\"):
                 data_path = config.data_file
            elif "/" in config.data_file or "\\" in config.data_file:
                 # Assume it's a relative path provided by user (like temp_lib_gen/...)
                 data_path = f"chain_data/{config.data_file}"
            else:
                 # Default location
                 data_path = f"chain_data/{config.data_file}"
            
            # Normalize slashes for LAMMPS
            data_path = data_path.replace("\\", "/")
            header_vars.append(f"variable data_file string {data_path}")

        if header_vars:
            content = "\n".join(header_vars) + "\n\n" + content

        with open(output_path, 'w') as f:
            f.write(content)
        print(f"Generated input script: {output_path}")

    def _prepare_directories(self, config: SimulationConfig, clean: bool = False):
        """Creates the necessary output directories."""
        outdir = config.output_dir
        
        # Clean directory if requested
        if clean and os.path.exists(outdir):
            print(f"Cleaning output directory: {outdir}")
            shutil.rmtree(outdir)
            
        # Create main output dir and subdirectories
        subdirs = ["bond", "angle", "restart", "chain"]
        
        # Create base dir
        os.makedirs(outdir, exist_ok=True)
        
        for subdir in subdirs:
            os.makedirs(os.path.join(outdir, subdir), exist_ok=True)
            
        print(f"Prepared output directories in: {outdir}")

    def run(self, config: SimulationConfig, verbose: bool = True, clean_dir: bool = True):
        """
        Runs the simulation using the provided configuration.
        If template_path is provided, generates a new input script.
        Otherwise, runs config.input_script directly (assuming it's ready).
        """
        # Prepare directories first
        self._prepare_directories(config, clean=clean_dir)
        
        if config.input_script:
            script_to_run = config.input_script
            # save the script to output dir for record-keeping
            dest_script = os.path.join(config.output_dir, os.path.basename(script_to_run))
            subprocess.run(["cp", script_to_run, dest_script], check=True)
            script_to_run = dest_script

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
            # Run directly to allow real-time output to terminal
            subprocess.run(cmd, check=True)
        else:
            # Suppress output
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print("Simulation completed.")
