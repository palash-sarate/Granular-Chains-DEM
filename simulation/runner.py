import subprocess
import os
import sys
import re
from typing import List
from .config import SimulationConfig

import shutil

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

        if config.lepton_file:
            lepton_path = config.lepton_file.replace("\\", "/")
            header_vars.append(f"variable lepton_inc string {lepton_path}")

        # Resume functionality: inject resume specific variables
        # resume_file
        if config.resume_file:
            resume_path = config.resume_file.replace("\\", "/")
            header_vars.append(f"variable resume_file string {resume_path}")
            
        if header_vars:
            content = "\n".join(header_vars) + "\n\n" + content

        # --- INSERT REGIONS/WALLS FROM CONFIG ---
        # placeholder to look for in templates
        placeholder = "# <REGIONS_WALLS>"

        blocks: List[str] = []
        # raw region lines (user-provided)
        if config.regions:
            blocks.extend(config.regions)
        # raw wall/fix blocks (user-provided)
        if config.wall_blocks:
            blocks.extend(config.wall_blocks)

        # also support structured `walls` dict (backwards compatible)
        if config.walls:
            # format structured walls into fix lines (x/y/z if present)
            def _format_wall(axis: str, walls: dict) -> str:
                mat = walls.get("material", {})
                tang = walls.get("tangential", {})
                roll = walls.get("rolling", {})
                twist = walls.get("twisting", "marshall")
                bounds = walls.get(axis, {})
                lo = bounds.get("min", 0.0)
                hi = bounds.get("max", 0.0)
                return (
                    f"fix {axis}walls all wall/gran granular "
                    f"hertz/material {mat.get('kn',1.0e8)} {mat.get('nu',0.3)} {mat.get('rest',5e-4)} &\n"
                    f"    tangential {tang.get('style','linear_history')} {tang.get('kt',300.0)} {tang.get('gamma',1.0)} {tang.get('mu',0.1)} &\n"
                    f"    rolling {roll.get('style','sds')} {roll.get('kr',200.0)} {roll.get('gamma',100.0)} {roll.get('mu',0.1)} &\n"
                    f"    twisting {twist} &\n"
                    f"    {axis}plane {lo} {hi}\n\n"
                )
            for ax in ("x", "y", "z"):
                if ax in config.walls:
                    blocks.append(_format_wall(ax, config.walls))

        if blocks:
            block_text = "\n".join(blocks) + "\n\n"
            if placeholder in content:
                content = content.replace(placeholder, block_text)
            else:
                # place under the WALL DEFINITIONS section if it exists
                marker = "# --- WALL DEFINITIONS ---"
                if marker in content:
                    content = content.replace(marker, marker + "\n" + block_text)
                else:
                    # fallback: append at end
                    content += "\n# --- GENERATED REGIONS/WALLS ---\n" + block_text
        else:
            # remove placeholder if present
            content = content.replace(placeholder, "")

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

    def run(self, config: SimulationConfig, verbose: bool = True, clean_dir: bool = True, prep_dirs: bool = True):
        """
        Runs the simulation using the provided configuration.
        If template_path is provided, generates a new input script.
        Otherwise, runs config.input_script directly (assuming it's ready).
        """
        # Prepare directories first
        if prep_dirs:
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
  
        cmd = [self.lammps_exe, "-in", script_to_run, "-log", f"{config.output_dir}/lammps.log"]
        
        # We can still pass variables via command line as a backup or for variables not in the script
        # vars_dict = config.to_lammps_vars()
        # for key, value in vars_dict.items():
        #     cmd.extend(["-var", key, value])
            
        self._execute(cmd, verbose)

    def resume(self, config: SimulationConfig, verbose: bool = True):
        """
        Resumes the simulation from a restart file.
        """

        if config.template:
            template_path = f"simulation_templates/{config.template}"
            # Generate a temporary or specific input script
            script_to_run = f"{config.output_dir}/in.{config.simulation}.resume"
            self.generate_input_script(config, template_path, script_to_run)            
        else:
            raise ValueError("Resume template must be provided in the config.")
        
        cmd = [self.lammps_exe, "-in", script_to_run, "-log", f"{config.output_dir}/lammps_resume.log"]
        
        # We can still pass variables via command line as a backup or for variables not in the script
        # vars_dict = config.to_lammps_vars()
        # for key, value in vars_dict.items():
        #     cmd.extend(["-var", key, value])
            
        self._execute(cmd, verbose)

    def _execute(self, cmd: List[str], verbose: bool):
        # If a log_file is provided, write both stdout and stderr to it.
        # When verbose is True and no log_file is provided, stream to terminal.
        cmd_str = ' '.join(cmd)

        if verbose:
            print(f"Executing: {cmd_str}")
            subprocess.run(cmd, check=True)
        else:
            # Discard output
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print("Simulation completed.")
