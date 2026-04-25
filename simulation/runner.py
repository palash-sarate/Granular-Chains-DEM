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
        self.output_dir = ""
        self._features = None # Cache for lmp -h features

    def _get_lammps_features(self) -> List[str]:
        """Detect supported packages from lmp -h."""
        if self._features is not None:
            return self._features
        try:
            result = subprocess.run([self.lammps_exe, "-h"], 
                                 capture_output=True, text=True, timeout=5)
            # Find the 'Installed packages:' section
            packages = []
            capture = False
            for line in result.stdout.splitlines():
                if "Installed packages:" in line:
                    capture = True
                    continue
                if capture:
                    if line.startswith("List of"):
                        break
                    if not line.strip():
                        continue
                    packages.extend(line.split())
            self._features = packages
            return packages
        except Exception:
            self._features = []
            return []

    def _has_gpu(self) -> bool:
        """Detect if local machine has an NVIDIA GPU via nvidia-smi."""
        return shutil.which("nvidia-smi") is not None

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

        if config.dump_file:
            dump_path = config.dump_file.replace("\\", "/")
            header_vars.append(f"variable dump_inc string {dump_path}")

        # Resume functionality: inject resume specific variables
        if config.resume_file:
            resume_path = config.resume_file.replace("\\", "/")
            header_vars.append(f"variable resume_file string {resume_path}")
            
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
            print(f"Surgical cleanup of output directory: {outdir}")
            # Instead of deleting everything, we delete specific subdirs and files
            # but preserve sim_metadata.json
            for item in os.listdir(outdir):
                item_path = os.path.join(outdir, item)
                if item == "sim_metadata.json":
                    continue
                try:
                    if os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                    else:
                        os.remove(item_path)
                except Exception as e:
                    print(f"Warning: Could not remove {item_path}: {e}")
            
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
        self.output_dir = config.output_dir
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
        
        self._apply_acceleration(cmd, config)
        self._execute(cmd, config, verbose)

    def resume(self, config: SimulationConfig, verbose: bool = True, prep_dirs: bool = True):
        """
        Resumes the simulation from a restart file.
        """
        self.output_dir = config.output_dir
        if prep_dirs:
            # For resumption, we typically don't want to clean the directory
            self._prepare_directories(config, clean=False)

        if config.template:
            template_path = f"simulation_templates/{config.template}"
            # Generate a temporary or specific input script
            script_to_run = f"{config.output_dir}/in.{config.simulation}.resume"
            self.generate_input_script(config, template_path, script_to_run)            
        else:
            raise ValueError("Resume template must be provided in the config.")
        
        cmd = [self.lammps_exe, "-in", script_to_run, "-log", f"{config.output_dir}/lammps_resume.log"]
        
        self._apply_acceleration(cmd, config)
        self._execute(cmd, config, verbose)

    def _apply_acceleration(self, cmd: List[str], config: SimulationConfig):
        """Detects hardware and applies optimal acceleration flags to the command."""
        supported_pkgs = self._get_lammps_features()
        gpu_available = self._has_gpu()

        # Try KOKKOS
        if config.use_kokkos and "KOKKOS" in supported_pkgs:
            kokkos_cmd = ["-k", "on"]
            if gpu_available:
                kokkos_cmd.extend(["g", "1"])
                print(f"Applying KOKKOS acceleration (GPU=Detected, threads={config.num_threads})")
            else:
                print(f"WARNING: GPU not detected via nvidia-smi. Falling back to KOKKOS-CPU (threads={config.num_threads})")
            
            kokkos_cmd.extend([
                "t", str(config.num_threads),
                "-sf", "kk",
                "-pk", "kokkos", "newton", "on", "neigh", "half"
            ])
            cmd.extend(kokkos_cmd)
            return

        # Try INTEL
        if config.use_intel and "INTEL" in supported_pkgs:
            cmd.extend(["-sf", "intel"])
            print(f"Applying INTEL acceleration (threads={config.num_threads})")
            return
        
        # Fallback to OPENMP
        if config.num_threads > 1 and "OPENMP" in supported_pkgs:
            cmd.extend(["-sf", "omp", "-pk", "omp", str(config.num_threads)])
            print(f"Applying OPENMP acceleration (threads={config.num_threads})")
        else:
            print(f"Standard execution (no acceleration) with {config.num_threads} threads.")

    def _execute(self, cmd: List[str], config: SimulationConfig, verbose: bool):
        # 1. Wrap command with mpiexec if parallelism is requested
        nprocs = config.num_procs if config.num_procs is not None else 1
        if nprocs and nprocs > 1:
            if config.num_threads > 1:
                # Hybrid MPI+OpenMP: Bind each MPI process to a set of cores equal to num_threads
                cmd = ["mpiexec", "--map-by", f"socket:PE={config.num_threads}", "--bind-to", "core", "-n", str(nprocs)] + cmd
            else:
                cmd = ["mpiexec", "--bind-to", "core", "--map-by", "socket", "-n", str(nprocs)] + cmd
            
        # 2. Prepare environment (OpenMP tuning)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(config.num_threads)
        env["OMP_PROC_BIND"] = "spread"
        env["OMP_PLACES"] = "cores"

        cmd_str = ' '.join(cmd)

        if verbose:
            print(f"Executing: {cmd_str}")
            subprocess.run(cmd, check=True, env=env)
        else:
            # Discard output
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
            print("Simulation completed.")
