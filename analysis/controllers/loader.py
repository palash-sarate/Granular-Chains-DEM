import os
import glob
import subprocess
import tempfile
import shutil
from tkinter import filedialog, messagebox, simpledialog
import tkinter as tk
from typing import Optional, Dict, Callable
import pandas as pd

class SimulationLoader:
    """Orchestrates file/folder loading logic and OS dialogs."""
    def __init__(self, root: tk.Misc, data_ctrl, playback_ctrl, vtk_ctrl, renderer_ctrl, ui_callbacks: Dict[str, Callable]):
        self.root = root
        self.data_ctrl = data_ctrl
        self.playback_ctrl = playback_ctrl
        self.vtk_ctrl = vtk_ctrl
        self.renderer_ctrl = renderer_ctrl
        self.ui_callbacks = ui_callbacks

    def open_dump_folder(self):
        folder = filedialog.askdirectory(title='Select simulation folder')
        if not folder: return
        self.load_simulation_folder(folder)

    def load_simulation_folder(self, folder: str, force_reload: bool = False):
        if not folder: return
        self.playback_ctrl.pause()
        pref = self.ui_callbacks.get('enable_preloading', True)
        if hasattr(pref, 'get'): 
            pref = pref.get()
            
        ok, err, num_queued = self.data_ctrl.load_folder(folder, force_reload=force_reload, enable_preloading=pref)
        if not ok:
            # If it failed to load as a LAMMPS sim, check if it's a folder of VTKs
            vtks = glob.glob(os.path.join(folder, "*.vtk"))
            if vtks:
                for v in sorted(vtks):
                    self.ui_callbacks['add_vtk'](v)
                return
            messagebox.showerror('Error', err)
            return
        
        self.renderer_ctrl.clear()
        self.renderer_ctrl.update_persistent_bounds()
        self.ui_callbacks['on_load_success'](num_queued)

    def open_dump_files(self):
        files = filedialog.askopenfilenames(title='Select dump files', filetypes=[('Dump files', 'dump*'), ('All', '*.*')])
        if not files: return
        self.playback_ctrl.pause()
        ok, err, num_queued = self.data_ctrl.load_dump_files(list(files))
        if not ok:
            messagebox.showerror('Error', err)
        else:
            self.renderer_ctrl.clear()
            self.renderer_ctrl.update_persistent_bounds()
            self.ui_callbacks['on_load_success']()

    def open_data_file(self):
        f = filedialog.askopenfilename(title='Select data file', filetypes=[('Data files', '*.data'), ('All', '*.*')])
        if not f: return
        self.open_data_file_path(f)

    def open_data_file_path(self, path: str):
        self.playback_ctrl.pause()
        
        manual_map = None
        if messagebox.askyesno("Column Mapping", "Do you want to manually specify column indices for this data file?"):
            manual_map = self.ui_callbacks['ask_column_mapping'](path)
            if manual_map is None: # User cancelled
                return

        ok, err, num_queued = self.data_ctrl.load_data_file(path, manual_map=manual_map)
        if not ok:
            messagebox.showerror('Error', err)
        else:
            self.renderer_ctrl.clear()
            self.renderer_ctrl.update_persistent_bounds()
            self.ui_callbacks['on_load_success'](clear_vtk=False)

    def handle_dropped_files(self, filenames: list):
        if not filenames: return
        
        vtks = [f for f in filenames if f.lower().endswith('.vtk')]
        others = [f for f in filenames if not f.lower().endswith('.vtk')]
        
        if len(filenames) == 1 and os.path.isdir(filenames[0]):
            self.load_simulation_folder(filenames[0])
            return

        if others:
            if len(others) == 1 and others[0].endswith('.data'):
                self.open_data_file_path(others[0])
            else:
                self.playback_ctrl.pause()
                ok, err, num_queued = self.data_ctrl.load_dump_files(others)
                if not ok: messagebox.showerror('Error', err)
                else:
                    self.renderer_ctrl.clear()
                    self.renderer_ctrl.update_persistent_bounds()
                    self.ui_callbacks['on_load_success'](num_queued, clear_vtk=False)
        
        for v in vtks:
            self.ui_callbacks['add_vtk'](v)
            
        restarts = [f for f in filenames if f.lower().endswith('.bin')]
        if restarts:
            self.open_restart_editor(restarts[0])

    def open_restart_editor(self, path: str):
        """UI entry point for the Restart Editor."""
        self.playback_ctrl.pause()
        self.ui_callbacks['open_restart_editor'](path)

    def run_lammps_script(self, script_content: str, working_dir: str, num_procs: Optional[int] = None):
        """Helper to run a LAMMPS script using the project's SimulationRunner logic with MPI support.
        """
        import os
        from simulation import SimulationRunner
        in_file = os.path.join(working_dir, "in.temp_mod")
        with open(in_file, "w") as f:
            f.write(script_content)
            
        runner = SimulationRunner(lammps_executable="lmp")
        
        # Determine processor count
        nprocs = num_procs if num_procs is not None else os.cpu_count()
        
        # Build command with MPI if needed
        base_cmd = [runner.lammps_exe, "-in", "in.temp_mod"]
        cmd = base_cmd
        if nprocs and nprocs > 1:
            cmd = ["mpiexec", "-n", str(nprocs)] + base_cmd

        # Optimization: Hiding console on Windows
        startupinfo = None
        if os.name == 'nt':
            import subprocess
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            
        try:
            # We use a direct subprocess.run here to capture the output and use startupinfo
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = "1" # Stick to 1 thread for MPI utilities
            
            result = subprocess.run(cmd, cwd=working_dir, capture_output=True, text=True, startupinfo=startupinfo, timeout=60, env=env)
            return result.returncode == 0, result.stdout + result.stderr
        except Exception as e:
            return False, str(e)

    def visualize_restart_binary(self, restart_path: str):
        """Converts a binary restart to a temporary dump for preview."""
        temp_dir = tempfile.mkdtemp()
        dump_file = os.path.join(temp_dir, "preview.dump")
        
        # Avoid backslashes in f-string expressions for Python 3.11 compatibility
        p_safe = restart_path.replace('\\', '/')
        d_safe = dump_file.replace('\\', '/')
        
        script = f"""
read_restart "{p_safe}"
# Use write_dump to export instantly without setting up neighbors/fixes
write_dump all custom "{d_safe}" id mol type x y z diameter
"""
        ok, log = self.run_lammps_script(script, temp_dir)
        if not ok:
            shutil.rmtree(temp_dir)
            return None, f"LAMMPS Error: {log}"
        
        if not os.path.exists(dump_file):
            shutil.rmtree(temp_dir)
            return None, "LAMMPS failed to generate preview dump."
            
        # Parse it
        from analysis.data_manager import parse_dump_file
        df = parse_dump_file(dump_file)
        shutil.rmtree(temp_dir)
        return df, None

    def apply_restart_modification(self, in_path: str, out_path: str, mol_ids: list, reset_ids: bool = False):
        """Performs deletions and optionally resets IDs before saving a new restart binary."""
        temp_dir = tempfile.mkdtemp()
        
        # Avoid backslashes in f-string expressions
        i_safe = in_path.replace('\\', '/')
        o_safe = out_path.replace('\\', '/')
        
        # Logic to delete if any IDs provided
        del_cmd = ""
        if mol_ids:
            logic = " || ".join([f"mol == {mid}" for mid in mol_ids])
            del_cmd = f"""
variable to_delete atom "{logic}"
group d_grp variable to_delete
delete_atoms group d_grp
"""
        
        # Logic to reset IDs if requested
        reset_cmd = ""
        if reset_ids:
            reset_cmd = """
reset_atoms id
reset_atoms mol all compress yes
"""
            
        script = f"""
read_restart "{i_safe}"
{del_cmd}
{reset_cmd}
write_restart "{o_safe}"
"""
        ok, log = self.run_lammps_script(script, temp_dir)
        shutil.rmtree(temp_dir)
        return ok, log
