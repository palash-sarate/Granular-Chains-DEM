import os
import glob
from tkinter import filedialog, messagebox
import tkinter as tk
from typing import Optional, Dict, Callable

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
        ok, err, num_queued = self.data_ctrl.load_data_file(path)
        if not ok:
            messagebox.showerror('Error', err)
        else:
            self.renderer_ctrl.clear()
            self.renderer_ctrl.update_persistent_bounds()
            self.ui_callbacks['on_load_success']()

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
                    self.ui_callbacks['on_load_success'](num_queued)
        
        for v in vtks:
            self.ui_callbacks['add_vtk'](v)
