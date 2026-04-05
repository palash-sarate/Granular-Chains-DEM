import os
import tkinter as tk
from tkinter import messagebox
from typing import Dict, Optional, Callable
from analysis.ts_windows import (BondPlotWindow, AnglePlotWindow, AtomPlotWindow,
                                   LeptonPlotWindow)

class AnalysisToolManager:
    """Manages the lifecycle and visibility of secondary analysis windows."""
    def __init__(self, root: tk.Misc, data_ctrl, get_ui_vars_cb: Callable):
        self.root = root
        self.data_ctrl = data_ctrl
        self.get_ui_vars_cb = get_ui_vars_cb
        
        self.windows: Dict[str, Optional[tk.Toplevel]] = {
            'bond': None,
            'angle': None,
            'atom': None,
            'lepton': None
        }
        self.saved_ranges: Dict[str, str] = {
            'bond': '',
            'angle': '',
            'atom': ''
        }

    def open_bond_win(self):
        if self._check_data():
            if self.windows['bond'] is None or not self.windows['bond'].winfo_exists():
                vars = self.get_ui_vars_cb()
                self.windows['bond'] = BondPlotWindow(
                    self.root, self.data_ctrl.df_bonds, vars['chain_size'], self.saved_ranges['bond']
                )
            else: self.windows['bond'].lift()

    def open_angle_win(self):
        if self._check_data():
            if self.windows['angle'] is None or not self.windows['angle'].winfo_exists():
                vars = self.get_ui_vars_cb()
                self.windows['angle'] = AnglePlotWindow(
                    self.root, self.data_ctrl.df_angles, vars['chain_size'], self.saved_ranges['angle']
                )
            else: self.windows['angle'].lift()

    def open_atom_win(self):
        if self._check_data():
            if self.windows['atom'] is None or not self.windows['atom'].winfo_exists():
                vars = self.get_ui_vars_cb()
                self.windows['atom'] = AtomPlotWindow(
                    self.root, self.data_ctrl.df_mi, vars['chain_size'], self.saved_ranges['atom']
                )
            else: self.windows['atom'].lift()

    def open_lepton_win(self):
        if self.windows['lepton'] is not None and self.windows['lepton'].winfo_exists():
            self.windows['lepton'].lift()
            return
        path = self._resolve_lepton_path()
        if path:
            try: 
                self.windows['lepton'] = LeptonPlotWindow(self.root, path)
            except Exception as e: 
                messagebox.showerror("Error", f"Failed to open Lepton window: {e}")
        else: 
            messagebox.showerror("Error", "Could not locate lepton.inc in simulation folder or templates.")

    def _check_data(self) -> bool:
        if self.data_ctrl.df_mi is None:
            messagebox.showinfo('No data', 'Load a simulation first.')
            return False
        return True

    def _resolve_lepton_path(self) -> Optional[str]:
        if self.data_ctrl.current_sim_folder:
            p = os.path.join(self.data_ctrl.current_sim_folder, 'lepton.inc')
            if os.path.exists(p): return p
        p = os.path.join(os.getcwd(), 'simulation_templates', 'lepton.inc')
        return p if os.path.exists(p) else None

    def refresh_windows(self):
        if self.windows['bond'] and self.windows['bond'].winfo_exists():
            self.windows['bond'].refresh_df(self.data_ctrl.df_bonds)
        if self.windows['angle'] and self.windows['angle'].winfo_exists():
            self.windows['angle'].refresh_df(self.data_ctrl.df_angles)
        if self.windows['atom'] and self.windows['atom'].winfo_exists():
            self.windows['atom'].refresh_df(self.data_ctrl.df_mi)
