import tkinter as tk
from tkinter import ttk, filedialog
import threading
import inspect
import sys
import os
from typing import Dict, Any, Callable
from simulation.orchestrator import SimulationOrchestrator

class SimulationLauncherController:
    def __init__(self, parent_frame: tk.Frame, orchestrator: SimulationOrchestrator, on_simulation_started: Callable = None):
        self.parent = parent_frame
        self.orchestrator = orchestrator
        self.on_simulation_started = on_simulation_started
        
        self.widgets: Dict[str, tk.Widget] = {}
        self.arg_vars: Dict[str, tk.Variable] = {}
        self.current_method = None
        
        self._setup_ui()
        
    def _setup_ui(self):
        tk.Label(self.parent, text='=== Simulation Control ===', font=('Arial', 10, 'bold')).pack(fill=tk.X, pady=(0, 10))
        
        # Method Selection
        tk.Label(self.parent, text='Select Simulation:').pack(anchor='w')
        self.method_var = tk.StringVar()
        methods = self._get_discoverable_methods()
        self.method_combo = ttk.Combobox(self.parent, textvariable=self.method_var, values=list(methods.keys()), state='readonly')
        self.method_combo.pack(fill=tk.X, pady=(0, 10))
        self.method_combo.bind('<<ComboboxSelected>>', self._on_method_selected)
        
        # Dynamic Arguments Container
        self.args_frame = tk.LabelFrame(self.parent, text="Parameters")
        self.args_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Run Button
        self.run_btn = tk.Button(self.parent, text='🚀 Launch Simulation', bg='#e1f5fe', command=self.run_simulation)
        self.run_btn.pack(fill=tk.X, pady=(10, 0))
        
        # Status
        self.status_label = tk.Label(self.parent, text='Ready', fg='gray')
        self.status_label.pack(fill=tk.X)

    def _get_discoverable_methods(self):
        """Finds all methods in the orchestrator excluding private ones."""
        methods = {}
        for name, func in inspect.getmembers(self.orchestrator, predicate=inspect.ismethod):
            if not name.startswith('_'):
                methods[name] = func
        return methods

    def _on_method_selected(self, event=None):
        method_name = self.method_var.get()
        method = getattr(self.orchestrator, method_name)
        self._build_form(method)
        
    def _build_form(self, method):
        # Clear existing
        for widget in self.args_frame.winfo_children():
            widget.destroy()
        self.arg_vars.clear()
        
        sig = inspect.signature(method)
        
        for name, param in sig.parameters.items():
            if name == 'self': continue
            
            row = tk.Frame(self.args_frame)
            row.pack(fill=tk.X, pady=2)
            
            # Label
            label_text = name.replace('_', ' ').title()
            tk.Label(row, text=f"{label_text}:", width=12, anchor='w').pack(side=tk.LEFT)
            
            # Determine Widget by type
            default_val = param.default if param.default is not inspect.Parameter.empty else ""
            
            # Best effort type detection
            if isinstance(default_val, bool) or param.annotation == bool:
                var = tk.BooleanVar(value=bool(default_val))
                widget = tk.Checkbutton(row, variable=var)
                widget.pack(side=tk.LEFT)
            elif name == 'num_procs' and (default_val is None or default_val == inspect.Parameter.empty):
                # Default to max cores for the UI, can be dialed down
                detected_cores = os.cpu_count() or 1
                var = tk.StringVar(value=str(detected_cores))
                widget = tk.Entry(row, textvariable=var, width=20)
                widget.pack(side=tk.LEFT, fill=tk.X, expand=True)
            elif name.endswith('_dir') or name.endswith('_path'):
                var = tk.StringVar(value=str(default_val))
                widget = tk.Entry(row, textvariable=var, width=15)
                widget.pack(side=tk.LEFT, fill=tk.X, expand=True)
                btn = tk.Button(row, text="...", width=2, command=lambda v=var: self._browse(v))
                btn.pack(side=tk.LEFT)
            else:
                var = tk.StringVar(value=str(default_val))
                widget = tk.Entry(row, textvariable=var, width=20)
                widget.pack(side=tk.LEFT, fill=tk.X, expand=True)
            
            self.arg_vars[name] = var

    def _browse(self, var):
        path = filedialog.askdirectory() if 'dir' in var.get() or not var.get() else filedialog.askopenfilename()
        if path:
            var.set(path)

    def run_simulation(self):
        method_name = self.method_var.get()
        if not method_name: return
        
        method = getattr(self.orchestrator, method_name)
        
        # Collect args and coerce types based on defaults
        args = {}
        sig = inspect.signature(method)
        for name, var in self.arg_vars.items():
            val = var.get()
            # Coerce to match signature default type
            default = sig.parameters[name].default
            if default is not inspect.Parameter.empty and default is not None:
                try:
                    if isinstance(default, int): val = int(val)
                    elif isinstance(default, float): val = float(val)
                    elif isinstance(default, list): val = [x.strip() for x in val.split(',')]
                except: pass
            args[name] = val
            
        self.run_btn.config(state=tk.DISABLED, text="⌛ Running...")
        self.status_label.config(text=f"Running {method_name}...", fg='blue')
        
        # Use threading to keep UI alive
        thread = threading.Thread(target=self._execute_thread, args=(method, args))
        thread.daemon = True
        thread.start()

    def _execute_thread(self, method, args):
        try:
            # Capture stdout/stderr? Maybe later. For now just run.
            method(**args)
            self.parent.after(0, lambda: self._on_finish(True, "Simulation Completed"))
        except Exception as e:
            self.parent.after(0, lambda: self._on_finish(False, str(e)))

    def _on_finish(self, success, msg):
        self.run_btn.config(state=tk.NORMAL, text="🚀 Launch Simulation")
        if success:
            self.status_label.config(text=msg, fg='green')
            if self.on_simulation_started:
                self.on_simulation_started()
        else:
            self.status_label.config(text=f"Error: {msg}", fg='red')
            tk.messagebox.showerror("Simulation Error", msg)
