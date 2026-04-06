import os
import shutil
import tempfile
import tkinter as tk
from tkinter import messagebox, filedialog
from typing import Optional, Callable

class RestartEditorController:
    """Manages the UI and logic for modifying LAMMPS binary restart files."""
    def __init__(self, root: tk.Misc, data_ctrl, loader_ctrl, renderer, on_close_cb: Callable, sync_cb: Callable):
        self.root = root
        self.data_ctrl = data_ctrl
        self.loader_ctrl = loader_ctrl
        self.renderer = renderer
        self.on_close_cb = on_close_cb
        self.sync_cb = sync_cb
        
        self.source_path: Optional[str] = None
        self.temp_path: Optional[str] = None
        self.selected_mols: set[int] = set()
        
        # UI Elements
        self.panel: Optional[tk.LabelFrame] = None
        self.mol_entry: Optional[tk.Entry] = None
        self._setup_ui()

    def _setup_ui(self):
        self.panel = tk.LabelFrame(self.root, text="Restart Editor", padx=10, pady=10)
        
        header = tk.Frame(self.panel)
        header.pack(fill=tk.X, pady=(0, 5))
        tk.Label(header, text="Molecules to Delete:", font=('Segoe UI', 9, 'bold')).pack(side=tk.LEFT)
        
        self.mol_entry = tk.Entry(self.panel, font=('Consolas', 10))
        self.mol_entry.pack(fill=tk.X, pady=2)
        
        # Sync entry changes back to selected_mols if typed manually
        self.mol_entry.bind('<KeyRelease>', self._on_entry_changed)
        
        tk.Label(self.panel, text="e.g. 5, 12, 18-20", fg='gray', font=('Segoe UI', 8)).pack(anchor='e')
        
        btn_row = tk.Frame(self.panel)
        btn_row.pack(fill=tk.X, pady=(10, 0))
        
        tk.Button(btn_row, text="Preview", command=self.preview, bg='#e3f2fd').pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(btn_row, text="Reset", command=self.reset).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)
        tk.Button(btn_row, text="Save", command=self.save, bg='#e8f5e9').pack(side=tk.LEFT, fill=tk.X, expand=True)
        tk.Button(btn_row, text="Close", command=self.close).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

    def _on_entry_changed(self, event=None):
        """Update internal selection set based on manual typing."""
        from analysis.utilities import parse_range_spec
        spec = self.mol_entry.get().strip()
        try:
            ids = parse_range_spec(spec)
            self.selected_mols = set(ids)
            # Sync with highlighter
            self.sync_cb(spec)
        except: pass

    def toggle_molecule(self, mol_id: int):
        """Toggle a molecule's inclusion in the deletion list."""
        if mol_id in self.selected_mols:
            self.selected_mols.remove(mol_id)
        else:
            self.selected_mols.add(mol_id)
        
        # Update entry
        sorted_ids = sorted(list(self.selected_mols))
        spec = ", ".join(map(str, sorted_ids))
        self.mol_entry.delete(0, tk.END)
        self.mol_entry.insert(0, spec)
        
        # Sync with highlighter
        self.sync_cb(spec)

    def open(self, path: str):
        if not path: return
        self.source_path = path
        self.selected_mols.clear()
        self.mol_entry.delete(0, tk.END)
        self.panel.config(text=f"Editing: {os.path.basename(path)}")
        self.panel.pack(fill=tk.X, pady=5)
        
        # Initial preview
        self.data_ctrl.is_preview_mode = True
        df, err = self.loader_ctrl.visualize_restart_binary(path)
        if err:
            messagebox.showerror("Error", err)
            self.close()
            return
            
        self.data_ctrl.preview_df = df
        self.renderer.show_timestep(0)

    def preview(self):
        if not self.source_path: return
        mol_str = self.mol_entry.get().strip()
        if not mol_str: return
        
        from analysis.utilities import parse_range_spec
        try:
            ids = parse_range_spec(mol_str)
        except ValueError:
            messagebox.showerror("Error", "Invalid IDs. Use: 5, 12, 18-21")
            return
            
        if not ids: return
        
        # Create temp file
        fd, temp_path = tempfile.mkstemp(suffix=".bin")
        os.close(fd)
        self.temp_path = temp_path
        
        ok, log = self.loader_ctrl.apply_restart_modification(self.source_path, temp_path, ids)
        if not ok:
            messagebox.showerror("LAMMPS Error", log)
            return
            
        # Visualize result
        df, err = self.loader_ctrl.visualize_restart_binary(temp_path)
        if err:
            messagebox.showerror("Error", err)
            return
            
        self.data_ctrl.preview_df = df
        self.renderer.show_timestep(0)
        
        # Clear selection after successful deletion
        self.selected_mols.clear()
        self.mol_entry.delete(0, tk.END)
        self.sync_cb("")
        
        messagebox.showinfo("Preview", f"Successfully deleted {len(ids)} molecules. Check viewer.")

    def reset(self):
        if self.temp_path and os.path.exists(self.temp_path):
            try: os.remove(self.temp_path)
            except: pass
        self.mol_entry.delete(0, tk.END)
        self.selected_mols.clear()
        # Trigger highlighter clear
        self.sync_cb("")
        self.open(self.source_path)

    def save(self):
        """Finalizes all deletions and re-indexes the output file for a clean result."""
        # Use temp_path if it exists (meaning the user already clicked Preview)
        # otherwise use source_path (editing the original directly).
        source = self.temp_path if self.temp_path else self.source_path
        
        mol_str = self.mol_entry.get().strip()
        from analysis.utilities import parse_range_spec
        try:
            new_ids = parse_range_spec(mol_str)
        except:
            new_ids = []

        default_name = os.path.basename(self.source_path).replace(".bin", "_mod.bin")
        out_path = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(self.source_path),
            initialfile=default_name,
            filetypes=[("Restart Files", "*.bin")]
        )
        if not out_path: return
        
        # We apply any NEW additions from the box on top of whatever 'source' is.
        # This handles the case where Preview was clicked (clearing the box) 
        # as well as the case where the user goes straight to Save.
        ok, log = self.loader_ctrl.apply_restart_modification(
            source, 
            out_path, 
            new_ids, 
            reset_ids=True
        )
        
        if ok:
            messagebox.showinfo("Success", f"Saved modified and re-indexed restart binary:\n{os.path.basename(out_path)}")
            self.close()
        else:
            messagebox.showerror("Export Error", f"LAMMPS failed to save:\n{log}")

    def close(self):
        self.data_ctrl.is_preview_mode = False
        self.data_ctrl.preview_df = None
        self.panel.pack_forget()
        self.source_path = None
        if self.temp_path and os.path.exists(self.temp_path):
            try: os.remove(self.temp_path)
            except: pass
        self.temp_path = None
        self.on_close_cb()
