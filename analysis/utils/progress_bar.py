import tkinter as tk
from tkinter import ttk

class ProgressBar(tk.Toplevel):
    """
    A reusable progress bar window that can be modal (blocking interaction with parent).
    """
    def __init__(self, parent, title="Processing", message="Please wait...", modal=True):
        super().__init__(parent)
        self.title(title)
        self.parent = parent
        self.modal = modal
        
        # Center the window
        self.geometry("300x120")
        self.resizable(False, False)
        self._center_window()
        
        # Setup UI
        self.protocol("WM_DELETE_WINDOW", lambda: None) # Prevent closing
        
        main_frame = tk.Frame(self, padx=20, pady=20)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        self.label = tk.Label(main_frame, text=message, wraplength=260)
        self.label.pack(pady=(0, 10))
        
        self.progress = ttk.Progressbar(main_frame, orient=tk.HORIZONTAL, length=260, mode='determinate')
        self.progress.pack(pady=5)
        
        self.status_label = tk.Label(main_frame, text="0%", font=('TkDefaultFont', 8), fg='gray')
        self.status_label.pack()
        
        if self.modal:
            self.transient(parent)
            self._set_parent_state(parent, tk.DISABLED)
            self.grab_set()
            # On some systems we might need to wait for visibility before grab
            
    def _set_parent_state(self, widget, state):
        """Recursively set the state of all child widgets."""
        for child in widget.winfo_children():
            # Skip the progress bar itself if it's a child (though it's Toplevel here)
            if child == self:
                continue
            
            try:
                # Try to set the state. Some widgets like Frames/Labels might fail 
                # or not have a visible 'disabled' state, but we try anyway.
                child.configure(state=state)
            except tk.TclError:
                pass
            
            # Recurse into children of this widget (e.g. Frames)
            self._set_parent_state(child, state)

    def _center_window(self):
        self.withdraw()
        self.update_idletasks()
        
        parent_x = self.parent.winfo_rootx()
        parent_y = self.parent.winfo_rooty()
        parent_width = self.parent.winfo_width()
        parent_height = self.parent.winfo_height()
        
        my_width = self.winfo_width()
        my_height = self.winfo_height()
        
        x = parent_x + (parent_width // 2) - (my_width // 2)
        y = parent_y + (parent_height // 2) - (my_height // 2)
        
        self.geometry(f"+{x}+{y}")
        self.deiconify()

    def set_progress(self, value, maximum=100):
        """Update progress bar value (0 to maximum)."""
        self.progress['maximum'] = maximum
        self.progress['value'] = value
        percent = int((value / maximum) * 100) if maximum > 0 else 0
        self.status_label.config(text=f"{percent}%")
        self.update()

    def set_text(self, text):
        """Update the message text."""
        self.label.config(text=text)
        self.update()

    def finish(self):
        """Close the progress bar and release grab."""
        if self.modal:
            self._set_parent_state(self.parent, tk.NORMAL)
            self.grab_release()
        self.destroy()
