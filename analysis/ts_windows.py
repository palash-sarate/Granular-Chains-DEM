"""
analysis/ts_windows.py
Time-series plot windows for the Chains Simulation Viewer.

Three independently closable Toplevel windows:
  BondPlotWindow  – distance vs timestep for selected bonds
  AnglePlotWindow – angle vs timestep for selected angles
  AtomPlotWindow  – x/y/z/vx/vy/vz/fx/fy/fz vs timestep for selected atoms
"""
import tkinter as tk
import pandas as pd
import re
import numpy as np
import os
# Embedded matplotlib figure
import matplotlib.figure as _mpf
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                NavigationToolbar2Tk)

# ── Palette & channel constants ────────────────────────────────────────────────

_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
           '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

ATOM_CHANNELS  = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz']
BOND_CHANNELS  = ['dist', 'energy', 'force']
ANGLE_CHANNELS = ['theta', 'energy']

_CH_LS = {
    'x': '-',  'y': '--',  'z': ':',
    'vx': '-', 'vy': '--', 'vz': ':',
    'fx': '-', 'fy': '--', 'fz': ':',
    # Bonds
    'dist': '-', 'energy': '--', 'force': ':',
    # Angles
    'theta': '-', 'energy': '--'
}

from analysis.utilities import parse_range_spec

# ── Shared base window ─────────────────────────────────────────────────────────

class _TsWin(tk.Toplevel):
    """Base class for Bond/Angle/Atom time-series windows.

    Subclasses must set SAVE_ATTR, TITLE, YLABEL, SHOW_CS and
    override _plot_series(ids) → (plotted_list, skipped_list).
    """
    SAVE_ATTR = ''
    TITLE     = 'Time Series'
    YLABEL    = 'Value'
    SHOW_CS   = True   # whether to show the Chain-size entry

    def __init__(self, parent, df_mi: pd.DataFrame,
                 cs_var: tk.StringVar, initial: str = ''):
        super().__init__(parent)
        self.parent   = parent
        self.df_mi    = df_mi
        self.cs_var   = cs_var
        self.title(self.TITLE)
        self.resizable(True, True)
        self.protocol('WM_DELETE_WINDOW', self._on_close)
        self.range_var = tk.StringVar(value=initial)
        self._build_ui()
        if initial.strip():
            self._do_plot()

    # ── UI construction ────────────────────────────────────────────────────────
    def _build_ui(self):
        # Top row (chain-size for Bond/Angle; channel checkboxes for Atom)
        top = tk.Frame(self)
        top.pack(fill=tk.X, padx=8, pady=(8, 0))
        self._top_controls(top)

        # Range entry row
        rf = tk.Frame(self)
        rf.pack(fill=tk.X, padx=8, pady=(4, 0))
        tk.Label(rf, text='Range:').pack(side=tk.LEFT)
        e = tk.Entry(rf, textvariable=self.range_var, width=24)
        e.pack(side=tk.LEFT, padx=(4, 0))
        e.bind('<Return>', lambda _e: self._do_plot())
        tk.Button(rf, text='Plot',  command=self._do_plot).pack(side=tk.LEFT, padx=(6, 0))
        tk.Button(rf, text='Clear', command=self._clear).pack(side=tk.LEFT, padx=(4, 0))

        # Active-series listbox (informational)
        lf = tk.Frame(self)
        lf.pack(fill=tk.X, padx=8, pady=(4, 0))
        tk.Label(lf, text='Active series:', anchor='w').pack(fill=tk.X)
        li = tk.Frame(lf)
        li.pack(fill=tk.X)
        self.lb = tk.Listbox(li, height=3, width=52, activestyle='none', exportselection=False)
        sb = tk.Scrollbar(li, orient=tk.VERTICAL, command=self.lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lb.config(yscrollcommand=sb.set)
        self.lb.pack(fill=tk.X)


        self.fig = _mpf.Figure(figsize=(8, 4), tight_layout=True)
        self.ax  = self.fig.add_subplot(111)
        self.ax.set_xlabel('Timestep')
        self.ax.set_ylabel(self.YLABEL)
        self.ax.grid(True, alpha=0.3)

        cf = tk.Frame(self)
        cf.pack(fill=tk.BOTH, expand=True, padx=8, pady=(4, 0))
        self.canvas = FigureCanvasTkAgg(self.fig, master=cf)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        tbf = tk.Frame(self)
        tbf.pack(fill=tk.X, padx=8)
        NavigationToolbar2Tk(self.canvas, tbf).update()

        # Status bar
        self.sv = tk.StringVar(value="Enter range (e.g. 1-3,5) and click Plot.")
        tk.Label(self, textvariable=self.sv, fg='gray',
                 anchor='w').pack(fill=tk.X, padx=8, pady=(2, 6))

    def _top_controls(self, fr):
        """Override to add mode-specific top widgets."""
        if self.SHOW_CS:
            tk.Label(fr, text='Chain size:').pack(side=tk.LEFT)
            tk.Entry(fr, textvariable=self.cs_var,
                     width=5).pack(side=tk.LEFT, padx=(4, 0))

    # ── Plot ───────────────────────────────────────────────────────────────────

    def _clear(self):
        self.ax.clear()
        self.ax.set_xlabel('Timestep')
        self.ax.set_ylabel(self.YLABEL)
        self.ax.grid(True, alpha=0.3)
        self.lb.delete(0, tk.END)
        self.canvas.draw()
        self.sv.set('Cleared.')

    def _do_plot(self):
        spec = self.range_var.get().strip()
        if not spec:
            return
        try:
            ids = parse_range_spec(spec)
        except (ValueError, TypeError) as exc:
            self.sv.set(f'Range parse error: {exc}')
            return

        self.ax.clear()
        self.ax.set_xlabel('Timestep')
        self.ax.set_ylabel(self.YLABEL)
        self.ax.set_title(self.TITLE)
        self.ax.grid(True, alpha=0.3)
        self.lb.delete(0, tk.END)

        ok, skip = self._plot_series(ids)
        if ok:
            self.ax.legend(fontsize='small', loc='best')
        msg = f'Plotted {len(ok)} series.'
        if skip:
            msg += f'  Skipped {skip} (no data in df).'
        self.sv.set(msg)
        self.canvas.draw()

    def _plot_series(self, ids: list) -> tuple:
        """Override in subclass. Draw onto self.ax; update self.lb.
        Return (plotted_ids, skipped_ids).
        """
        raise NotImplementedError

    def _col(self, i: int) -> str:
        return _COLORS[i % len(_COLORS)]

    # ── Lifecycle ──────────────────────────────────────────────────────────────

    def refresh_df(self, df_mi: pd.DataFrame):
        """Called by ViewerApp when simulation data reloads."""
        self.df_mi = df_mi
        if self.range_var.get().strip():
            self._do_plot()

    def _on_close(self):
        """Save current range string to parent so it is restored on reopen."""
        if self.SAVE_ATTR:
            setattr(self.parent, self.SAVE_ATTR, self.range_var.get())
        self.destroy()


# ── Bond Distances ─────────────────────────────────────────────────────────────

class BondPlotWindow(_TsWin):
    SAVE_ATTR = '_bond_range_spec'
    TITLE     = 'Bond Properties'
    YLABEL    = 'Value'

    def __init__(self, parent, df_mi, cs_var, initial=''):
        self._ch_vars: dict = {ch: tk.BooleanVar(value=(ch == 'dist')) for ch in BOND_CHANNELS}
        self._lines: dict   = {} # {bond_id: {channel: Line2D}}
        super().__init__(parent, df_mi, cs_var, initial)

    def _top_controls(self, fr):
        super()._top_controls(fr)
        tk.Label(fr, text='  Channels:').pack(side=tk.LEFT)
        for ch in BOND_CHANNELS:
            tk.Checkbutton(fr, text=ch, variable=self._ch_vars[ch],
                           command=self._sync_vis).pack(side=tk.LEFT, padx=2)

    def _sync_vis(self):
        for bond_lines in self._lines.values():
            for ch, ln in bond_lines.items():
                ln.set_visible(self._ch_vars[ch].get())
        self.canvas.draw_idle()

    def _plot_series(self, ids):
        self._lines = {}
        if self.df_mi is None or self.df_mi.empty:
            return [], ids
            
        avail = [ch for ch in BOND_CHANNELS if ch in self.df_mi.columns]
        ok, skip = [], []
        
        for i, bid in enumerate(ids):
            try:
                rows = self.df_mi.xs(bid, level='id')
            except KeyError:
                skip.append(bid)
                continue
            
            if rows.empty:
                skip.append(bid)
                continue
                
            ts = rows.index.values
            col = self._col(i)
            self._lines[bid] = {}
            
            for ch in avail:
                vis = self._ch_vars[ch].get()
                ln, = self.ax.plot(ts, rows[ch].values,
                                   label=f'Bond {bid} {ch}',
                                   color=col,
                                   linestyle=_CH_LS.get(ch, '-'),
                                   visible=vis)
                self._lines[bid][ch] = ln
            
            self.lb.insert(tk.END, f'Bond {bid} ({len(avail)} channels)')
            ok.append(bid)
            
        return ok, skip


# ── Angle Values ───────────────────────────────────────────────────────────────

class AnglePlotWindow(_TsWin):
    SAVE_ATTR = '_angle_range_spec'
    TITLE     = 'Angle Properties'
    YLABEL    = 'Value'

    def __init__(self, parent, df_mi, cs_var, initial=''):
        self._ch_vars: dict = {ch: tk.BooleanVar(value=(ch == 'theta')) for ch in ANGLE_CHANNELS}
        self._lines: dict   = {} # {angle_id: {channel: Line2D}}
        super().__init__(parent, df_mi, cs_var, initial)

    def _top_controls(self, fr):
        super()._top_controls(fr)
        tk.Label(fr, text='  Channels:').pack(side=tk.LEFT)
        for ch in ANGLE_CHANNELS:
            tk.Checkbutton(fr, text=ch, variable=self._ch_vars[ch],
                           command=self._sync_vis).pack(side=tk.LEFT, padx=2)

    def _sync_vis(self):
        for angle_lines in self._lines.values():
            for ch, ln in angle_lines.items():
                ln.set_visible(self._ch_vars[ch].get())
        self.canvas.draw_idle()

    def _plot_series(self, ids):
        self._lines = {}
        if self.df_mi is None or self.df_mi.empty:
            return [], ids
            
        avail = [ch for ch in ANGLE_CHANNELS if ch in self.df_mi.columns]
        ok, skip = [], []
        
        for i, aid in enumerate(ids):
            try:
                rows = self.df_mi.xs(aid, level='id')
            except KeyError:
                skip.append(aid)
                continue
            
            if rows.empty:
                skip.append(aid)
                continue
                
            ts = rows.index.values
            col = self._col(i)
            self._lines[aid] = {}
            
            for ch in avail:
                vis = self._ch_vars[ch].get()
                ln, = self.ax.plot(ts, rows[ch].values,
                                   label=f'Angle {aid} {ch}',
                                   color=col,
                                   linestyle=_CH_LS.get(ch, '-'),
                                   visible=vis)
                self._lines[aid][ch] = ln
            
            self.lb.insert(tk.END, f'Angle {aid} ({len(avail)} channels)')
            ok.append(aid)
            
        return ok, skip


# ── Atom Properties ────────────────────────────────────────────────────────────

class AtomPlotWindow(_TsWin):
    SAVE_ATTR = '_atom_range_spec'
    TITLE     = 'Atom Properties'
    YLABEL    = 'Value'
    SHOW_CS   = False   # atom IDs come directly from the range string

    def __init__(self, parent, df_mi, cs_var, initial=''):
        # Channel BooleanVars must exist BEFORE _build_ui (called by super)
        self._ch_vars: dict = {ch: tk.BooleanVar(value=True) for ch in ATOM_CHANNELS}
        self._lines: dict   = {}   # {atom_id: {channel: Line2D}}
        super().__init__(parent, df_mi, cs_var, initial)

    def _top_controls(self, fr):
        super()._top_controls(fr)
        tk.Label(fr, text='  Channels:').pack(side=tk.LEFT)
        for ch in ATOM_CHANNELS:
            tk.Checkbutton(fr, text=ch, variable=self._ch_vars[ch],
                           command=self._sync_vis).pack(side=tk.LEFT, padx=2)

    def _sync_vis(self):
        """Show/hide individual lines when a checkbox is toggled — no replot."""
        for atom_lines in self._lines.values():
            for ch, ln in atom_lines.items():
                ln.set_visible(self._ch_vars[ch].get())
        self.canvas.draw_idle()

    def _plot_series(self, ids):
        self._lines = {}
        avail = [ch for ch in ATOM_CHANNELS if ch in self.df_mi.columns]
        ok, skip = [], []
        for i, aid in enumerate(ids):
            try:
                rows = self.df_mi.xs(aid, level='id')
            except KeyError:
                skip.append(aid)
                continue
            if rows.empty:
                skip.append(aid)
                continue
            ts  = rows.index.values
            col = self._col(i)
            self._lines[aid] = {}
            for ch in avail:
                if ch not in rows.columns:
                    continue
                vis = self._ch_vars[ch].get()
                ln, = self.ax.plot(ts, rows[ch].values,
                                   label=f'Atom {aid} {ch}',
                                   color=col,
                                   linestyle=_CH_LS.get(ch, '-'),
                                   visible=vis)
                self._lines[aid][ch] = ln
            self.lb.insert(tk.END, f'Atom {aid}  ({len(avail)} channels)')
            ok.append(aid)
        return ok, skip

# ── Lepton Potential Visualization ─────────────────────────────────────────────

class LeptonPlotWindow(tk.Toplevel):
    """
    Parses a lepton.inc file and plots the specified bond and angle potentials.
    Allows real-time editing of both variables and equation strings.
    """
    TITLE = 'Lepton Potentials'

    def __init__(self, parent, lepton_inc_path: str):
        super().__init__(parent)
        self.parent = parent
        self.path = lepton_inc_path
        self.title(f"{self.TITLE} - {os.path.basename(lepton_inc_path) if lepton_inc_path else ''}")
        self.geometry("1100x750")
        self.resizable(True, True)

        # Persistence for original file state
        self._orig_vars = {}
        self._orig_bond_expr = ""
        self._orig_angle_expr = ""

        # UI Vars
        self.var_entries = {}       # {name: (Label, Entry, StringVar)}
        self.bond_expr_var = tk.StringVar()
        self.angle_expr_var = tk.StringVar()

        self._parse_lepton_inc(lepton_inc_path)
        self._build_ui()
        
        # Initial traces - will trigger _do_plot
        self._set_traces(True)
        self._do_plot()

    def _parse_lepton_inc(self, path):
        if not path or not os.path.exists(path):
            return
        
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except:
            return

        var_pattern = re.compile(r'^\s*variable\s+(\w+)\s+equal\s+([\d\.e\-+]+|v_\w+)')
        bond_pattern = re.compile(r'^\s*bond_coeff\s+(\d+)\s+([\d\.e\-+]+)\s+"([^"]+)"')
        angle_pattern = re.compile(r'^\s*angle_coeff\s+(\d+)\s+([\d\.e\-+]+)\s+"([^"]+)"')

        # 1. First pass: Collect global variables
        global_vars = {}
        for line in lines:
            line = line.split('#')[0].strip()
            if not line: continue
            m_var = var_pattern.match(line)
            if m_var:
                name, val_str = m_var.groups()
                if val_str.startswith('v_'):
                    ref_name = val_str[2:]
                    val = global_vars.get(ref_name, 0.0)
                else:
                    try: val = float(val_str)
                    except: val = 0.0
                global_vars[name] = val
        
        # 2. Second pass: Collect expressions and their associated mappings/aliases
        for line in lines:
            line = line.split('#')[0].strip()
            if not line: continue

            m_bond = bond_pattern.match(line)
            if m_bond and not self._orig_bond_expr:
                full_expr = m_bond.group(3)
                parts = full_expr.split(';')
                self._orig_bond_expr = parts[0].strip()
                # Extract aliases: k=v_k
                for p in parts[1:]:
                    if '=' in p:
                        lhs, rhs = [x.strip() for x in p.split('=', 1)]
                        # If rhs is a global var (v_name), link it
                        val = global_vars.get(rhs[2:] if rhs.startswith('v_') else rhs, 0.0)
                        self._orig_vars[lhs] = val
                continue

            m_angle = angle_pattern.match(line)
            if m_angle and not self._orig_angle_expr:
                full_expr = m_angle.group(3)
                parts = full_expr.split(';')
                self._orig_angle_expr = parts[0].strip()
                for p in parts[1:]:
                    if '=' in p:
                        lhs, rhs = [x.strip() for x in p.split('=', 1)]
                        val = global_vars.get(rhs[2:] if rhs.startswith('v_') else rhs, 0.0)
                        self._orig_vars[lhs] = val
                continue

        # Also include any globals that weren't aliased but might be used
        for k, v in global_vars.items():
            if k not in self._orig_vars:
                self._orig_vars[k] = v

        # Set UI vars
        self.bond_expr_var.set(self._orig_bond_expr)
        self.angle_expr_var.set(self._orig_angle_expr)

    def _set_traces(self, enabled: bool):
        if enabled:
            self._trace_ids = []
            tid = self.bond_expr_var.trace_add("write", lambda *a: self._on_change())
            self._trace_ids.append((self.bond_expr_var, tid))
            tid = self.angle_expr_var.trace_add("write", lambda *a: self._on_change())
            self._trace_ids.append((self.angle_expr_var, tid))
            for v_name, (lbl, ent, sv) in self.var_entries.items():
                tid = sv.trace_add("write", lambda *a: self._on_change())
                self._trace_ids.append((sv, tid))
        else:
            for sv, tid in self._trace_ids:
                sv.trace_remove("write", tid)
            self._trace_ids = []

    def _on_change(self):
        # Debounced or immediate update
        self._do_plot()

    def _reset_to_file(self):
        self._set_traces(False)
        for name, val in self._orig_vars.items():
            if name in self.var_entries:
                self.var_entries[name][2].set(str(val))
        self.bond_expr_var.set(self._orig_bond_expr)
        self.angle_expr_var.set(self._orig_angle_expr)
        self._set_traces(True)
        self._do_plot()

    def _build_ui(self):
        # Main PanedWindow for Plot (top) and Editor (bottom)
        pw = tk.PanedWindow(self, orient=tk.VERTICAL)
        pw.pack(fill=tk.BOTH, expand=True)

        # Plot Frame
        plot_f = tk.Frame(pw)
        pw.add(plot_f, height=450)
        
        self.fig = _mpf.Figure(figsize=(10, 5), tight_layout=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_f)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        tbf = tk.Frame(plot_f)
        tbf.pack(fill=tk.X, padx=5, pady=(0, 5))
        NavigationToolbar2Tk(self.canvas, tbf).update()

        # Editor Frame
        edit_f = tk.Frame(pw)
        pw.add(edit_f)

        # Expression rows
        expr_f = tk.LabelFrame(edit_f, text="Equation Strings")
        expr_f.pack(fill=tk.X, padx=10, pady=5)
        
        tk.Label(expr_f, text="Bond Expression:").grid(row=0, column=0, sticky='w', pady=2, padx=5)
        tk.Entry(expr_f, textvariable=self.bond_expr_var, width=80).grid(row=0, column=1, sticky='ew', pady=2, padx=5)
        
        tk.Label(expr_f, text="Angle Expression:").grid(row=1, column=0, sticky='w', pady=2, padx=5)
        tk.Entry(expr_f, textvariable=self.angle_expr_var, width=80).grid(row=1, column=1, sticky='ew', pady=2, padx=5)
        expr_f.columnconfigure(1, weight=1)

        # Variable grid
        var_f = tk.LabelFrame(edit_f, text="Lepton Parameters / Variables")
        var_f.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Scrollable area for variables if many
        v_canvas = tk.Canvas(var_f, height=120)
        v_scroll = tk.Scrollbar(var_f, orient=tk.VERTICAL, command=v_canvas.yview)
        v_grid = tk.Frame(v_canvas)
        v_canvas.create_window((0,0), window=v_grid, anchor='nw')
        v_canvas.configure(yscrollcommand=v_scroll.set)
        
        v_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        v_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        row, col = 0, 0
        names = sorted(self._orig_vars.keys())
        for name in names:
            val = self._orig_vars[name]
            l = tk.Label(v_grid, text=f"{name}:")
            l.grid(row=row, column=col*2, padx=(10,2), pady=2, sticky='e')
            sv = tk.StringVar(value=str(val))
            e = tk.Entry(v_grid, textvariable=sv, width=12)
            e.grid(row=row, column=col*2+1, padx=2, pady=2, sticky='w')
            self.var_entries[name] = (l, e, sv)
            
            col += 1
            if col > 3:
                col = 0
                row += 1
        
        v_grid.update_idletasks()
        v_canvas.config(scrollregion=v_canvas.bbox('all'))

        # Footer
        footer = tk.Frame(edit_f)
        footer.pack(fill=tk.X, padx=10, pady=5)
        tk.Button(footer, text="Reset to File Values", command=self._reset_to_file).pack(side=tk.LEFT)
        self.status_lbl = tk.Label(footer, text="Ready", fg="blue")
        self.status_lbl.pack(side=tk.RIGHT)

    def _to_latex(self, expr):
        """Very basic heuristic to make standard equations look better in titles."""
        if not expr: return ""
        # 1. Neutralize common math formatting characters to avoid partial rendering
        # Just replace ** with ^ (matplotlib handles simple ^ for single-digit exponents)
        ltx = expr.replace('**', '^')
        # We don't add { } because we don't know where they end without a real parser
        # 2. Wrap in $ only if we think it's mostly safe, otherwise return raw
        # If it has characters that mathtext doesn't like outside of math mode, 
        # it will fail in set_title if we don't wrap it.
        return rf"${ltx}$"

    def _evaluate_expression(self, expr, x_vals, x_name):
        # 1. Clean expression from LAMMPSisms
        clean_expr = expr.replace('^', '**')
        # 2. Get current values from UI entries
        current_vars = {}
        for name, (lbl, ent, sv) in self.var_entries.items():
            try: current_vars[name] = float(sv.get())
            except: current_vars[name] = 0.0
        
        # 3. Environment
        safe_names = {
            'abs': np.abs, 'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
            'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
            'asin': np.arcsin, 'acos': np.arccos, 'atan': np.arctan,
            'sinh': np.sinh, 'cosh': np.cosh, 'tanh': np.tanh,
            'pow': np.power, 'PI': np.pi
        }
        env = {**safe_names, **current_vars, x_name: x_vals}
        # Allow use of 'v_name' directly for convenience
        for name, val in current_vars.items():
            env[f"v_{name}"] = val
        
        try:
            res = eval(clean_expr, {"__builtins__": None}, env)
            return res, None
        except Exception as e:
            return np.zeros_like(x_vals), str(e)

    def _do_plot(self):
        self.fig.clear()
        ax_bond = self.fig.add_subplot(1, 2, 1)
        ax_angle = self.fig.add_subplot(1, 2, 2)
        
        err_msg = ""

        # Bond Plot
        b_expr = self.bond_expr_var.get().strip()
        if b_expr:
            # Need a range. Try to find r_max in current vars.
            # Look for 'r_max' aliased or directly
            r_max_val = 0.001
            for key in ['r_max', 'v_r_max']:
                if key in self.var_entries:
                    try: 
                        r_max_val = float(self.var_entries[key][2].get())
                        break
                    except: pass
            if r_max_val <= 0: r_max_val = 0.001
            
            r = np.linspace(1e-12, r_max_val * 2.5, 300)
            u, err = self._evaluate_expression(b_expr, r, 'r')
            if err: err_msg += f"Bond: {err}\n"
            
            ax_bond.plot(r * 1e3, u, color='blue', lw=2)
            # Display equation formatted, but fall back if it fails
            t_str = self._to_latex(b_expr)
            try:
                ax_bond.set_title(f"Bond: {t_str}", fontsize=9)
            except:
                ax_bond.set_title(f"Bond: {b_expr}", fontsize=9)
        else:
            ax_bond.text(0.5, 0.5, "Enter Bond Expression", ha='center')

        ax_bond.set_xlabel("r (mm)")
        ax_bond.set_ylabel("U")
        ax_bond.grid(True, alpha=0.3)

        # Angle Plot
        a_expr = self.angle_expr_var.get().strip()
        if a_expr:
            theta = np.linspace(0, np.pi, 200)
            u, err = self._evaluate_expression(a_expr, theta, 'theta')
            if err: err_msg += f"Angle: {err}\n"
            
            ax_angle.plot(np.degrees(theta), u, color='green', lw=2)
            t_str = self._to_latex(a_expr)
            try:
                ax_angle.set_title(f"Angle: {t_str}", fontsize=9)
            except:
                ax_angle.set_title(f"Angle: {a_expr}", fontsize=9)
        else:
            ax_angle.text(0.5, 0.5, "Enter Angle Expression", ha='center')

        ax_angle.set_xlabel("theta (deg)")
        ax_angle.set_ylabel("U")
        ax_angle.grid(True, alpha=0.3)

        if err_msg:
            self.status_lbl.config(text=f"Error: {err_msg.strip()}", fg="red")
        else:
            self.status_lbl.config(text="Plots updated", fg="blue")

        self.canvas.draw()
