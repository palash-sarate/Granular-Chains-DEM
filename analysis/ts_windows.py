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

from analysis.utilities import get_distance_series, get_angle_series, get_xyz_series


# ── Palette & channel constants ────────────────────────────────────────────────

_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
           '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

ATOM_CHANNELS = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz']

_CH_LS = {
    'x': '-',  'y': '--',  'z': ':',
    'vx': '-', 'vy': '--', 'vz': ':',
    'fx': '-', 'fy': '--', 'fz': ':',
}


# ── Range-spec parser ──────────────────────────────────────────────────────────

def parse_range_spec(spec: str) -> list:
    """Parse '1-3,5,7-9,123' → sorted unique list of ints.

    Same syntax as PDF page-range selection.
    Raises ValueError on bad input.
    """
    result = set()
    for part in spec.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            a, b = part.split('-', 1)
            result.update(range(int(a.strip()), int(b.strip()) + 1))
        else:
            result.add(int(part))
    return sorted(result)


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
        self.lb = tk.Listbox(li, height=3, width=52, activestyle='none')
        sb = tk.Scrollbar(li, orient=tk.VERTICAL, command=self.lb.yview)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.lb.config(yscrollcommand=sb.set)
        self.lb.pack(fill=tk.X)

        # Embedded matplotlib figure
        import matplotlib.figure as _mpf
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                        NavigationToolbar2Tk)
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
    TITLE     = 'Bond Distances'
    YLABEL    = 'Distance (m)'

    def _resolve(self, n: int):
        """Map bond number n → (a1, a2, label) using chain-size convention."""
        cs = max(2, int(self.cs_var.get()))
        nb = cs - 1
        ci = (n - 1) // nb
        bi = (n - 1) % nb
        cf = ci * cs + 1
        a1 = cf + bi
        a2 = a1 + 1
        return a1, a2, f'Bond {n}  (a{a1}\u2013a{a2})'

    def _plot_series(self, ids):
        ok, skip = [], []
        for i, n in enumerate(ids):
            try:
                a1, a2, lbl = self._resolve(n)
                data = get_distance_series(self.df_mi, a1, a2)
                if not data:
                    skip.append(n)
                    continue
                ts, vs = zip(*data)
                self.ax.plot(ts, vs, label=lbl, color=self._col(i))
                self.lb.insert(tk.END, lbl)
                ok.append(n)
            except Exception:
                skip.append(n)
        return ok, skip


# ── Angle Values ───────────────────────────────────────────────────────────────

class AnglePlotWindow(_TsWin):
    SAVE_ATTR = '_angle_range_spec'
    TITLE     = 'Angle Values'
    YLABEL    = 'Angle (\u00b0)'

    def _resolve(self, n: int):
        """Map angle number n → (a1, a2, a3, label) using chain-size convention."""
        cs = max(3, int(self.cs_var.get()))
        na = cs - 2
        ci = (n - 1) // na
        ai = (n - 1) % na
        cf = ci * cs + 1
        a1 = cf + ai
        a2 = a1 + 1
        a3 = a1 + 2
        return a1, a2, a3, f'Angle {n}  (a{a1}\u2013a{a2}\u2013a{a3})'

    def _plot_series(self, ids):
        ok, skip = [], []
        for i, n in enumerate(ids):
            try:
                a1, a2, a3, lbl = self._resolve(n)
                data = get_angle_series(self.df_mi, a1, a2, a3)
                if not data:
                    skip.append(n)
                    continue
                ts, vs = zip(*data)
                self.ax.plot(ts, vs, label=lbl, color=self._col(i))
                self.lb.insert(tk.END, lbl)
                ok.append(n)
            except Exception:
                skip.append(n)
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
        tk.Label(fr, text='Channels:').pack(side=tk.LEFT)
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
