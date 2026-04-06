import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Set
from vedo import Plotter, Spheres, Lines, Axes, Box, Cylinder, Cone, Plane, Mesh

class HighlightController:
    """Manages highlight state and rendering of highlighted atoms."""
    _HL_COLORS = {'atom': 'yellow', 'bond': 'orange', 'angle': 'violet', 'chain': 'lime'}
    
    def __init__(self, plotter: Plotter, get_data_cb, get_cs_cb):
        self.plotter = plotter
        self.get_data_cb = get_data_cb
        self.get_cs_cb = get_cs_cb
        self._highlighted_ids: Set = set()
        self._highlight_color: str = 'yellow'
        self._highlight_actors: List = []

    def clear(self):
        self._highlighted_ids = set()
        for act in self._highlight_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._highlight_actors = []

    def apply(self, mode: str, n: int):
        df, ts = self.get_data_cb()
        if df is None or ts is None:
            return False, "Data not available"
        
        subset = df[df['timestep'] == ts]
        cs = self.get_cs_cb()
        
        if mode == 'atom':
            ids, label = {n}, f'Atom {n}'
        elif mode == 'bond':
            nb = cs - 1
            if nb <= 0: return False, "Chain size must be ≥ 2"
            ci, binc = (n-1)//nb, (n-1)%nb
            cf = ci * cs + 1
            a1, a2 = cf + binc, cf + binc + 1
            ids, label = {a1, a2}, f'Bond {n} (atoms {a1}, {a2})'
        elif mode == 'angle':
            na = cs - 2
            if na <= 0: return False, "Chain size must be ≥ 3"
            ci, ainc = (n-1)//na, (n-1)%na
            cf = ci * cs + 1
            a1, a2, a3 = cf + ainc, cf + ainc + 1, cf + ainc + 2
            ids, label = {a1, a2, a3}, f'Angle {n} (atoms {a1}, {a2}, {a3})'
        elif mode == 'chain':
            if 'mol' in subset.columns:
                ids = set(subset[subset['mol'] == n]['id'].values.tolist())
            else:
                start = (n-1)*cs + 1
                ids = set(range(start, start + cs))
            label = f'Chain {n} ({len(ids)} atoms)'
        else: return False, "Unknown mode"

        if not ids: return False, "No IDs found"
        self._highlighted_ids = ids
        self._highlight_color = self._HL_COLORS.get(mode, 'yellow')
        return True, label

    def render_highlights(self, subset_df: pd.DataFrame):
        for act in self._highlight_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._highlight_actors = []
        
        if not self._highlighted_ids or subset_df.empty: return
        
        hl_data = subset_df[subset_df['id'].isin(self._highlighted_ids)]
        if hl_data.empty: return
        
        pos = hl_data[['x','y','z']].values
        dias = hl_data['diameter'].values * 1.1 # slightly larger
        
        sph = Spheres(pos, r=dias/2, c=self._highlight_color, res=12)
        self.plotter.add(sph)
        self._highlight_actors.append(sph)


class SimulationRenderer:
    """Manages the 3D scene, rendering of particles, chains, axes, and geometry actors."""
    def __init__(self, plotter: Plotter, data_ctrl, vtk_ctrl, hl_ctrl, get_ui_settings_cb):
        self.plotter = plotter
        self.data_ctrl = data_ctrl
        self.vtk_ctrl = vtk_ctrl
        self.hl_ctrl = hl_ctrl
        self.get_ui_settings_cb = get_ui_settings_cb
        self._dynamic_actors = []
        self._persistent_view_bounds = None

    def update_persistent_bounds(self):
        """Update the max extents of the axes to include particles and all loaded VTKs, then keep them fixed."""
        self._persistent_view_bounds = self.vtk_ctrl.recompute_bounds(self.data_ctrl._init_limits, only_visible=False)

    def fit_view(self):
        if self.plotter:
            self.plotter.reset_camera()
            self.plotter.render()

    def reset_view(self):
        self.fit_view()

    def clear(self):
        for act in self._dynamic_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._dynamic_actors = []

    def show_timestep(self, timestep):
        # 1. Fetch data for this specific frame
        frame_data = self.data_ctrl.get_atom_data_at_timestep(timestep)
        
        # 2. Check if data is available
        if frame_data is None or frame_data.empty:
            print(f"Timestep {timestep} not loaded. Requesting batch...")
            self.data_ctrl.request_batch_for_timestep(timestep)
            return

        # 3. Clear existing dynamic visuals
        for act in self._dynamic_actors: 
            try: self.plotter.remove(act)
            except Exception: pass
        self._dynamic_actors = []

        # 4. Generate and display new actors
        actors = self.plot_chain_data(frame_data)
        
        if self._persistent_view_bounds is None:
            self.update_persistent_bounds()
        
        b = self._persistent_view_bounds
        axes = Axes(xrange=(b[0], b[1]), yrange=(b[2], b[3]), zrange=(b[4], b[5]),
                    xtitle='X', ytitle='Y', ztitle='Z', c='black')

        self._dynamic_actors.extend(actors)
        self._dynamic_actors.append(axes)
        
        for act in self._dynamic_actors:
            self.plotter.add(act)

        # Geometry
        settings = self.get_ui_settings_cb()
        if settings.get('show_geometry'):
            geo_actors = self._build_geometry_actors(self._persistent_view_bounds)
            for g_act in geo_actors:
                self.plotter.add(g_act)
                self._dynamic_actors.append(g_act)
        
        if self.hl_ctrl:
            self.hl_ctrl.render_highlights(frame_data)
        self.plotter.render()

    def plot_chain_data(self, subset_df: pd.DataFrame):
        if subset_df.empty: return []
        actors = []
        pos = subset_df[['x','y','z']].values
        dias = subset_df['diameter'].values
        
        spheres = Spheres(pos, r=dias/2, c='gray', res=10)
        actors.append(spheres)
        return actors

    def _build_geometry_actors(self, bounds):
        if not self.data_ctrl.geometry_data: return []
        regions = self.data_ctrl.geometry_data.get('regions', {})
        box_id = self.data_ctrl.geometry_data.get('box_region')
        actors = []
        colors = ['red', 'green', 'blue', 'gold', 'cyan', 'magenta']
        color_idx = 0
        max_span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        plane_size = max_span if max_span > 0 else 1.0

        for r_id, r_data in regions.items():
            style = r_data.get('style')
            params = r_data.get('params', {})
            if style == 'union': continue

            is_box = (r_id == box_id)
            color = 'black' if is_box else colors[color_idx % len(colors)]
            if not is_box: color_idx += 1

            try:
                if style == 'block':
                    xlo, xhi = params['xlo'], params['xhi']
                    ylo, yhi = params['ylo'], params['yhi']
                    zlo, zhi = params['zlo'], params['zhi']
                    center = ((xlo + xhi) / 2.0, (ylo + yhi) / 2.0, (zlo + zhi) / 2.0)
                    dims = (xhi - xlo, yhi - ylo, zhi - zlo)
                    actor = Box(pos=center, length=dims[0], width=dims[1], height=dims[2]).c(color).alpha(0.10)
                    if is_box: actor = actor.wireframe().lw(1).alpha(1.0)
                    actors.append(actor)
                elif style == 'cylinder':
                    dim, c1, c2 = params['dim'], params['c1'], params['c2']
                    radius, lo, hi = params['radius'], params['lo'], params['hi']
                    if dim == 'x': pos, axis = ((lo + hi) / 2.0, c1, c2), (1, 0, 0)
                    elif dim == 'y': pos, axis = (c1, (lo + hi) / 2.0, c2), (0, 1, 0)
                    else: pos, axis = (c1, c2, (lo + hi) / 2.0), (0, 0, 1)
                    actors.append(Cylinder(pos=pos, r=radius, height=abs(hi - lo), axis=axis).c(color).alpha(0.10))
                elif style == 'cone':
                    dim, c1, c2 = params['dim'], params['c1'], params['c2']
                    radlo, radhi = params['radlo'], params['radhi']
                    lo, hi = params['lo'], params['hi']
                    avg_r = max((radlo + radhi) / 2.0, 1e-6)
                    if dim == 'x': pos, axis = ((lo + hi) / 2.0, c1, c2), (1, 0, 0)
                    elif dim == 'y': pos, axis = (c1, (lo + hi) / 2.0, c2), (0, 1, 0)
                    else: pos, axis = (c1, c2, (lo + hi) / 2.0), (0, 0, 1)
                    actors.append(Cone(pos=pos, axis=axis, r=avg_r, height=abs(hi - lo)).c(color).alpha(0.10))
                elif style == 'plane':
                    px, py, pz = params['px'], params['py'], params['pz']
                    nx, ny, nz = params['nx'], params['ny'], params['nz']
                    actors.append(Plane(pos=(px, py, pz), normal=(nx, ny, nz), s=(plane_size, plane_size)).c(color).alpha(0.10))
            except Exception:
                continue
        return actors
