import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Set
from vedo import Plotter, Spheres, Lines, Axes, Box, Cylinder, Cone, Plane, Mesh


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
        self._axes = None
        
        # Add interactive picker callback
        self.plotter.add_callback('LeftButtonPress', self._on_mouse_click)

    def _on_mouse_click(self, event):
        """Identify atom under the mouse click."""
        if not event.actor: return
        
        # We only care about clicks on the atoms (spheres)
        # Vedo's event.actor points to the Spheres object.
        # We find the closest point in that actor to the click.
        p = event.picked3d
        if p is None: return
        
        # Get data from controller
        ts = self.data_ctrl.current_timestep if hasattr(self.data_ctrl, 'current_timestep') else 0
        df = self.data_ctrl.get_atom_data_at_timestep(ts)
        if df.empty: return
        
        # Find nearest atom in the current frame data
        pos = df[['x','y','z']].values
        dist = np.linalg.norm(pos - p, axis=1)
        idx = np.argmin(dist)
        
        # if too far, ignore
        if dist[idx] > (df.iloc[idx].get('diameter', 0.1) * 2):
            return
            
        atom_info = df.iloc[idx]
        atom_id = df.index[idx] if 'id' not in df.columns else atom_info['id']
        mol_id = atom_info.get('mol', 'N/A')
        
        print(f"\n>>> Picked Atom ID: {atom_id}")
        print(f"    Molecule/Chain ID: {mol_id}")
        print(f"    Type: {atom_info.get('type', 'N/A')}")
        print(f"    Position: ({atom_info['x']:.3f}, {atom_info['y']:.3f}, {atom_info['z']:.3f})")
        
        # 2. Trigger Pick Callback for UI
        if hasattr(self, 'on_pick_cb') and self.on_pick_cb:
            try:
                # Handle molecular ID as integer for logic, float for safety
                m_id = int(float(mol_id)) if str(mol_id).replace('.','').isdigit() else mol_id
                info = {
                    'id': atom_id,
                    'mol': m_id,
                    'type': atom_info.get('type', 'N/A'),
                    'pos': (atom_info['x'], atom_info['y'], atom_info['z'])
                }
                self.on_pick_cb(info)
            except: pass

    def update_persistent_bounds(self):
        """Update the max extents of the axes to include particles and all loaded VTKs, then keep them fixed."""
        self._persistent_view_bounds = self.vtk_ctrl.recompute_bounds(self.data_ctrl._init_limits, only_visible=False)
        self.update_axes()

    def update_axes(self):
        """Redraw axes based on current persistent bounds."""
        if self._axes:
            try: self.plotter.remove(self._axes)
            except Exception: pass
            self._axes = None
            
        if self._persistent_view_bounds:
            b = self._persistent_view_bounds
            # Create standard coordinate axes
            self._axes = Axes(
                xrange=(b[0], b[1]), 
                yrange=(b[2], b[3]), 
                zrange=(b[4], b[5]),
                xtitle='X', ytitle='Y', ztitle='Z',
                c='black'
            )
            self.plotter.add(self._axes)
            self.plotter.render()

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

    def show_timestep(self, ts: int):
        """Render a specific timestep from batches or the preview dataset."""
        # Check for Ghost/Preview Mode first
        if self.data_ctrl.is_preview_mode and self.data_ctrl.preview_df is not None:
            # For restarts, we treat the whole df as the single 'frame'
            frame_data = self.data_ctrl.preview_df
            current_ts = 0 # Dummy TS for preview
        else:
            frame_data = self.data_ctrl.get_atom_data_at_timestep(ts)
            current_ts = ts
            
        if frame_data is None or frame_data.empty:
            return
            
        # 1. Clear previous dynamic actors
        for act in self._dynamic_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._dynamic_actors = []
        
        # 2. Render Atoms
        pos = frame_data[['x', 'y', 'z']].values
        dias = frame_data['diameter'].values
        
        # Colors based on type
        types = frame_data['type'].values
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow', 'black']
        
        # Group by (type, diameter) to ensure constant r and c for each vedo actor
        # This works around the vedo limitation where c and r cannot both be sequences.
        for (t, d), group in frame_data.groupby(['type', 'diameter']):
            group_pos = group[['x', 'y', 'z']].values
            r_val = d / 2
            c_val = colors[(int(t)-1)%len(colors)]
            
            spheres = Spheres(group_pos, r=r_val, c=c_val)
            self.plotter.add(spheres)
            self._dynamic_actors.append(spheres)
        
        # 3. Update Highlights
        self.hl_ctrl.render_highlights(frame_data)
        
        # Sync current rendering state back to the app if needed
        # (Usually handled via callbacks in the controller)
        self.plotter.render()

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
        if bounds is None or len(bounds) < 6:
            return []
            
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
