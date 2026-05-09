import pyvista as pv
import pandas as pd
import numpy as np
import os
from matplotlib.colors import ListedColormap

class UnifiedRenderer:
    """
    A unified rendering controller for PyVista simulations.
    Supports both interactive Streamlit previews and high-quality off-screen movie generation.
    """
    
    GEOM_COLORS = ["#ADD8E6", "#90EE90", "#FFB6C1", "#FFFFE0", "#E6E6FA", "#F08080", "#AFEEEE", "#FFE4E1"]
    PARTICLE_COLORS = ["#FF5733", "#33FF57", "#3357FF", "#F333FF", "#FF33A8", "#33FFF3", "#F3FF33", "#FF8C00"]

    @staticmethod
    def setup_plotter(window_size=[800, 600], off_screen=False):
        """Creates a standardized PyVista plotter."""
        import pyvista as pv
        plotter = pv.Plotter(window_size=window_size, off_screen=off_screen)
        plotter.background_color = "white"
        return plotter
    @staticmethod
    def parse_lammps_dump(dump_path):
        """Parses a LAMMPS dump file and returns a pandas DataFrame."""
        skip = 0
        cols = []
        try:
            with open(dump_path, 'r') as f:
                for i, line in enumerate(f):
                    if line.startswith("ITEM: ATOMS"):
                        cols = line.strip().split()[2:]
                        skip = i + 1
                        break
            
            # Optimization: Only load the columns we actually need for visualization
            needed_cols = [c for c in cols if c in ['x', 'y', 'z', 'mol', 'diameter']]
            df = pd.read_csv(dump_path, skiprows=skip, sep=r'\s+', names=cols, usecols=needed_cols, engine='c')
            return df
        except Exception:
            return None

    @staticmethod
    def apply_scene(plotter, vtk_files=None, dump_path=None, 
                    show_geometry=True, geom_opacity=0.4,
                    show_particles=True, color_by_mol=True,
                    show_axes=True, show_grid=False,
                    grid_color='gray',
                    camera_state=None,
                    offset=None,
                    zoom=1.0):
        """
        Populates a PyVista plotter with geometry, particles, and scene decorators.
        
        Args:
            plotter: pv.Plotter instance
            vtk_files: List of paths to .vtk geometry files
            dump_path: Path to .dump particle file
            camera_state: Optional camera position/focal point
            zoom: Zoom factor (applied via camera.zoom)
        """
        
        # 1. Add Geometry
        if show_geometry and vtk_files:
            for i, vtk_file in enumerate(vtk_files):
                if os.path.exists(vtk_file):
                    mesh = pv.read(vtk_file)
                    if offset:
                        mesh.translate(offset)
                    color = UnifiedRenderer.GEOM_COLORS[i % len(UnifiedRenderer.GEOM_COLORS)]
                    plotter.add_mesh(mesh, color=color, opacity=geom_opacity, show_edges=False)

        # 2. Add Particles
        if show_particles and dump_path and os.path.exists(dump_path):
            df = UnifiedRenderer.parse_lammps_dump(dump_path)
            if df is not None and not df.empty:
                if all(k in df.columns for k in ['x', 'y', 'z']):
                    points = df[['x', 'y', 'z']].values
                    if offset:
                        points += np.array(offset)
                    
                    pdata = pv.PolyData(points)
                    
                    if 'diameter' in df.columns:
                        pdata['radius'] = df['diameter'].values / 2.0
                    else:
                        pdata['radius'] = np.array([0.001] * len(points))
                    
                    sphere = pv.Sphere(radius=1.0, theta_resolution=8, phi_resolution=8)
                    
                    if color_by_mol and 'mol' in df.columns:
                        m_ids = df['mol'].values
                        pdata['color_idx'] = m_ids % len(UnifiedRenderer.PARTICLE_COLORS)
                        particles = pdata.glyph(scale="radius", geom=sphere, orient=False)
                        my_cmap = ListedColormap(UnifiedRenderer.PARTICLE_COLORS)
                        
                        plotter.add_mesh(
                            particles, 
                            scalars="color_idx", 
                            cmap=my_cmap, 
                            smooth_shading=True,
                            show_scalar_bar=False,
                            categories=True
                        )
                    else:
                        particles = pdata.glyph(scale="radius", geom=sphere, orient=False)
                        plotter.add_mesh(particles, color="#FF5733", smooth_shading=True)

        # 3. Decorators
        if show_axes:
            plotter.add_axes()
            
        if show_grid:
            plotter.show_grid(color=grid_color, font_size=10)

        # 4. Camera
        if camera_state:
            try:
                # Normalize camera state to list of lists (PyVista/stpyvista compatibility)
                norm_cam = [list(x) if isinstance(x, (list, tuple)) else x for x in camera_state]
                plotter.camera_position = norm_cam
            except Exception as e:
                # print(f"DEBUG: Camera restore failed: {e}")
                plotter.view_yz() # Fallback to +X
        else:
            plotter.view_yz() # Default to +X
            
        if zoom != 1.0:
            plotter.camera.zoom(zoom)
            
        return plotter
