import pyvista as pv
import numpy as np
import pandas as pd
import os
from .lammps_parser import LammpsParser
import math

class VTKExporter:
    def __init__(self, dataframe, output_dir, lammps_script=None):
        self.df = dataframe.copy()
        self.output_dir = output_dir
        self.lammps_script = lammps_script
        self.geometry = None
        self.variables = {}
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        if self.lammps_script and os.path.exists(self.lammps_script):
            parser = LammpsParser(self.lammps_script)
            self.geometry = parser.get_geometry()
            self.variables = parser.variables
            print(f"Loaded geometry: {list(self.geometry['regions'].keys())}")
            print(f"Loaded variables: {self.variables.keys()}")

    def export_series(self, start_frame=0, end_frame=None):
        """
        Export the simulation data to a series of .vtm files.
        """
        # Ensure index is reset
        if isinstance(self.df.index, pd.MultiIndex):
            self.df.reset_index(inplace=True)
        elif self.df.index.name in ['timestep', 'id']:
             self.df.reset_index(inplace=True)
             
        timesteps = sorted(self.df['timestep'].unique())
        
        if end_frame is None:
            end_frame = len(timesteps)
            
        selected_timesteps = timesteps[start_frame:end_frame]
        
        print(f"Exporting {len(selected_timesteps)} frames to VTK...")
        
        # Get simulation parameters for dynamic geometry
        amp = self.variables.get('amp', 0.0)
        freq = self.variables.get('freq', 0.0)
        dt = self.variables.get('dt', 1e-6) # Default dt if not found
        
        print(f"Dynamic parameters: amp={amp}, freq={freq}, dt={dt}")

        for i, ts in enumerate(selected_timesteps):
            frame_df = self.df[self.df['timestep'] == ts]
            
            # Calculate current time
            current_time = ts * dt
            
            # Calculate oscillation offset
            # x = amp * sin(2 * pi * freq * time)
            x_offset = 0.0
            if amp != 0 and freq != 0:
                x_offset = amp * math.sin(2 * math.pi * freq * current_time)
            
            self._export_frame(frame_df, ts, x_offset)
            
            if i % 10 == 0:
                print(f"Exported frame {i}/{len(selected_timesteps)} (Timestep {ts})")

    def _export_frame(self, frame_df, timestep, x_offset):
        """
        Export a single frame as a MultiBlock dataset (.vtm)
        """
        multi_block = pv.MultiBlock()
        
        # 1. Particles
        points = frame_df[['x', 'y', 'z']].values
        cloud = pv.PolyData(points)
        
        # Add attributes
        for col in frame_df.columns:
            if col not in ['x', 'y', 'z', 'timestep']:
                cloud.point_data[col] = frame_df[col].values
        
        multi_block.append(cloud, name="Particles")
        
        # 2. Geometry
        if self.geometry:
            regions = self.geometry['regions']
            box_id = self.geometry['box_region']
            
            for r_id, r_data in regions.items():
                style = r_data['style']
                params = r_data['params']
                
                mesh = None
                
                # Apply dynamic movement to hopper parts
                # In in.hopper_flow, hopper_cone and hopper_cyl have 'move'
                current_x_offset = 0.0
                if r_id in ['hopper_cone', 'hopper_cyl', 'hopper_union']:
                    current_x_offset = x_offset
                
                if style == 'block':
                    mesh = self._create_block(params)
                elif style == 'cylinder':
                    mesh = self._create_cylinder(params)
                elif style == 'cone':
                    mesh = self._create_cone(params)
                elif style == 'plane':
                    mesh = self._create_plane(params)
                
                if mesh:
                    # Apply translation
                    if current_x_offset != 0:
                        mesh.translate([current_x_offset, 0, 0], inplace=True)
                    
                    multi_block.append(mesh, name=f"Region_{r_id}")

        # Save
        filename = os.path.join(self.output_dir, f"simulation_{timestep}.vtm")
        multi_block.save(filename)

    def _create_block(self, params):
        bounds = [params['xlo'], params['xhi'], params['ylo'], params['yhi'], params['zlo'], params['zhi']]
        return pv.Cube(bounds=bounds)

    def _create_cylinder(self, params):
        # dim c1 c2 radius lo hi
        # pyvista cylinder is along X axis by default, centered at center
        radius = params['radius']
        height = params['hi'] - params['lo']
        center_axis = (params['lo'] + params['hi']) / 2
        
        # Create cylinder along Z (default is X, so we rotate)
        # direction=(0,0,1) aligns it with Z
        
        if params['dim'] == 'z':
            center = [params['c1'], params['c2'], center_axis]
            direction = [0, 0, 1]
        elif params['dim'] == 'y':
            center = [params['c1'], center_axis, params['c2']]
            direction = [0, 1, 0]
        else: # x
            center = [center_axis, params['c1'], params['c2']]
            direction = [1, 0, 0]
            
        return pv.Cylinder(center=center, direction=direction, radius=radius, height=height, resolution=30)

    def _create_cone(self, params):
        # dim c1 c2 radlo radhi lo hi
        # pyvista Cone is defined by center, direction, height, radius (base)
        # It tapers to a point usually, but we can use a truncated cone logic if needed?
        # PyVista's Cone source creates a cone with a tip. 
        # For a truncated cone (hopper), we might need to use a Cylinder with different top/bottom radii 
        # BUT PyVista Cylinder doesn't support different radii easily in the basic source.
        # We can use pv.Cone if radhi is 0.
        # If it's a truncated cone, we can construct it using pv.CylinderStructured or extruding a circle?
        # Actually, pv.Cone has a 'radius' (base) and 'radius_top' (tip) in newer versions? 
        # Let's check documentation... standard VTK ConeSource has Radius and Resolution.
        # Wait, PyVista's `Cone` function is a helper for `vtkConeSource`.
        # `vtkConeSource` does NOT support truncated cones directly (it makes a pointy cone).
        
        # Better approach for truncated cone: Use `pv.Cylinder`? No.
        # Use `pv.Plotter().add_mesh(pv.Cone(...))`?
        # Let's construct a mesh manually for the cone to be safe and correct.
        
        res = 30
        h = params['hi'] - params['lo']
        
        # Create two circles
        z_lo = params['lo']
        z_hi = params['hi']
        r_lo = params['radlo']
        r_hi = params['radhi']
        
        # Generate points
        theta = np.linspace(0, 2*np.pi, res, endpoint=False)
        x_lo = r_lo * np.cos(theta) + params['c1']
        y_lo = r_lo * np.sin(theta) + params['c2']
        z_lo_arr = np.full_like(theta, z_lo)
        
        x_hi = r_hi * np.cos(theta) + params['c1']
        y_hi = r_hi * np.sin(theta) + params['c2']
        z_hi_arr = np.full_like(theta, z_hi)
        
        # Combine points
        points_lo = np.column_stack((x_lo, y_lo, z_lo_arr))
        points_hi = np.column_stack((x_hi, y_hi, z_hi_arr))
        points = np.vstack((points_lo, points_hi))
        
        # Create faces (quads connecting the two circles)
        # Connectivity
        faces = []
        for i in range(res):
            # Points indices
            p1 = i
            p2 = (i + 1) % res
            p3 = p2 + res
            p4 = p1 + res
            faces.extend([4, p1, p2, p3, p4])
            
        mesh = pv.PolyData(points, faces)
        return mesh

    def _create_plane(self, params):
        # px py pz nx ny nz
        # Create a large plane centered at p, normal n
        center = [params['px'], params['py'], params['pz']]
        normal = [params['nx'], params['ny'], params['nz']]
        return pv.Plane(center=center, direction=normal, i_size=2, j_size=2)

