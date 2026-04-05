import os
from typing import Optional, Dict
from vedo import load as vload, Plotter

class VtkOverlayController:
    """Manages VTK mesh loading, visibility, and scene bounds."""
    def __init__(self, plotter: Plotter):
        self.plotter = plotter
        self.vtk_meshes: Dict[str, dict] = {}
        self.vtk_color_idx = 0
        self._scene_bounds = None
        self._colors = ['red', 'green', 'blue', 'yellow', 'cyan', 'magenta']
        self.master_visible = True

    def set_master_visibility(self, visible: bool):
        self.master_visible = visible
        for name, data in self.vtk_meshes.items():
            act = data['actor']
            if not visible:
                try: self.plotter.remove(act)
                except Exception: pass
            elif data['visible']:
                self.plotter.add(act)
        self.plotter.render()

    def add_mesh(self, path: str):
        try:
            name = os.path.basename(path)
            mesh = vload(path)
            if mesh is None: return False, f"Failed to load {name}"
            color = self._colors[self.vtk_color_idx % len(self._colors)]
            self.vtk_color_idx += 1
            mesh.c(color).alpha(0.5)
            
            if self.master_visible:
                self.plotter.add(mesh)
            
            self.vtk_meshes[name] = {'actor': mesh, 'visible': True, 'path': path}
            return True, name
        except Exception as e:
            return False, str(e)

    def set_visibility(self, name: str, visible: bool):
        if name in self.vtk_meshes:
            self.vtk_meshes[name]['visible'] = visible
            act = self.vtk_meshes[name]['actor']
            if self.master_visible and visible:
                self.plotter.add(act)
            else:
                try: self.plotter.remove(act)
                except Exception: pass
            self.plotter.render()

    def remove_meshes(self, names: list):
        for name in names:
            if name in self.vtk_meshes:
                act = self.vtk_meshes[name]['actor']
                try: self.plotter.remove(act)
                except Exception: pass
                del self.vtk_meshes[name]

    def clear(self):
        for mesh_data in self.vtk_meshes.values():
            try: self.plotter.remove(mesh_data['actor'])
            except Exception: pass
        self.vtk_meshes.clear()
        self.vtk_color_idx = 0
        self._scene_bounds = None

    def recompute_bounds(self, init_limits: Optional[tuple], only_visible=True):
        if init_limits is not None:
            xmin, xmax = init_limits[0]
            ymin, ymax = init_limits[1]
            zmin, zmax = init_limits[2]
        else:
            xmin, xmax = float('inf'), float('-inf')
            ymin, ymax = float('inf'), float('-inf')
            zmin, zmax = float('inf'), float('-inf')

        for mesh_data in self.vtk_meshes.values():
            is_vis = self.master_visible and mesh_data['visible']
            if not only_visible or is_vis:
                try:
                    bnds = mesh_data['actor'].bounds()
                    if len(bnds) == 6:
                        xmin = min(xmin, bnds[0]); xmax = max(xmax, bnds[1])
                        ymin = min(ymin, bnds[2]); ymax = max(ymax, bnds[3])
                        zmin = min(zmin, bnds[4]); zmax = max(zmax, bnds[5])
                except Exception: pass

        if xmin == float('inf'):
            xmin, xmax, ymin, ymax, zmin, zmax = -1, 1, -1, 1, -1, 1

        x_pad = (xmax - xmin) * 0.05 if xmax > xmin else 0.1
        y_pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.1
        z_pad = (zmax - zmin) * 0.05 if zmax > zmin else 0.1
        
        self._scene_bounds = [
            xmin - x_pad, xmax + x_pad,
            ymin - y_pad, ymax + y_pad,
            zmin - z_pad, zmax + z_pad
        ]
        return self._scene_bounds
