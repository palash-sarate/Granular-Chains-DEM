import pandas as pd
import numpy as np
from typing import Optional, List, Tuple, Set
from vedo import Plotter, Spheres
from analysis.utilities import parse_range_spec

class HighlightController:
    """Manages highlight state and rendering of highlighted atoms."""
    _HL_COLORS = {'atom': 'red', 'bond': '#39FF14', 'angle': 'violet', 'chain': 'black'}
    
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

    def apply(self, mode: str, spec: str):
        df, ts = self.get_data_cb()
        if df is None or ts is None:
            return False, "Data not available"
        
        try:
            target_ids = parse_range_spec(spec)
        except Exception:
            return False, f"Invalid range spec: {spec}"

        if not target_ids:
            self.clear()
            return True, "Selection cleared"

        try:
            # Handle MultiIndex (timestep, id)
            if 'timestep' in df.index.names:
                subset = df.xs(ts, level='timestep')
            else:
                subset = df[df['timestep'] == ts]
        except Exception:
            subset = df
            
        cs = self.get_cs_cb()
        ids = set()
        
        for n in target_ids:
            if mode == 'atom':
                ids.add(int(n))
            elif mode == 'bond':
                nb = cs - 1
                if nb <= 0: continue
                ci, binc = (int(n)-1)//nb, (int(n)-1)%nb
                cf = ci * cs + 1
                a1, a2 = cf + binc, cf + binc + 1
                ids.update({a1, a2})
            elif mode == 'angle':
                na = cs - 2
                if na <= 0: continue
                ci, ainc = (int(n)-1)//na, (int(n)-1)%na
                cf = ci * cs + 1
                a1, a2, a3 = cf + ainc, cf + ainc + 1, cf + ainc + 2
                ids.update({a1, a2, a3})
            elif mode == 'chain':
                if 'mol' in subset.columns:
                    m_matches = subset[subset['mol'].astype(float).astype(int) == int(n)]
                    if not m_matches.empty:
                        if 'id' in m_matches.columns:
                            ids.update(m_matches['id'].astype(float).astype(int).values.tolist())
                        else:
                            if hasattr(m_matches.index, 'names') and 'id' in m_matches.index.names:
                                raw_ids = m_matches.index.get_level_values('id')
                            else:
                                raw_ids = m_matches.index
                            ids.update(raw_ids.astype(float).astype(int).tolist())
                else:
                    start = (int(n)-1)*cs + 1
                    ids.update(range(start, start + cs))

        if not ids:
            return False, f"No atoms found for {mode}(s) {spec}"
            
        self._highlighted_ids = ids
        self._highlight_color = self._HL_COLORS.get(mode, 'yellow')
        label = f"{mode.capitalize()}s {spec} ({len(ids)} atoms)"
        return True, label

    def render_highlights(self, subset_df: pd.DataFrame):
        for act in self._highlight_actors:
            try: self.plotter.remove(act)
            except Exception: pass
        self._highlight_actors = []
        
        if not self._highlighted_ids or subset_df.empty: return
        
        # Match IDs in dataframe (columns or index)
        if 'id' in subset_df.columns:
            hl_data = subset_df[subset_df['id'].isin(self._highlighted_ids)]
        elif subset_df.index.name == 'id':
            hl_data = subset_df[subset_df.index.isin(self._highlighted_ids)]
        else:
            try:
                hl_data = subset_df[subset_df.index.get_level_values('id').isin(self._highlighted_ids)]
            except:
                hl_data = subset_df[subset_df.index.isin(self._highlighted_ids)]
                
        if hl_data.empty: return
        
        pos = hl_data[['x','y','z']].values
        dias = hl_data.get('diameter', 0.01).values * 1.15
        
        sph = Spheres(pos, r=dias/2, c=self._highlight_color, res=12)
        sph.name = "highlight"
        self.plotter.add(sph)
        self._highlight_actors.append(sph)
