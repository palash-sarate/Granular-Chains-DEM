import os
import re
from typing import Optional, Set

class LammpsParser:
    def __init__(self, script_path):
        self.script_path = script_path
        self.regions = {}
        self.variables = {}
        self.box_region_id = None
        self.boundary = ['p', 'p', 'p'] # Default
        self.parse()

    def parse(self):
        visited: Set[str] = set()
        self._parse_file(self.script_path, visited)

    def _parse_file(self, script_path: str, visited: Set[str]) -> None:
        if not script_path:
            return

        # Use absolute paths for cycle protection.
        abs_path = os.path.abspath(script_path)
        if abs_path in visited:
            return
        visited.add(abs_path)

        if not os.path.isfile(abs_path):
            print(f"LammpsParser: include target not found: {script_path}")
            return

        including_dir = os.path.dirname(abs_path)

        with open(abs_path, 'r') as f:
            lines = f.readlines()

        for line in lines:
            line = line.split('#')[0].strip()  # Remove comments
            if not line:
                continue

            parts = line.split()
            command = parts[0]

            if command == 'region':
                self._parse_region(parts)
            elif command == 'create_box':
                self._parse_create_box(parts)
            elif command == 'boundary':
                self.boundary = parts[1:4]
            # elif command == 'variable':
            #     self._parse_variable(parts)
            # elif command == 'include':
            #     include_target = self._resolve_include_target(parts[1:], including_dir)
            #     if include_target:
            #         self._parse_file(include_target, visited)

    def _resolve_include_target(self, include_tokens, including_dir: str) -> Optional[str]:
        """Resolve `include` path tokens to a concrete file path.

        Supports:
        - `include ${var}` where `var` is defined via `variable ... string|equal ...`
        - relative paths (try relative to including file dir, then CWD)
        - quoted paths (single/double)
        """
        if not include_tokens:
            return None

        raw_target = " ".join(include_tokens).strip()
        if len(raw_target) >= 2 and raw_target[0] == raw_target[-1] and raw_target[0] in ("'", '"'):
            raw_target = raw_target[1:-1]

        # Expand ${var} references using current known variables.
        # If a variable is unknown, we skip this include rather than guessing.
        pattern = re.compile(r"\$\{([^}]+)\}")
        unresolved: Set[str] = set()

        def repl(match: re.Match) -> str:
            name = match.group(1)
            if name not in self.variables:
                unresolved.add(name)
                return ""
            return str(self.variables[name])

        resolved = pattern.sub(repl, raw_target)
        if unresolved:
            print(f"LammpsParser: unresolved include variable(s) {sorted(unresolved)} in: {raw_target}")
            return None

        # If still empty, skip.
        resolved = resolved.strip()
        if not resolved:
            return None

        # If absolute and exists, use it.
        if os.path.isabs(resolved):
            if os.path.isfile(resolved):
                return resolved
            return None

        # Relative: first try relative to including file directory.
        candidate = os.path.abspath(os.path.join(including_dir, resolved))
        if os.path.isfile(candidate):
            return candidate

        # Fallback: relative to current working directory.
        candidate_cwd = os.path.abspath(os.path.join(os.getcwd(), resolved))
        if os.path.isfile(candidate_cwd):
            return candidate_cwd

        # Nothing found.
        return None

    def _parse_variable(self, parts):
        # variable name style args...
        # e.g. variable amp equal 0.005
        # e.g. variable amp string 0.005
        if len(parts) >= 4:
            name = parts[1]
            style = parts[2]
            value = " ".join(parts[3:])
            
            # Try to convert to float if possible
            try:
                if style in ['equal', 'string']:
                    # Remove quotes if present
                    value = value.replace('"', '').replace("'", "")
                    # Check if it's a number
                    float_val = float(value)
                    self.variables[name] = float_val
                else:
                    self.variables[name] = value
            except ValueError:
                self.variables[name] = value

    def _parse_region(self, parts):
        # region ID style args...
        r_id = parts[1]
        style = parts[2]
        args = parts[3:]
        
        # Handle 'move' or 'rotate' keywords which might appear at the end
        # We will ignore dynamic movement for static visualization for now, 
        # or just take the initial position.
        # LAMMPS syntax: region ID style args keyword arg ...
        
        # Basic parsing based on style
        params = {}
        if style == 'block':
            # xlo xhi ylo yhi zlo zhi
            params = {
                'xlo': float(args[0]), 'xhi': float(args[1]),
                'ylo': float(args[2]), 'yhi': float(args[3]),
                'zlo': float(args[4]), 'zhi': float(args[5])
            }
        elif style == 'cylinder':
            # dim c1 c2 radius lo hi
            params = {
                'dim': args[0],
                'c1': float(args[1]), 'c2': float(args[2]),
                'radius': float(args[3]),
                'lo': float(args[4]), 'hi': float(args[5])
            }
        elif style == 'cone':
            # dim c1 c2 radlo radhi lo hi
            params = {
                'dim': args[0],
                'c1': float(args[1]), 'c2': float(args[2]),
                'radlo': float(args[3]), 'radhi': float(args[4]),
                'lo': float(args[5]), 'hi': float(args[6])
            }
        elif style == 'plane':
            # px py pz nx ny nz
            params = {
                'px': float(args[0]), 'py': float(args[1]), 'pz': float(args[2]),
                'nx': float(args[3]), 'ny': float(args[4]), 'nz': float(args[5])
            }
        elif style == 'sphere':
            # x y z radius
            params = {
                'x': float(args[0]), 'y': float(args[1]), 'z': float(args[2]),
                'radius': float(args[3])
            }
        elif style == 'union':
            # n reg-ID1 reg-ID2 ...
            count = int(args[0])
            sub_regions = args[1:1+count]
            params = {'sub_regions': sub_regions}
        
        self.regions[r_id] = {'style': style, 'params': params}

    def _parse_create_box(self, parts):
        # create_box N region-ID
        # Find the region ID argument
        # It's usually the 3rd argument: create_box 1 my_region
        if len(parts) >= 3:
            self.box_region_id = parts[2]

    def get_geometry(self):
        return {
            'regions': self.regions,
            'box_region': self.box_region_id,
            'boundary': self.boundary
        }
