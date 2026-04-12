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

            if command == 'create_box':
                self._parse_create_box(parts)
            elif command == 'boundary':
                self.boundary = parts[1:4]
            elif command == 'timestep':
                self._parse_timestep(parts)

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

    def _parse_timestep(self, parts):
        # timestep dt
        if len(parts) >= 2:
            try:
                self.dt = float(parts[1])
            except ValueError:
                self.dt = 0.001 # Default

    def _parse_create_box(self, parts):
        # create_box N region-ID
        # Find the region ID argument
        # It's usually the 3rd argument: create_box 1 my_region
        if len(parts) >= 3:
            self.box_region_id = parts[2]

    def get_geometry(self):
        return {
            'box_region': self.box_region_id,
            'boundary': self.boundary,
            'dt': getattr(self, 'dt', 0.001)
        }
