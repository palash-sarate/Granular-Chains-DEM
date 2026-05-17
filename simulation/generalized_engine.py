import os
import re
import json
import random
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

from simulation.config import SimulationConfig
from simulation.runner import SimulationRunner
from utilities.env_loader import load_env

# Check for PyYAML availability
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


def parse_simple_yaml(content: str) -> dict:
    """
    A lightweight, robust YAML-to-dict parser in vanilla Python.
    Serves as a reliable fallback if PyYAML is not pre-installed on the host system.
    Handles nested key-values, lists, and list items.
    """
    lines = content.splitlines()
    result = {}
    stack = [result]
    indent_stack = [-1]
    
    for line in lines:
        # Strip comments and trailing whitespace
        line_no_comment = re.sub(r'#.*$', '', line)
        stripped = line_no_comment.strip()
        if not stripped:
            continue
            
        indent = len(line_no_comment) - len(line_no_comment.lstrip())
        
        # Adjust stack based on indentation
        while indent <= indent_stack[-1] and len(stack) > 1:
            stack.pop()
            indent_stack.pop()
            
        current = stack[-1]
        
        # Check for list item e.g. "- target: x" or "- chain_*.dump"
        if stripped.startswith("-"):
            val_str = stripped[1:].strip()
            # Clean enclosing quotes
            if val_str.startswith('"') and val_str.endswith('"'): val_str = val_str[1:-1]
            elif val_str.startswith("'") and val_str.endswith("'"): val_str = val_str[1:-1]
            
            # If item is key-value e.g. "- target: x"
            if ":" in val_str:
                k, v = val_str.split(":", 1)
                k, v = k.strip(), v.strip()
                if v.startswith('"') and v.endswith('"'): v = v[1:-1]
                elif v.startswith("'") and v.endswith("'"): v = v[1:-1]
                
                item_dict = {k: v}
                if isinstance(current, list):
                    current.append(item_dict)
                    stack.append(item_dict)
                    indent_stack.append(indent + 2)
            else:
                if isinstance(current, list):
                    current.append(val_str)
        elif ":" in stripped:
            k, v = stripped.split(":", 1)
            k, v = k.strip(), v.strip()
            if v.startswith('"') and v.endswith('"'): v = v[1:-1]
            elif v.startswith("'") and v.endswith("'"): v = v[1:-1]
            
            # Auto-convert basic types
            if v.lower() == "true":
                v = True
            elif v.lower() == "false":
                v = False
            elif v == "":
                # Check lookahead to see if next block is a list or dict
                is_list = False
                for next_line in lines[lines.index(line)+1:]:
                    n_stripped = re.sub(r'#.*$', '', next_line).strip()
                    if n_stripped:
                        if n_stripped.startswith("-"):
                            is_list = True
                        break
                v = [] if is_list else {}
                if isinstance(current, dict):
                    current[k] = v
                elif isinstance(current, list):
                    current.append({k: v})
                stack.append(v)
                indent_stack.append(indent)
            else:
                try:
                    if "." in v: v = float(v)
                    else: v = int(v)
                except ValueError:
                    pass
                if isinstance(current, dict):
                    current[k] = v
                elif isinstance(current, list):
                    current.append({k: v})
                    
    return result


class SchemaValidationError(Exception):
    """Exception raised when a simulation schema or template fails pre-flight validation."""
    pass


class SchemaValidator:
    """Proactive linter that validates YAML schema layout and matches placeholders in templates."""
    
    @staticmethod
    def lint_template(content: str, bound_variables: set, file_label: str) -> List[str]:
        """Finds any missing or unbound template variables."""
        errors = []
        placeholders = set(re.findall(r'\{\{\s*(\w+)\s*\}\}', content))
        for token in placeholders:
            if token not in bound_variables:
                errors.append(f"Missing parameter '{token}' required by {file_label}")
        return errors

    @classmethod
    def validate_preflight(cls, 
                           schema: dict, 
                           schema_path: str,
                           stage_id: str, 
                           param_overrides: Optional[dict] = None,
                           parent_dir: Optional[str] = None) -> None:
        """
        Validates the entire schema topology and performs a deep lint on target stage templates.
        Raises SchemaValidationError on failure.
        """
        # 1. Topological checks
        sim_type = schema.get("simulation_type")
        if not sim_type:
            raise SchemaValidationError(f"Schema at {schema_path} is missing 'simulation_type'.")

        # Normalize stages (handles both dict and list-of-dicts)
        stages_raw = schema.get("stages", [])
        stages_list = []
        if isinstance(stages_raw, dict):
            for k, v in stages_raw.items():
                if isinstance(v, dict):
                    v_copy = v.copy()
                    v_copy["id"] = k
                    stages_list.append(v_copy)
        elif isinstance(stages_raw, list):
            stages_list = stages_raw
        else:
            raise SchemaValidationError("The 'stages' definition must be a list or a dictionary.")

        # Find target stage
        stage = None
        for s in stages_list:
            if s.get("id") == stage_id:
                stage = s
                break

        if not stage:
            raise SchemaValidationError(f"Stage '{stage_id}' is not defined inside schema at {schema_path}.")

        # 2. Check template declarations
        templates_spec = stage.get("templates", {})
        if not templates_spec:
            # Fallback to flatter keys e.g. template_path, includes
            input_script_tpl = stage.get("template_path")
            includes_list = stage.get("includes", {})
        else:
            input_script_tpl = templates_spec.get("input_script")
            includes_list = templates_spec.get("includes", [])

        if not input_script_tpl:
            raise SchemaValidationError(f"Stage '{stage_id}' does not declare an 'input_script' template.")

        # 3. Consolidate bound variables list
        bound_vars = {
            "outdir", "seed", "parent_restart_path", "resume_file"
        }

        # Global parameters (supports both dict and list format)
        global_params_raw = schema.get("global_parameters", schema.get("global_params", []))
        if isinstance(global_params_raw, dict):
            bound_vars.update(global_params_raw.keys())
        elif isinstance(global_params_raw, list):
            for p in global_params_raw:
                if isinstance(p, dict) and "name" in p:
                    bound_vars.add(p["name"])

        # Stage parameters (supports both dict and list format)
        stage_params_raw = stage.get("parameters", stage.get("params", []))
        if isinstance(stage_params_raw, dict):
            bound_vars.update(stage_params_raw.keys())
        elif isinstance(stage_params_raw, list):
            for p in stage_params_raw:
                if isinstance(p, dict) and "name" in p:
                    bound_vars.add(p["name"])

        # Add manual parameter overrides
        if param_overrides:
            bound_vars.update(param_overrides.keys())

        # Include files trigger automatically bound variables in LAMMPS e.g. includes_lepton_inc
        # Add those beforehand so they are registered as bound variables
        if isinstance(includes_list, dict):
            for k in includes_list.keys():
                bound_vars.add(f"includes_{k}")
        elif isinstance(includes_list, list):
            for inc in includes_list:
                if isinstance(inc, dict):
                    target = inc.get("target") or inc.get("name")
                    if target:
                        var_key = "includes_" + target.replace("/", "_").replace("\\", "_").replace(".", "_")
                        bound_vars.add(var_key)

        # 4. Perform Resume/Branch validation
        mode = stage.get("mode")
        if mode in ["resume", "branch"] or stage.get("parent_stage"):
            if not parent_dir:
                raise SchemaValidationError(f"Stage '{stage_id}' runs in '{mode}' mode but is missing a '--parent_dir'.")

        # 5. Locate and Lint Input Template
        schema_dir = Path(schema_path).resolve().parent
        
        # Clean prefix if written as simulation_templates/
        input_script_tpl_clean = input_script_tpl
        if input_script_tpl.startswith("simulation_templates/"):
            input_script_tpl_clean = input_script_tpl[len("simulation_templates/"):]
        elif input_script_tpl.startswith("simulation_templates\\"):
            input_script_tpl_clean = input_script_tpl[len("simulation_templates\\"):]
            
        input_tpl_path = schema_dir / input_script_tpl_clean
        if not input_tpl_path.exists():
            input_tpl_path = Path("simulation_templates") / input_script_tpl_clean
            if not input_tpl_path.exists():
                raise SchemaValidationError(f"Input script template not found on disk: {input_script_tpl}")

        with open(input_tpl_path, "r", encoding="utf-8") as f:
            input_content = f.read()

        errors = cls.lint_template(input_content, bound_vars, f"input script template '{input_script_tpl}'")

        # 6. Locate and Lint Include Templates
        if isinstance(includes_list, dict):
            for inc_name, inc_tpl in includes_list.items():
                # Allow placeholders inside include paths
                if "{{" in str(inc_tpl):
                    continue
                inc_tpl_str = str(inc_tpl)
                inc_tpl_clean = inc_tpl_str
                if inc_tpl_str.startswith("simulation_templates/"):
                    inc_tpl_clean = inc_tpl_str[len("simulation_templates/"):]
                elif inc_tpl_str.startswith("simulation_templates\\"):
                    inc_tpl_clean = inc_tpl_str[len("simulation_templates\\"):]
                
                inc_tpl_path = schema_dir / inc_tpl_clean
                if not inc_tpl_path.exists():
                    inc_tpl_path = Path("simulation_templates") / inc_tpl_clean
                    
                if inc_tpl_path.exists():
                    with open(inc_tpl_path, "r", encoding="utf-8") as f:
                        inc_content = f.read()
                    errors.extend(cls.lint_template(inc_content, bound_vars, f"include template '{inc_tpl}'"))
        elif isinstance(includes_list, list):
            for inc in includes_list:
                if isinstance(inc, dict):
                    inc_tpl = inc.get("template")
                    if inc_tpl:
                        if "{{" in str(inc_tpl):
                            continue
                        inc_tpl_str = str(inc_tpl)
                        inc_tpl_clean = inc_tpl_str
                        if inc_tpl_str.startswith("simulation_templates/"):
                            inc_tpl_clean = inc_tpl_str[len("simulation_templates/"):]
                        elif inc_tpl_str.startswith("simulation_templates\\"):
                            inc_tpl_clean = inc_tpl_str[len("simulation_templates\\"):]
                        
                        inc_tpl_path = schema_dir / inc_tpl_clean
                        if not inc_tpl_path.exists():
                            inc_tpl_path = Path("simulation_templates") / inc_tpl_clean
                            
                        if inc_tpl_path.exists():
                            with open(inc_tpl_path, "r", encoding="utf-8") as f:
                                inc_content = f.read()
                            errors.extend(cls.lint_template(inc_content, bound_vars, f"include template '{inc_tpl}'"))

        # Raise critical validation failures
        if errors:
            raise SchemaValidationError(
                f"\n❌ Schema Pre-Flight Linter Failed for Stage '{stage_id}':\n" + 
                "\n".join([f"  - {err}" for err in errors])
            )
            
        print("✓ Pre-flight schema validation completed successfully (0 errors found).")


class GeneralizedEngine:
    def __init__(self):
        # Load environment configurations
        self.env = load_env()
        self.lammps_exe = self.env.get("LAMMPS_EXECUTABLE", "lmp")
        self.dumping_yard = self.env.get("DUMPING_YARD", "dumping_yard")
        
    def load_schema(self, schema_path: str) -> dict:
        """Loads a YAML or JSON simulation schema dynamically."""
        path = Path(schema_path)
        if not path.exists():
            raise FileNotFoundError(f"Simulation schema not found at: {schema_path}")
            
        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
            
        if path.suffix.lower() in [".yaml", ".yml"]:
            if HAS_YAML:
                try:
                    return yaml.safe_load(content)
                except Exception as e:
                    print(f"Warning: PyYAML failed parsing {schema_path}, falling back to custom parser. Error: {e}")
                    return parse_simple_yaml(content)
            else:
                return parse_simple_yaml(content)
        else:
            # Assume JSON
            return json.loads(content)

    def render_template(self, template_content: str, variables: Dict[str, Any]) -> str:
        """Replaces {{placeholder}} blocks inside a template with active parameter values."""
        def replace_match(match):
            var_name = match.group(1).strip()
            if var_name in variables:
                val = variables[var_name]
                # Format floats cleanly for LAMMPS
                if isinstance(val, float):
                    return f"{val:.12g}"
                return str(val)
            return match.group(0) # Keep placeholder if not bound
            
        pattern = re.compile(r'\{\{\s*(\w+)\s*\}\}')
        return pattern.sub(replace_match, template_content)

    def run_stage(self, 
                  schema_path: str, 
                  stage_id: str, 
                  run_name: Optional[str] = None, 
                  parent_dir: Optional[str] = None, 
                  param_overrides: Optional[dict] = None, 
                  num_procs: int = 1, 
                  num_threads: int = 1, 
                  use_kokkos: bool = True, 
                  use_intel: bool = True, 
                  inplace: bool = False) -> str:
        """
        Executes a specific pipeline stage of a simulation schema.
        Interpolates parameters, creates output directories, renders files, and launches LAMMPS.
        """
        schema = self.load_schema(schema_path)
        simulation_type = schema.get("simulation_type", "Generalized_Simulation")
        
        # 1. Run Pre-Flight Validation Linter Proactively!
        SchemaValidator.validate_preflight(
            schema=schema,
            schema_path=schema_path,
            stage_id=stage_id,
            param_overrides=param_overrides,
            parent_dir=parent_dir
        )
        
        # 2. Normalize and Locate Stage (handles both dict and list-of-dicts)
        stages_raw = schema.get("stages", [])
        stages_list = []
        if isinstance(stages_raw, dict):
            for k, v in stages_raw.items():
                if isinstance(v, dict):
                    v_copy = v.copy()
                    v_copy["id"] = k
                    stages_list.append(v_copy)
        elif isinstance(stages_raw, list):
            stages_list = stages_raw
            
        stage = None
        for s in stages_list:
            if s.get("id") == stage_id:
                stage = s
                break
                
        if not stage:
            raise ValueError(f"Stage '{stage_id}' not found in simulation schema: {schema_path}")
            
        # 3. Consolidate Parameters (Global defaults -> Stage defaults -> Overrides)
        variables = {}
        
        # Load global parameters
        global_params_raw = schema.get("global_parameters", schema.get("global_params", []))
        if isinstance(global_params_raw, dict):
            variables.update(global_params_raw)
        elif isinstance(global_params_raw, list):
            for p in global_params_raw:
                if isinstance(p, dict) and "name" in p:
                    variables[p["name"]] = p.get("default")
                
        # Load stage parameters
        stage_params_raw = stage.get("parameters", stage.get("params", []))
        if isinstance(stage_params_raw, dict):
            variables.update(stage_params_raw)
        elif isinstance(stage_params_raw, list):
            for p in stage_params_raw:
                if isinstance(p, dict) and "name" in p:
                    variables[p["name"]] = p.get("default")
                
        # Apply parameter overrides
        if param_overrides:
            for k, v in param_overrides.items():
                variables[k] = v
                
        # Ensure a seed exists
        if "seed" not in variables:
            variables["seed"] = int(time.time() * 1000) % 1000000
            
        # Resolve nested placeholders inside the parameters themselves (e.g. data_file containing {{N}})
        # We run it in a loop to support up to 3 levels of nested dependencies!
        for _ in range(3):
            resolved_any = False
            for k, v in list(variables.items()):
                if isinstance(v, str) and "{{" in v:
                    rendered = self.render_template(v, variables)
                    if rendered != v:
                        variables[k] = rendered
                        resolved_any = True
            if not resolved_any:
                break

        # 4. Handle Branching & Resumption Lineages
        parent_real_path = None
        if stage.get("mode") in ["branch", "resume"] or stage.get("parent_stage"):
            if not parent_dir:
                raise ValueError(f"Stage '{stage_id}' requires a '--parent_dir' path to resume/branch.")
                
            parent_real_path = os.path.abspath(parent_dir)
            if not os.path.exists(parent_real_path):
                raise FileNotFoundError(f"Parent directory not found: {parent_real_path}")
                
            # Locate parent stage's outputs from parent run metadata
            parent_metadata_file = os.path.join(parent_real_path, "metadata.json")
            parent_restart_file = None
            
            if os.path.exists(parent_metadata_file):
                with open(parent_metadata_file, "r") as f:
                    parent_meta = json.load(f)
                
                parent_outputs = parent_meta.get("outputs", {})
                restart_rel = parent_outputs.get("restart_file")
                if restart_rel:
                    parent_restart_file = os.path.join(parent_real_path, restart_rel)
            
            # Fallback restart search
            if not parent_restart_file or not os.path.exists(parent_restart_file):
                restart_dir = os.path.join(parent_real_path, "restart")
                if os.path.exists(restart_dir):
                    restarts = sorted(
                        [os.path.join(restart_dir, f) for f in os.listdir(restart_dir) if f.endswith(".bin")],
                        key=os.path.getmtime
                    )
                    if restarts:
                        parent_restart_file = restarts[-1]
                        
            if not parent_restart_file or not os.path.exists(parent_restart_file):
                raise FileNotFoundError(f"Could not locate a valid restart binary file inside: {parent_real_path}")
                
            parent_restart_normalized = parent_restart_file.replace("\\", "/")
            variables["parent_restart_path"] = parent_restart_normalized
            variables["resume_file"] = parent_restart_normalized
            
        # 5. Determine Run Directory
        if inplace and parent_real_path:
            output_dir = parent_real_path
            run_name = os.path.basename(output_dir)
        else:
            if not run_name:
                run_name = f"{simulation_type}_{stage_id}_S{variables['seed']}"
            output_dir = os.path.abspath(os.path.join(self.dumping_yard, simulation_type, run_name))
            
        variables["outdir"] = output_dir.replace("\\", "/")
        
        # 6. Create Stage Directories
        os.makedirs(output_dir, exist_ok=True)
        
        outputs_spec = stage.get("outputs", {})
        
        # Always create logs/restart/dump dirs
        for sub in ["restart", "log"]:
            os.makedirs(os.path.join(output_dir, sub), exist_ok=True)
            
        # Parse dump files outputs to auto-create directories
        for dump_pat in outputs_spec.get("dump_files", []):
            if "/" in dump_pat or "\\" in dump_pat:
                sub_dir = os.path.dirname(dump_pat)
                os.makedirs(os.path.join(output_dir, sub_dir), exist_ok=True)
                
        # 7. Render Input Script and Includes Templates
        templates_spec = stage.get("templates", {})
        if not templates_spec:
            input_script_tpl = stage.get("template_path")
            includes_list = stage.get("includes", {})
        else:
            input_script_tpl = templates_spec.get("input_script")
            includes_list = templates_spec.get("includes", [])
            
        # Resolve templates path
        schema_dir = Path(schema_path).resolve().parent
        
        # Clean prefix if written as simulation_templates/
        input_script_tpl_clean = input_script_tpl
        if input_script_tpl.startswith("simulation_templates/"):
            input_script_tpl_clean = input_script_tpl[len("simulation_templates/"):]
        elif input_script_tpl.startswith("simulation_templates\\"):
            input_script_tpl_clean = input_script_tpl[len("simulation_templates\\"):]
            
        input_tpl_path = schema_dir / input_script_tpl_clean
        
        if not input_tpl_path.exists():
            input_tpl_path = Path("simulation_templates") / input_script_tpl_clean
            if not input_tpl_path.exists():
                raise FileNotFoundError(f"Input script template not found: {input_tpl_path}")
                
        with open(input_tpl_path, "r", encoding="utf-8") as f:
            tpl_content = f.read()
            
        # Interpolate variables inside template paths (if any)
        # Render main input script
        rendered_script = self.render_template(tpl_content, variables)
        
        # Generate script path
        script_suffix = "resume" if stage.get("mode") == "resume" else "run"
        temp_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(output_dir))), "temp")
        os.makedirs(temp_dir, exist_ok=True)
        rendered_script_path = os.path.join(temp_dir, f"in.{stage_id}.{script_suffix}")
        
        # Render target includes files
        if isinstance(includes_list, dict):
            for inc_target_raw, inc_tpl_raw in includes_list.items():
                # Interpolate placeholders in include names and templates (e.g. m{{N}})
                inc_target = self.render_template(str(inc_target_raw), variables)
                inc_tpl = self.render_template(str(inc_tpl_raw), variables)
                
                inc_tpl_clean = inc_tpl
                if inc_tpl.startswith("simulation_templates/"):
                    inc_tpl_clean = inc_tpl[len("simulation_templates/"):]
                elif inc_tpl.startswith("simulation_templates\\"):
                    inc_tpl_clean = inc_tpl[len("simulation_templates\\"):]
                
                inc_tpl_path = schema_dir / inc_tpl_clean
                if not inc_tpl_path.exists():
                    inc_tpl_path = Path("simulation_templates") / inc_tpl_clean
                    if not inc_tpl_path.exists():
                        print(f"Warning: Include template not found: {inc_tpl}")
                        continue
                        
                with open(inc_tpl_path, "r", encoding="utf-8") as f:
                    inc_content = f.read()
                    
                rendered_inc = self.render_template(inc_content, variables)
                target_dest = os.path.join(output_dir, inc_target)
                
                os.makedirs(os.path.dirname(target_dest), exist_ok=True)
                
                with open(target_dest, "w", encoding="utf-8") as f:
                    f.write(rendered_inc)
                print(f"Rendered include target: {inc_target}")
                
                var_key = "includes_" + inc_target.replace("/", "_").replace("\\", "_").replace(".", "_")
                variables[var_key] = target_dest.replace("\\", "/")
                
        elif isinstance(includes_list, list):
            for inc in includes_list:
                if isinstance(inc, dict):
                    inc_target = self.render_template(str(inc.get("target")), variables)
                    inc_tpl = self.render_template(str(inc.get("template")), variables)
                    
                    inc_tpl_clean = inc_tpl
                    if inc_tpl.startswith("simulation_templates/"):
                        inc_tpl_clean = inc_tpl[len("simulation_templates/"):]
                    elif inc_tpl.startswith("simulation_templates\\"):
                        inc_tpl_clean = inc_tpl[len("simulation_templates\\"):]
                    
                    inc_tpl_path = schema_dir / inc_tpl_clean
                    if not inc_tpl_path.exists():
                        inc_tpl_path = Path("simulation_templates") / inc_tpl_clean
                        if not inc_tpl_path.exists():
                            print(f"Warning: Include template not found: {inc_tpl}")
                            continue
                            
                    with open(inc_tpl_path, "r", encoding="utf-8") as f:
                        inc_content = f.read()
                        
                    rendered_inc = self.render_template(inc_content, variables)
                    target_dest = os.path.join(output_dir, inc_target)
                    
                    os.makedirs(os.path.dirname(target_dest), exist_ok=True)
                    
                    with open(target_dest, "w", encoding="utf-8") as f:
                        f.write(rendered_inc)
                    print(f"Rendered include target: {inc_target}")
                    
                    var_key = "includes_" + inc_target.replace("/", "_").replace("\\", "_").replace(".", "_")
                    variables[var_key] = target_dest.replace("\\", "/")
                
        # Re-render main script in case include variables were bound
        rendered_script = self.render_template(rendered_script, variables)
        
        with open(rendered_script_path, "w", encoding="utf-8") as f:
            f.write(rendered_script)
        print(f"Rendered main script to: {rendered_script_path}")
        
        # 8. Write Standardized metadata.json
        metadata = {
            "name": run_name,
            "simulation_type": simulation_type,
            "simulation": simulation_type,
            "stage_id": stage_id,
            "parent": parent_real_path,
            "params": variables,
            "seed": variables["seed"],
            "steps_target": variables.get("relax_steps") or variables.get("run_steps") or 1000000,
            "schema_path": os.path.abspath(schema_path),
            "outputs": outputs_spec,
            "timestamp": datetime.now().isoformat()
        }
        
        metadata_path = os.path.join(output_dir, "metadata.json")
        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=4)
        print(f"Wrote standardized metadata: {metadata_path}")
        
        # 9. Construct SimulationConfig and execute via SimulationRunner
        config = SimulationConfig(
            input_script=rendered_script_path,
            outdir_override=output_dir,
            num_procs=num_procs,
            num_threads=num_threads,
            use_kokkos=use_kokkos,
            use_intel=use_intel,
            simulation=simulation_type,
            run=run_name
        )
        
        runner = SimulationRunner(lammps_executable=self.lammps_exe)
        
        print(f"--- Starting Stage Execution: {simulation_type} => {stage_id} ---")
        try:
            runner.run(config, verbose=True, clean_dir=False, prep_dirs=False)
            print(f"--- Stage Completed Successfully: {stage_id} ---")
        except Exception as e:
            print(f"Error during stage execution: {e}")
            raise e
            
        return output_dir
