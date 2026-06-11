import os
import sys
import pandas as pd
from typing import List, Dict, Optional
from Pulse.pulse_core import PBSManager

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class PipelineEngine:
    # Key: Stage Name
    # Value: Dict of dependencies, description, and relative output file path
    RUN_PIPELINE = {
        "load_atoms": {
            "dependencies": [],
            "description": "Loads raw dumps and creates simplified atomic coordinate cache.",
            "output_file": "load_atoms.parquet"
        },
        "time_duration": {
            "dependencies": ["load_atoms"],
            "description": "Finds the timestep and duration at which the hopper becomes completely empty.",
            "output_file": "time_duration.parquet"
        },
        "mass_flow_rate": {
            "dependencies": ["load_atoms"],
            "description": "Calculates number of beads and total mass in the hopper per timestep.",
            "output_file": "mass_flow_rate.parquet"
        }
    }
    
    GLOBAL_PIPELINE = {
        "pool_time_durations": {
            "dependencies": ["time_duration"], # Depends on time_duration stage of selected runs
            "description": "Pools run-level time duration outcomes into a global parameter sweep.",
            "output_file": "time_durations.parquet"
        }
    }
    
    @staticmethod
    def get_status_path(run_path: str) -> str:
        return os.path.join(run_path, "results_pipeline", "pipeline_status.parquet")
        
    @staticmethod
    def load_status(run_path: str) -> pd.DataFrame:
        status_path = PipelineEngine.get_status_path(run_path)
        stages = list(PipelineEngine.RUN_PIPELINE.keys())
        if os.path.exists(status_path):
            try:
                df = pd.read_parquet(status_path)
                df.set_index('stage', inplace=True)
                for s in stages:
                    if s not in df.index:
                        df.loc[s] = ["not_started", "", ""]
                df = df.loc[df.index.intersection(stages)]
                return df
            except:
                pass
                
        df = pd.DataFrame([
            {"stage": s, "status": "not_started", "job_id": "", "error": ""}
            for s in stages
        ]).set_index('stage')
        return df

    @staticmethod
    def save_status(df: pd.DataFrame, run_path: str):
        status_path = PipelineEngine.get_status_path(run_path)
        os.makedirs(os.path.dirname(status_path), exist_ok=True)
        df.reset_index().to_parquet(status_path, index=False)

    @staticmethod
    def get_global_status_path() -> str:
        return os.path.join(ROOT_DIR, "Results_Pipeline", "pipeline_status.parquet")

    @staticmethod
    def load_global_status() -> pd.DataFrame:
        status_path = PipelineEngine.get_global_status_path()
        stages = list(PipelineEngine.GLOBAL_PIPELINE.keys())
        if os.path.exists(status_path):
            try:
                df = pd.read_parquet(status_path)
                df.set_index('stage', inplace=True)
                for s in stages:
                    if s not in df.index:
                        df.loc[s] = ["not_started", "", ""]
                df = df.loc[df.index.intersection(stages)]
                return df
            except:
                pass
                
        df = pd.DataFrame([
            {"stage": s, "status": "not_started", "job_id": "", "error": ""}
            for s in stages
        ]).set_index('stage')
        return df

    @staticmethod
    def save_global_status(df: pd.DataFrame):
        status_path = PipelineEngine.get_global_status_path()
        os.makedirs(os.path.dirname(status_path), exist_ok=True)
        df.reset_index().to_parquet(status_path, index=False)

    @staticmethod
    def get_run_pipeline_status(run_path: str, active_pbs_job_ids: set) -> Dict:
        """Loads status and detects dead PBS jobs for a run."""
        df = PipelineEngine.load_status(run_path)
        changed = False
        stages = list(PipelineEngine.RUN_PIPELINE.keys())
        
        status_dict = {}
        for stage in stages:
            st = df.loc[stage, "status"]
            job_id = df.loc[stage, "job_id"]
            
            if st in ["pending", "running"] and pd.notna(job_id) and str(job_id).strip():
                clean_jid = str(job_id).strip()
                if clean_jid not in active_pbs_job_ids:
                    st = "failed"
                    df.loc[stage, "status"] = "failed"
                    df.loc[stage, "error"] = "PBS job died or ended without updating status."
                    changed = True
                    
            status_dict[stage] = {
                "status": st,
                "job_id": "" if pd.isna(job_id) else str(job_id),
                "error": "" if pd.isna(df.loc[stage, "error"]) else str(df.loc[stage, "error"])
            }
            
        if changed:
            PipelineEngine.save_status(df, run_path)
            
        return status_dict

    @staticmethod
    def get_global_pipeline_status(active_pbs_job_ids: set) -> Dict:
        """Loads global status and detects dead PBS jobs."""
        df = PipelineEngine.load_global_status()
        changed = False
        stages = list(PipelineEngine.GLOBAL_PIPELINE.keys())
        
        status_dict = {}
        for stage in stages:
            st = df.loc[stage, "status"]
            job_id = df.loc[stage, "job_id"]
            
            if st in ["pending", "running"] and pd.notna(job_id) and str(job_id).strip():
                clean_jid = str(job_id).strip()
                if clean_jid not in active_pbs_job_ids:
                    st = "failed"
                    df.loc[stage, "status"] = "failed"
                    df.loc[stage, "error"] = "Global PBS job died or ended without updating status."
                    changed = True
                    
            status_dict[stage] = {
                "status": st,
                "job_id": "" if pd.isna(job_id) else str(job_id),
                "error": "" if pd.isna(df.loc[stage, "error"]) else str(df.loc[stage, "error"])
            }
            
        if changed:
            PipelineEngine.save_global_status(df)
            
        return status_dict

    @staticmethod
    def get_downstream_stages(start_stage: str) -> List[str]:
        """Traverses the DAG to find all stages depending directly or transitively on start_stage."""
        downstream = []
        queue = [start_stage]
        while queue:
            curr = queue.pop(0)
            for stage, config in PipelineEngine.RUN_PIPELINE.items():
                if curr in config["dependencies"] and stage not in downstream:
                    downstream.append(stage)
                    queue.append(stage)
        return downstream

    @staticmethod
    def submit_stage_job(run_path: str, stage_name: str, depend_job_id: str = None) -> str:
        """Submits a pipeline stage as a PBS job, supporting dependencies."""
        run_name = os.path.basename(run_path)
        job_name = f"P_{stage_name[:6]}_{run_name[-5:]}"
        pipeline_dir = os.path.join(run_path, "results_pipeline")
        os.makedirs(pipeline_dir, exist_ok=True)
        log_path = os.path.join(pipeline_dir, f"job_{stage_name}.log")
        
        if os.path.exists(log_path):
            try: os.remove(log_path)
            except: pass
            
        job_script = f"""#!/bin/bash
#PBS -N {job_name}
#PBS -q workq
#PBS -l nodes=master:ppn=2
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o {log_path}

cd $PBS_O_WORKDIR
export PYTHONPATH=$PYTHONPATH:$PBS_O_WORKDIR

if [ -f /home/guest/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/guest/miniconda3/etc/profile.d/conda.sh
    conda activate gchain
fi

python Pulse/pipeline/run_stage.py --run '{run_path}' --stage '{stage_name}'
"""
        script_path = os.path.join(pipeline_dir, f"temp_{stage_name}.pbs")
        with open(script_path, "w") as f:
            f.write(job_script)
            
        extra_args = []
        if depend_job_id:
            extra_args = ["-W", f"depend=afterok:{depend_job_id}"]
            
        try:
            job_id = PBSManager.submit_job(script_path, extra_args)
            
            df_status = PipelineEngine.load_status(run_path)
            df_status.loc[stage_name, "status"] = "pending"
            df_status.loc[stage_name, "job_id"] = job_id
            df_status.loc[stage_name, "error"] = ""
            PipelineEngine.save_status(df_status, run_path)
            
            return job_id
        except Exception as e:
            print(f"Failed to submit pipeline job: {e}")
            df_status = PipelineEngine.load_status(run_path)
            df_status.loc[stage_name, "status"] = "failed"
            df_status.loc[stage_name, "error"] = str(e)
            PipelineEngine.save_status(df_status, run_path)
            raise e
        finally:
            if os.path.exists(script_path):
                try: os.remove(script_path)
                except: pass

    @staticmethod
    def submit_global_stage_job(stage_name: str, run_paths: List[str]) -> str:
        """Submits the global pipeline stage as a PBS job."""
        global_dir = os.path.join(ROOT_DIR, "Results_Pipeline")
        os.makedirs(global_dir, exist_ok=True)
        log_path = os.path.join(global_dir, f"job_{stage_name}.log")
        
        if os.path.exists(log_path):
            try: os.remove(log_path)
            except: pass
            
        runs_arg = " ".join([f"'{p}'" for p in run_paths])
        
        job_script = f"""#!/bin/bash
#PBS -N Global_Pool
#PBS -q workq
#PBS -l nodes=master:ppn=2
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -o {log_path}

cd $PBS_O_WORKDIR
export PYTHONPATH=$PYTHONPATH:$PBS_O_WORKDIR

if [ -f /home/guest/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/guest/miniconda3/etc/profile.d/conda.sh
    conda activate gchain
fi

python Pulse/pipeline/run_global_stage.py --stage '{stage_name}' --runs {runs_arg}
"""
        script_path = os.path.join(global_dir, f"temp_{stage_name}.pbs")
        with open(script_path, "w") as f:
            f.write(job_script)
            
        try:
            job_id = PBSManager.submit_job(script_path)
            
            df_status = PipelineEngine.load_global_status()
            df_status.loc[stage_name, "status"] = "pending"
            df_status.loc[stage_name, "job_id"] = job_id
            df_status.loc[stage_name, "error"] = ""
            PipelineEngine.save_global_status(df_status)
            
            return job_id
        except Exception as e:
            print(f"Failed to submit Global PBS job: {e}")
            df_status = PipelineEngine.load_global_status()
            df_status.loc[stage_name, "status"] = "failed"
            df_status.loc[stage_name, "error"] = str(e)
            PipelineEngine.save_global_status(df_status)
            raise e
        finally:
            if os.path.exists(script_path):
                try: os.remove(script_path)
                except: pass

    @staticmethod
    def trigger_pipeline(run_path: str, start_stage: str) -> Dict[str, str]:
        """Triggers execution from the start_stage, automatically queueing dependencies topologically."""
        to_execute = [start_stage] + PipelineEngine.get_downstream_stages(start_stage)
        
        execution_order = []
        visited = set()
        
        def visit(stage):
            if stage in visited:
                return
            visited.add(stage)
            for dep in PipelineEngine.RUN_PIPELINE[stage]["dependencies"]:
                if dep in to_execute:
                    visit(dep)
            execution_order.append(stage)
            
        for s in to_execute:
            visit(s)
            
        job_ids = {}
        for stage in execution_order:
            dep_job_ids = []
            for dep in PipelineEngine.RUN_PIPELINE[stage]["dependencies"]:
                if dep in job_ids:
                    dep_job_ids.append(job_ids[dep])
            
            if dep_job_ids:
                depend_str = ":".join(dep_job_ids)
                job_id = PipelineEngine.submit_stage_job(run_path, stage, depend_job_id=depend_str)
            else:
                job_id = PipelineEngine.submit_stage_job(run_path, stage)
                
            job_ids[stage] = job_id
            
        return job_ids
