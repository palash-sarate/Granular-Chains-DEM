import subprocess
import re
from typing import List, Dict, Optional

class PBSManager:
    @staticmethod
    def get_jobs(user: Optional[str] = None) -> List[Dict]:
        """Gets a list of jobs, optionally filtered by user."""
        cmd = ["qstat", "-f"]
        if user:
            # Note: qstat -f doesn't always support -u directly in all versions
            # so we fetch all and filter in python for robustness
            pass
            
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                return []
            return PBSManager.parse_qstat_f(result.stdout, user)
        except Exception:
            return []

    @staticmethod
    def parse_qstat_f(output: str, filter_user: Optional[str] = None) -> List[Dict]:
        """Parses the output of qstat -f into a list of dictionaries."""
        jobs = []
        current_job = {}
        
        # Split by "Job Id:" to get individual job blocks
        blocks = re.split(r'Job Id:\s+', output)
        
        for block in blocks:
            if not block.strip():
                continue
                
            lines = block.splitlines()
            job_id = lines[0].strip()
            job_data = {"id": job_id}
            
            # Use a regex to catch "Attribute = Value"
            # Some values might span multiple lines if they are indented
            last_key = None
            for line in lines[1:]:
                match = re.match(r'^\s+([\w.]+)\s+=\s+(.*)$', line)
                if match:
                    key, value = match.groups()
                    job_data[key] = value
                    last_key = key
                elif last_key and line.startswith("    "):
                    # Continuation line
                    job_data[last_key] += line.strip()
            
            # Filter by user if requested
            if filter_user:
                owner = job_data.get("Job_Owner", "")
                if filter_user not in owner:
                    continue
                    
            jobs.append(job_data)
        return jobs

    @staticmethod
    def get_utilization(job_data: Dict) -> Dict:
        """Extracts utilization metrics from job data."""
        return {
            "cpu_percent": job_data.get("resources_used.cpupercent", "0"),
            "mem_used": job_data.get("resources_used.mem", "0kb"),
            "vmem_used": job_data.get("resources_used.vmem", "0kb"),
            "walltime": job_data.get("resources_used.walltime", "00:00:00")
        }

    @staticmethod
    def submit_job(script_path: str, extra_args: List[str] = None) -> str:
        """Submits a job and returns the Job ID."""
        cmd = ["qsub"]
        if extra_args:
            cmd.extend(extra_args)
        cmd.append(script_path)
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            raise Exception(f"Submission failed: {result.stderr}")

    @staticmethod
    def delete_job(job_id: str):
        """Deletes a job."""
        subprocess.run(["qdel", job_id], check=True)
