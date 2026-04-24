import typer
from rich.console import Console
from rich.table import Table
from rich import print as rprint
from .pulse_core import PBSManager
import os
import subprocess

app = typer.Typer(help="Pulse: Smart PBS Job Manager")
console = Console()

@app.command()
def status(user: str = os.environ.get("USER", "guest")):
    """Show the status of jobs in a beautiful table."""
    jobs = PBSManager.get_jobs(user=user)
    
    if not jobs:
        rprint("[yellow]No active jobs found.[/yellow]")
        return

    table = Table(title=f"Active Jobs for {user}")
    table.add_column("Job ID", style="cyan")
    table.add_column("Name", style="magenta")
    table.add_column("State", style="bold")
    table.add_column("Queue", style="green")
    table.add_column("Time", style="blue")
    table.add_column("Comment", style="white")

    for job in jobs:
        state = job.get("job_state", "?")
        state_style = "bold green" if state == "R" else "bold yellow" if state == "Q" else "bold red" if state == "H" else "white"
        
        table.add_row(
            job.get("id", "N/A"),
            job.get("Job_Name", "N/A"),
            f"[{state_style}]{state}[/{state_style}]",
            job.get("queue", "N/A"),
            job.get("resources_used.walltime", "00:00:00"),
            job.get("comment", "")[:50] # Truncate long comments
        )

    console.print(table)

@app.command()
def submit(script: str):
    """Submit a PBS script with pre-flight checks."""
    if not os.path.exists(script):
        rprint(f"[red]Error: Script {script} not found.[/red]")
        return
        
    try:
        job_id = PBSManager.submit_job(script)
        rprint(f"[green]Job submitted successfully! ID: {job_id}[/green]")
    except Exception as e:
        rprint(f"[red]{e}[/red]")

@app.command()
def delete(job_id: str):
    """Delete a specific job."""
    try:
        PBSManager.delete_job(job_id)
        rprint(f"[green]Job {job_id} deleted.[/green]")
    except Exception as e:
        rprint(f"[red]Error deleting job: {e}[/red]")

import sys

@app.command()
def dashboard():
    """Launch the Pulse Web Dashboard."""
    rprint("[bold blue]Launching Pulse Dashboard...[/bold blue]")
    rprint("[yellow]Note: You may need to use SSH port forwarding (e.g., -L 8501:localhost:8501) to view it locally.[/yellow]")
    
    # Launch streamlit as a subprocess using the current python interpreter
    dashboard_path = os.path.join(os.path.dirname(__file__), "dashboard.py")
    
    try:
        process = subprocess.Popen([sys.executable, "-m", "streamlit", "run", dashboard_path])
        process.wait()
    except KeyboardInterrupt:
        rprint("\n[yellow]Stopping Pulse Dashboard...[/yellow]")
        process.terminate()
    except Exception as e:
        rprint(f"[red]Error: {e}[/red]")
    finally:
        # Ensure the process is truly dead
        if 'process' in locals() and process.poll() is None:
            process.kill()

if __name__ == "__main__":
    app()
