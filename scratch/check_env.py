import shutil
import subprocess
import os

def check_env():
    print(f"PATH: {os.environ.get('PATH')}")
    print(f"lmp path: {shutil.which('lmp')}")
    print(f"nvidia-smi path: {shutil.which('nvidia-smi')}")
    
    try:
        res = subprocess.run(['lmp', '-h'], capture_output=True, text=True, timeout=5)
        print("LAMMPS -h output (Installed packages):")
        capture = False
        for line in res.stdout.splitlines():
            if "Installed packages:" in line:
                capture = True
                continue
            if capture:
                if not line.strip() or line.startswith("List of"):
                    break
                print(line.strip())
    except Exception as e:
        print(f"Error running lmp -h: {e}")

if __name__ == "__main__":
    check_env()
