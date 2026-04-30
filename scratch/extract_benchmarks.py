import os
import json
import re
import argparse

def parse_lammps_log(filepath):
    data = {}
    try:
        if not os.path.exists(filepath):
            return None
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
            for i, line in enumerate(lines):
                # Loop time and basics
                if "Loop time of" in line:
                    loop_match = re.search(r"Loop time of ([\d.]+) on (\d+) procs for (\d+) steps with (\d+) atoms", line)
                    if loop_match:
                        data['loop_time'] = float(loop_match.group(1))
                        data['mpi_procs'] = int(loop_match.group(2))
                        data['steps'] = int(loop_match.group(3))
                        data['atoms'] = int(loop_match.group(4))
                        data['time_10k'] = data['loop_time'] * (10000 / data['steps'])
                
                # CPU use
                if "CPU use with" in line:
                    cpu_match = re.search(r"([\d.]+)% CPU use", line)
                    if cpu_match: data['cpu_use_pct'] = float(cpu_match.group(1))

                # Timing breakdown
                if "Section |" in line and "min time" in line:
                    for j in range(i + 2, i + 15):
                        if j >= len(lines): break
                        table_line = lines[j].strip()
                        if not table_line or "|" not in table_line: break
                        parts = table_line.split('|')
                        if len(parts) < 2: continue
                        label = parts[0].strip()
                        val_str = parts[-1].strip()
                        if val_str:
                            try:
                                val = float(val_str)
                                if label == "Pair": data['pair_pct'] = val
                                elif label == "Bond": data['bond_pct'] = val
                                elif label == "Neigh": data['neigh_pct'] = val
                                elif label == "Comm": data['comm_pct'] = val
                                elif label == "Modify": data['modify_pct'] = val
                                elif label == "Other": data['other_pct'] = val
                            except: pass

                # Other stats
                if "Ave neighs/atom =" in line:
                    neigh_match = re.search(r"Ave neighs/atom = ([\d.]+)", line)
                    if neigh_match: data['ave_neighs'] = float(neigh_match.group(1))
                
                # Memory Usage
                if "Per MPI rank memory allocation" in line:
                    mem_match = re.search(r"max\) = ([\d.]+) \| ([\d.]+) \| ([\d.]+) Mbytes", line)
                    if mem_match:
                        data['mem_max_rank'] = float(mem_match.group(3))

    except Exception as e:
        pass
    return data

def main():
    parser = argparse.ArgumentParser(description="Compare LAMMPS benchmarks across directories")
    parser.add_argument("dirs", nargs="*", default=["dumping_yard"], help="Directories to search for runs")
    args = parser.parse_args()

    results = []
    
    for base_dir in args.dirs:
        if not os.path.isdir(base_dir):
            continue
            
        for root, dirs, files in os.walk(base_dir):
            if "lammps.log" in files and "metadata.json" in files:
                log_path = os.path.join(root, "lammps.log")
                meta_path = os.path.join(root, "metadata.json")
                
                try:
                    with open(meta_path, 'r') as f:
                        meta = json.load(f)
                        n_val = meta.get('N', 'Unknown')
                        if isinstance(n_val, list): n_val = n_val[0]
                except:
                    n_val = "Unknown"
                
                benchmark = parse_lammps_log(log_path)
                # If loop_time is missing, it's a "running" simulation
                if benchmark is None: continue
                
                benchmark['N'] = n_val
                # Use the parent directory name as the simulation tag
                benchmark['sim_tag'] = os.path.basename(os.path.dirname(root))
                benchmark['path'] = root
                results.append(benchmark)
    
    if not results:
        print("No completed benchmarks found.")
        return

    # Sort by N, then MPI, then Simulation Tag
    def sort_key(x):
        try:
            n = int(x.get('N', 0))
        except:
            n = 0
        return (n, x.get('mpi_procs', 0), x.get('sim_tag', ''))
    
    results.sort(key=sort_key)
    
    # Print unified table
    header = "| N | MPI | Simulation | Time 10k (s) | Perf | CPU% | Pair% | Bond% | Modify% | Other% | RAM/Rank |"
    sep = "|---|---|---|---|---|---|---|---|---|---|---|"
    print(header)
    print(sep)
    
    for r in results:
        # Only show rows with loop_time (completed benchmarks)
        if 'loop_time' not in r: continue
        
        perf = f"{r.get('steps', 0)/r.get('loop_time', 1):.1f}"
        cols = [
            str(r.get('N', '?')),
            str(r.get('mpi_procs', '?')),
            r.get('sim_tag', 'Unknown'),
            f"{r.get('time_10k', 0):.2f}",
            perf,
            f"{r.get('cpu_use_pct', 0):.1f}",
            f"{r.get('pair_pct', 0):.2f}",
            f"{r.get('bond_pct', 0):.2f}",
            f"{r.get('modify_pct', 0):.2f}",
            f"{r.get('other_pct', 0):.2f}",
            f"{r.get('mem_max_rank', 0):.1f}"
        ]
        print("| " + " | ".join(cols) + " |")

if __name__ == "__main__":
    main()
