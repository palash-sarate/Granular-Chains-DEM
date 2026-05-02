import sys
import os
sys.path.append("Pulse")
from pulse_core import PBSManager
print(f"PBSManager has load_lineage: {hasattr(PBSManager, 'load_lineage')}")
if hasattr(PBSManager, 'load_lineage'):
    print("Success")
else:
    print("Failure: load_lineage not found")
    print(f"Methods in PBSManager: {[m for m in dir(PBSManager) if not m.startswith('__')]}")
