import ovito
from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier
from pathlib import Path

print("Testing OVITO pipeline...")
# Create a dummy dataset or just verify the modifiers
mod = ConstructSurfaceModifier(radius=0.005)
print("Modifier created successfully!")
