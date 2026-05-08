import vedo
import os

# No allow_interaction in this version

plt = vedo.Plotter(offscreen=True, size=(800, 600))
s = vedo.Sphere().c('red')
plt.add(s)

plt.axes = 1
plt.add_scale_indicator(s=0.1)
plt.background('white')

plt.render()
plt.screenshot('test_viz.png')
plt.close()

if os.path.exists('test_viz.png'):
    print("SUCCESS: test_viz.png created")
    print(f"Size: {os.path.getsize('test_viz.png')} bytes")
else:
    print("FAILED: test_viz.png not created")
