import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1. Setup synthetic data based on your images
# (In a real scenario, you would load your .npy or .csv data here)
x = np.linspace(-5e6, 5e6, 200)
y = np.linspace(-2.5e6, 2.5e6, 100)
X, Y = np.meshgrid(x, y)

# Simulating the Topography anomaly (the circular peak)
topo = 20000 * np.exp(-((X + 1.5e6)**2 + Y**2) / (2 * 2.5e5**2))

# Simulating the Moho Depth (baseline -30km with a deep root under the mountain)
moho = -30000 - 50000 * np.exp(-((X + 1.5e6)**2 + Y**2) / (2 * 4e5**2))

# 2. Create the 3D Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot Topography
surf1 = ax.plot_surface(X, Y, topo, cmap='gist_earth', edgecolor='none', alpha=0.8)

# Plot Moho Undulation
surf2 = ax.plot_surface(X, Y, moho, cmap='RdBu_r', edgecolor='none', alpha=0.6)

# 3. Formatting
ax.set_title('3D Visualization: Topography vs. Moho Depth', fontsize=15)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Elevation/Depth (m)')

# Adjusting the Z-axis to show the depth clearly
ax.set_zlim(-80000, 30000)

plt.colorbar(surf1, shrink=0.5, aspect=10, label='Topography (m)')
plt.show()

