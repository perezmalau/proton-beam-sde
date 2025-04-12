import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.ndimage import gaussian_filter
from proton_path import proton_path
import seaborn
import sys

"""Initilise starting configoration"""
x0 = [0, 0, 0]
omega0 = [np.pi / 2, 0]
E0 = 100
dt = 0.05

N = 10000

"""Simulate paths"""

# Radius of the circle
radius = 0.05

# Generate random angles
theta = np.random.uniform(0, 2 * np.pi, N)
phi = np.random.uniform(0, np.pi, N)

# Points on the circle
y = radius * np.sin(phi) * np.cos(theta)
z = radius * np.sin(phi) * np.sin(theta)

# Points on the plane x = 0
x = np.zeros(N)

# Combine the coordinates
x = np.column_stack((x, y, z))

# Loop over points
x_range = (0, 9)
y_range = (-5, 5)
result = proton_path(x[0], omega0, E0, dt)
x_res = []
y_res = []
s_res = []
for i in range(len(result.points)):
    x_res.append(result.points[i][0][0])
    y_res.append(result.points[i][0][1])
    s_res.append(result.points[i][3])
heatmap, xedges, yedges = np.histogram2d(
    x_res, y_res, bins=90, range=[x_range, y_range], weights=s_res, density=False
)

for point in x[1:]:
    result = proton_path(point, omega0, E0, dt)
    x_res = []
    y_res = []
    s_res = []
    for i in range(len(result.points)):
        x_res.append(result.points[i][0][0])
        y_res.append(result.points[i][0][1])
        s_res.append(result.points[i][3])
    heatmap_tmp, xedges, yedges = np.histogram2d(
        x_res, y_res, bins=90, range=[x_range, y_range], weights=s_res, density=False
    )
    heatmap = heatmap + heatmap_tmp

np.divide(heatmap, N)

fig, ax = plt.subplots(figsize=(8, 8))  # Adjust the size as needed
seaborn.heatmap(
    np.transpose(heatmap),
    cmap="viridis",
    norm=mcolors.LogNorm(),
    yticklabels=[round(x, 1) for x in np.linspace(-5, 5, 90)],
    xticklabels=[round(x, 1) for x in np.linspace(0, 9, 90)],
)

plt.savefig("heatmap.pdf")


"""Monte-Carlo for Bragg peak"""

fig, ax = plt.subplots()
plt.plot(np.divide(np.sum(heatmap, axis=1), N))

plt.savefig("bragg-peak.pdf")
