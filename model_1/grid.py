import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import random

# Parameters
N = 150
num_initial_infected = 4
placement_mode = "clustered"  # Change to "random", "even", or "far"

# States
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

# Color map
cmap = colors.ListedColormap(['green', 'red', 'black'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Initialize grid
grid = np.zeros((N, N), dtype=int)

# Place infected individuals
if placement_mode == "random":
    positions = set()
    while len(positions) < num_initial_infected:
        i = random.randint(0, N - 1)
        j = random.randint(0, N - 1)
        positions.add((i, j))
    initial_infected_positions = list(positions)

elif placement_mode == "even":
    initial_infected_positions = []
    step = N // (int(np.sqrt(num_initial_infected)) + 1)
    for i in range(step, N - step, step):
        for j in range(step, N - step, step):
            if len(initial_infected_positions) < num_initial_infected:
                initial_infected_positions.append((i, j))

elif placement_mode == "far":
    initial_infected_positions = [
        (1, 1),
        (1, N-2),
        (N-2, 1),
        (N-2, N-2)
    ][:num_initial_infected]

elif placement_mode == "clustered":
        # Cluster the infected individuals in a tight group at the center
    cluster_radius = 2  # Size of the cluster (2x2 block in the center)
    center_i, center_j = N // 2, N // 2
    initial_positions = []

    for di in range(-cluster_radius, cluster_radius + 1):
        for dj in range(-cluster_radius, cluster_radius + 1):
            if len(initial_positions) < num_initial_infected:
                i = center_i + di
                j = center_j + dj
                if 0 <= i < N and 0 <= j < N:
                    initial_positions.append((i, j))

else:
    raise ValueError("Invalid placement_mode. Choose from 'random', 'even', or 'far'.")

# Infect selected positions
for i, j in initial_positions:
    grid[i, j] = INFECTED

# Plot initial placement
plt.figure(figsize=(6, 6))
plt.matshow(grid, cmap=cmap, norm=norm, fignum=1)
plt.title(f"Initial Infection Placement: {placement_mode.capitalize()}")
plt.axis('off')
plt.show()
