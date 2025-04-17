import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import colors
import random

# Grid size and parameters
N = 150
p = 0.3
q = 0.3
steps = 500
snapshots = list(np.linspace(0, steps, 8, dtype=int))


# States
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

# Color mapping
cmap = colors.ListedColormap(['green', 'red', 'black'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Create directory for screenshots
output_folder = "sir_screenshots"
os.makedirs(output_folder, exist_ok=True)

def initialize_grid(mode="random", num_initial_infected=4):
    grid = np.zeros((N, N), dtype=int)
    if mode == "random":
        positions = set()
        while len(positions) < num_initial_infected:
            i = random.randint(0, N - 1)
            j = random.randint(0, N - 1)
            positions.add((i, j))
        initial_infected_positions = list(positions)
    elif mode == "even":
        initial_infected_positions = []
        sqrt_n = int(np.ceil(np.sqrt(num_initial_infected)))
        xs = np.linspace(10, N - 10, sqrt_n, dtype=int)
        ys = np.linspace(10, N - 10, sqrt_n, dtype=int)
        for x in xs:
            for y in ys:
                if len(initial_infected_positions) < num_initial_infected:
                    initial_infected_positions.append((x, y))
    elif mode == "far":
        initial_infected_positions = [
            (1, 1), (1, N - 2), (N - 2, 1), (N - 2, N - 2),
            (N // 2, N // 2), (N // 4, N // 4), (N // 4, 3 * N // 4),
            (3 * N // 4, N // 4), (3 * N // 4, 3 * N // 4)
        ][:num_initial_infected]
    elif mode == "clustered":
        cx, cy = N // 2, N // 2
        initial_infected_positions = []
        while len(initial_infected_positions) < num_initial_infected:
            dx, dy = np.random.randint(-10, 10), np.random.randint(-10, 10)
            i, j = cx + dx, cy + dy
            if 0 <= i < N and 0 <= j < N:
                initial_infected_positions.append((i, j))
    else:
        raise ValueError("Invalid placement mode")

    for i, j in initial_infected_positions:
        grid[i, j] = INFECTED
    return grid

def simulate_and_save(mode):
    grid = initialize_grid(mode)
    for step in range(steps + 1):
        new_grid = grid.copy()
        for i in range(N):
            for j in range(N):
                if grid[i, j] == SUSCEPTIBLE:
                    neighbors = [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]
                    for ni, nj in neighbors:
                        if 0 <= ni < N and 0 <= nj < N:
                            if grid[ni, nj] == INFECTED and random.random() < p:
                                new_grid[i, j] = INFECTED
                                break
                elif grid[i, j] == INFECTED:
                    if random.random() < q:
                        new_grid[i, j] = RECOVERED
        grid = new_grid

        if step in snapshots:
            index = snapshots.index(step) + 1
            plt.figure(figsize=(6, 6))
            plt.imshow(grid, cmap=cmap, norm=norm)
            plt.title(f"{mode.title()} Placement - Step {step}")
            plt.axis('off')
            filename = f"{index}_{mode}.png"
            plt.savefig(os.path.join(output_folder, filename), bbox_inches='tight')
            plt.close()

# Run for all placements
for mode in ["random", "even", "far", "clustered"]:
    simulate_and_save(mode)

print("✅ Screenshots saved in 'sir_screenshots' folder.")
