import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors

# Parameters
N = 150  # grid size (NxN)
p = 0.5  # probability of infection
q = 0.2  # probability of recovery
steps = 100

# State definitions
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

# Set up color map
cmap = colors.ListedColormap(['green', 'red', 'black'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Initialize grid
grid = np.zeros((N, N), dtype=int)
grid[N//2, N//2] = INFECTED  # start with one infected cell in the center

# Function to update the grid
def update(frame):
    global grid
    new_grid = grid.copy()
    for i in range(N):
        for j in range(N):
            if grid[i, j] == SUSCEPTIBLE:
                # Check 4 neighbors
                neighbors = [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
                for ni, nj in neighbors:
                    if 0 <= ni < N and 0 <= nj < N:
                        if grid[ni, nj] == INFECTED and np.random.rand() < p:
                            new_grid[i, j] = INFECTED
                            break
            elif grid[i, j] == INFECTED:
                if np.random.rand() < q:
                    new_grid[i, j] = RECOVERED
    grid = new_grid
    mat.set_data(grid)
    return [mat]

# Set up the plot
fig, ax = plt.subplots()
mat = ax.matshow(grid, cmap=cmap, norm=norm)
plt.title("SIR Model Simulation")

ani = animation.FuncAnimation(fig, update, frames=steps, interval=200, blit=True)
plt.show()
