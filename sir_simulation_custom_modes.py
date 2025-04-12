import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import random

# Grid size
N = 150

# Infection and recovery probabilities
p = 0.5
q = 0.2

# Simulation steps
steps = 1000

# States
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

# Colors: green = S, red = I, black = R
cmap = colors.ListedColormap(['green', 'red', 'black'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Initialize grid
grid = np.zeros((N, N), dtype=int)

# -------------------------------
# Number of initially infected
# -------------------------------
num_initial_infected = 15  # 👈 change this value

# -------------------------------
# Placement mode: "random", "even", "far"
# -------------------------------
placement_mode = "far"

if placement_mode == "random":
    positions = set()
    while len(positions) < num_initial_infected:
        i = random.randint(0, N - 1)
        j = random.randint(0, N - 1)
        positions.add((i, j))
    initial_infected_positions = list(positions)

elif placement_mode == "even":
    initial_infected_positions = []
    step = N // int(np.sqrt(num_initial_infected) + 1)
    for i in range(step, N - step, step):
        for j in range(step, N - step, step):
            if len(initial_infected_positions) < num_initial_infected:
                initial_infected_positions.append((i, j))

elif placement_mode == "far":
    # Spread out: corners and midpoints
    initial_infected_positions = [
        (1, 1),
        (1, N-2),
        (N-2, 1),
        (N-2, N-2),
        (N//2, N//2)
    ][:num_initial_infected]

else:
    raise ValueError("Invalid placement_mode. Choose from 'random', 'even', or 'far'.")

# Infect selected positions
for i, j in initial_infected_positions:
    if 0 <= i < N and 0 <= j < N:
        grid[i, j] = INFECTED

# Track SIR counts over time
S_counts = []
I_counts = []
R_counts = []


# -------------------------------
# Update function
# -------------------------------
def update(frame):
    global grid
    new_grid = grid.copy()

    s_count = np.sum(grid == SUSCEPTIBLE)
    i_count = np.sum(grid == INFECTED)
    r_count = np.sum(grid == RECOVERED)

    S_counts.append(s_count)
    I_counts.append(i_count)
    R_counts.append(r_count)


    for i in range(N):
        for j in range(N):
            if grid[i, j] == SUSCEPTIBLE:
                neighbors = [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
                for ni, nj in neighbors:
                    if 0 <= ni < N and 0 <= nj < N:
                        if grid[ni, nj] == INFECTED and np.random.rand() < p:
                            new_grid[i, j] = INFECTED
                            break
            elif grid[i, j] == INFECTED:
                if np.random.rand() < q:
                    new_grid[i, j] = RECOVERED

    grid[:] = new_grid
    mat.set_data(grid)
    ax.set_title(f"SIR Model at t={frame}")
    return [mat]

# -------------------------------
# Plot and animate
# -------------------------------
fig, ax = plt.subplots()
mat = ax.matshow(grid, cmap=cmap, norm=norm)
ax.set_title("SIR Model at t=0")

ani = animation.FuncAnimation(fig, update, frames=steps, interval=200, blit=True)
plt.show()


# Plot the final SIR graph
total = N * N
plt.figure()
plt.plot(np.array(S_counts) / total, label='Susceptible', color='green')
plt.plot(np.array(I_counts) / total, label='Infected', color='red')
plt.plot(np.array(R_counts) / total, label='Recovered', color='black')
plt.xlabel("Time Step")
plt.ylabel("Proportion of Population")
plt.title("SIR Model Dynamics")
plt.legend()
plt.grid(True)
plt.show()
