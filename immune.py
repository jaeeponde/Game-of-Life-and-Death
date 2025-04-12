import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import random

# -------------------------------
# Parameters
# -------------------------------
N = 150              # Grid size (NxN)
p = 0.5             # Infection probability
q = 0.2             # Recovery probability
steps = 100        # Number of time steps
vaccination_rate = 0.1  # 👈 % of population vaccinated

# -------------------------------
# State definitions
# -------------------------------
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2
IMMUNE = 3  # Vaccinated

# Color map: green = S, red = I, black = R, blue = vaccinated
cmap = colors.ListedColormap(['green', 'red', 'black', 'blue'])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# -------------------------------
# Initialize grid
# -------------------------------
grid = np.zeros((N, N), dtype=int)

# Place infected individuals evenly
num_initial_infected = 5
infected_positions = []
step = N // int(np.sqrt(num_initial_infected) + 1)
for i in range(step, N - step, step):
    for j in range(step, N - step, step):
        if len(infected_positions) < num_initial_infected:
            infected_positions.append((i, j))
for i, j in infected_positions:
    grid[i, j] = INFECTED

# Place immune (vaccinated) individuals randomly
num_immune_cells = int(vaccination_rate * N * N)
immune_positions = set()
while len(immune_positions) < num_immune_cells:
    i = random.randint(0, N - 1)
    j = random.randint(0, N - 1)
    if grid[i, j] == SUSCEPTIBLE:
        grid[i, j] = IMMUNE
        immune_positions.add((i, j))

# -------------------------------
# Track SIR dynamics
# -------------------------------
S_counts, I_counts, R_counts = [], [], []

def update(frame):
    global grid
    new_grid = grid.copy()

    # Count current states
    S_counts.append(np.sum(grid == SUSCEPTIBLE))
    I_counts.append(np.sum(grid == INFECTED))
    R_counts.append(np.sum(grid == RECOVERED))

    # Update rules
    for i in range(N):
        for j in range(N):
            state = grid[i, j]
            if state == SUSCEPTIBLE:
                neighbors = [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
                for ni, nj in neighbors:
                    if 0 <= ni < N and 0 <= nj < N:
                        if grid[ni, nj] == INFECTED and np.random.rand() < p:
                            new_grid[i, j] = INFECTED
                            break
            elif state == INFECTED:
                if np.random.rand() < q:
                    new_grid[i, j] = RECOVERED
            # Recovered and Immune remain unchanged

    grid[:] = new_grid
    mat.set_data(grid)
    ax.set_title(f"Vaccination Simulation at t={frame}")
    return [mat]

# -------------------------------
# Animate the simulation
# -------------------------------
fig, ax = plt.subplots()
mat = ax.matshow(grid, cmap=cmap, norm=norm)
ani = animation.FuncAnimation(fig, update, frames=steps, interval=200, blit=True, repeat=False)
plt.show()

# -------------------------------
# Plot SIR graph after animation
# -------------------------------
total = N * N
plt.figure()
plt.plot(np.array(S_counts) / total, label='Susceptible', color='green')
plt.plot(np.array(I_counts) / total, label='Infected', color='red')
plt.plot(np.array(R_counts) / total, label='Recovered', color='black')
plt.xlabel("Time Step")
plt.ylabel("Proportion of Population")
plt.title(f"SIR Dynamics with {int(vaccination_rate*100)}% Vaccinated")
plt.legend()
plt.grid(True)
plt.show()
