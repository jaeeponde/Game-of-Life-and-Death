import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random
import os

# Parameters
GRID_SIZE = 150
TOTAL_POPULATION = GRID_SIZE * GRID_SIZE
BETA = 0.5
GAMMA = 0.2
R0 = BETA / GAMMA
HERD_IMMUNITY_THRESHOLD = 1 - (1 / R0)
vaccination_levels = [0.0, HERD_IMMUNITY_THRESHOLD - 0.2, HERD_IMMUNITY_THRESHOLD, HERD_IMMUNITY_THRESHOLD + 0.1]
timesteps = 100

# SIR color map: 0 = S (green), 1 = I (red), 2 = R (black), 3 = V (gray)
cmap = mcolors.ListedColormap(["green", "red", "black", "gray"])

def initialize_grid(vaccination_rate, initial_infected=4):
    grid = np.zeros((GRID_SIZE, GRID_SIZE))  # 0: S, 1: I, 2: R, 3: V
    num_vaccinated = int(TOTAL_POPULATION * vaccination_rate)
    vaccinated_positions = random.sample(range(TOTAL_POPULATION), num_vaccinated)

    for pos in vaccinated_positions:
        x, y = divmod(pos, GRID_SIZE)
        grid[x, y] = 3  # Vaccinated

    # Infect a given number of susceptible individuals
    susceptible_positions = [i for i in range(TOTAL_POPULATION) if grid[i // GRID_SIZE, i % GRID_SIZE] == 0]
    if len(susceptible_positions) >= initial_infected:
        infected_positions = random.sample(susceptible_positions, initial_infected)
        for pos in infected_positions:
            x, y = divmod(pos, GRID_SIZE)
            grid[x, y] = 1  # Infected
    else:
        print("Not enough susceptibles to infect the initial number.")

    return grid

def simulate(grid, vaccination_rate, tag):
    S, I, R = [], [], []

    folder_path = f"sir_grids/{tag}"
    os.makedirs(folder_path, exist_ok=True)

    for t in range(timesteps):
        new_grid = grid.copy()

        # Count current states
        unique, counts = np.unique(grid, return_counts=True)
        count_dict = dict(zip(unique, counts))
        S.append(count_dict.get(0, 0))
        I.append(count_dict.get(1, 0))
        R.append(count_dict.get(2, 0))

        # Save snapshot
        if t % (timesteps // 8) == 0 or t == timesteps - 1:
            plt.figure(figsize=(6, 6))
            plt.imshow(grid, cmap=cmap, vmin=0, vmax=3)
            plt.title(f"Vaccination Rate = {vaccination_rate:.2f}\nTimestep {t}")
            plt.axis('off')
            plt.savefig(f"{folder_path}/grid_t{t}.png")
            plt.close()

        # Disease spread
        for x in range(GRID_SIZE):
            for y in range(GRID_SIZE):
                if grid[x, y] == 1:  # Infected
                    if random.random() < GAMMA:
                        new_grid[x, y] = 2  # Recover
                    else:
                        for dx in [-1, 0, 1]:
                            for dy in [-1, 0, 1]:
                                if dx == 0 and dy == 0:
                                    continue
                                nx, ny = x + dx, y + dy
                                if 0 <= nx < GRID_SIZE and 0 <= ny < GRID_SIZE:
                                    if grid[nx, ny] == 0 and random.random() < BETA:
                                        new_grid[nx, ny] = 1
        grid = new_grid

    return S, I, R

# Run and collect infection curves
all_curves = {}

for v_rate in vaccination_levels:
    tag = f"v{int(v_rate * 100)}"
    print(f"Running simulation for vaccination rate: {v_rate:.2f}")
    grid = initialize_grid(v_rate, initial_infected=4)
    S, I, R = simulate(grid, v_rate, tag)
    all_curves[tag] = (S, I, R)

# Plot comparison of infected curves
plt.figure(figsize=(10, 6))
for tag, (S, I, R) in all_curves.items():
    plt.plot(I, label=f"Infected ({tag[1:]}%)")

plt.title("Infection Curve Comparison")
plt.xlabel("Time")
plt.ylabel("Number of Infected Individuals")
plt.legend()
plt.grid(False)
plt.set_cmap("viridis")
plt.tight_layout()
plt.savefig("sir_grids/infection_curve_comparison.png")
plt.show()
