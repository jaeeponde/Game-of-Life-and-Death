# quarantine_seir_simulation.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import random
import os

# -------------------------------
# Parameters
# -------------------------------
N = 150
p = 0.5  # Infection probability
q = 0.2  # Recovery probability
steps = 200
num_runs = 10
movement_mode = "random"  # 'none', 'random', 'restricted'

# -------------------------------
# States
# -------------------------------
SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

# Color map
cmap = colors.ListedColormap(['green', 'red', 'black'])
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

# -------------------------------
# Helper Functions
# -------------------------------
def init_grid(N, num_initial_infected):
    grid = np.zeros((N, N), dtype=int)
    positions = set()
    while len(positions) < num_initial_infected:
        i = random.randint(0, N - 1)
        j = random.randint(0, N - 1)
        positions.add((i, j))
    for i, j in positions:
        grid[i, j] = INFECTED
    return grid

def move_individuals(grid, mode):
    if mode == "none":
        return grid.copy()

    new_grid = grid.copy()
    indices = list(np.ndindex(N, N))
    random.shuffle(indices)

    for i, j in indices:
        if grid[i, j] == SUSCEPTIBLE or grid[i, j] == INFECTED:
            if mode == "random":
                di, dj = random.choice([(0,1),(0,-1),(1,0),(-1,0)])
                ni, nj = i + di, j + dj
                if 0 <= ni < N and 0 <= nj < N and new_grid[ni, nj] == SUSCEPTIBLE:
                    new_grid[ni, nj], new_grid[i, j] = new_grid[i, j], new_grid[ni, nj]
            elif mode == "restricted":
                if random.random() < 0.3:
                    di, dj = random.choice([(0,1),(0,-1),(1,0),(-1,0)])
                    ni, nj = i + di, j + dj
                    if 0 <= ni < N and 0 <= nj < N and new_grid[ni, nj] == SUSCEPTIBLE:
                        new_grid[ni, nj], new_grid[i, j] = new_grid[i, j], new_grid[ni, nj]
    return new_grid

def simulate_once(movement_mode, save_animation=False, run_id=0, mode_label=""):
    grid = init_grid(N, 10)
    S_series, I_series, R_series = [], [], []
    frames = []
    heatmaps = []

    fig, ax = plt.subplots()
    mat = ax.matshow(grid, cmap=cmap, norm=norm)

    for t in range(steps):
        grid = move_individuals(grid, movement_mode)
        new_grid = grid.copy()
        for i in range(N):
            for j in range(N):
                if grid[i, j] == SUSCEPTIBLE:
                    for ni, nj in [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]:
                        if 0 <= ni < N and 0 <= nj < N and grid[ni, nj] == INFECTED:
                            if random.random() < p:
                                new_grid[i, j] = INFECTED
                                break
                elif grid[i, j] == INFECTED:
                    if random.random() < q:
                        new_grid[i, j] = RECOVERED

        S_series.append(np.sum(grid == SUSCEPTIBLE))
        I_series.append(np.sum(grid == INFECTED))
        R_series.append(np.sum(grid == RECOVERED))

        heatmaps.append(grid.copy())
        grid = new_grid

        if save_animation:
            frame = ax.matshow(grid.copy(), cmap=cmap, norm=norm, animated=True)
            frames.append([frame])

    if save_animation:
        ani = animation.ArtistAnimation(fig, frames, interval=100, blit=True, repeat=False)
        filename = f"simulations3/simulation_{mode_label}_run{run_id}.mp4"
        ani.save(filename)
        plt.close(fig)

    return np.array(S_series), np.array(I_series), np.array(R_series), heatmaps

def analyze_runs(results):
    I_all = np.array([r[1] for r in results])
    mean_I = I_all.mean(axis=0)
    var_I = I_all.var(axis=0)

    peak_intensity = mean_I.max()
    peak_time = np.argmax(mean_I)
    above_thresh = mean_I > 0.75 * peak_intensity
    peak_duration = np.sum(above_thresh)

    threshold = 0.01 * (N*N)
    for t in range(len(mean_I)):
        if np.all(mean_I[t:] < threshold):
            equilibrium_time = t
            break
    else:
        equilibrium_time = steps

    return mean_I, var_I, peak_intensity, peak_time, peak_duration, equilibrium_time

def plot_results(mean_I, var_I, label):
    x = np.arange(len(mean_I))
    plt.figure()
    plt.plot(x, mean_I, color='purple', label='Mean Infected')
    plt.fill_between(x, mean_I - np.sqrt(var_I), mean_I + np.sqrt(var_I), alpha=0.3, color='purple')
    plt.title(f"Infection Curve - {label}")
    plt.xlabel("Time Step")
    plt.ylabel("Infected")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"infected_curve_{label}.png")
    plt.close()

def comparative_plot(metrics_dict):
    modes = list(metrics_dict.keys())
    for metric in ['peak', 'duration', 'equilibrium']:
        values = [metrics_dict[mode][metric] for mode in modes]
        plt.figure()
        bars = plt.bar(modes, values, color=plt.cm.viridis(np.linspace(0.2, 0.8, len(modes))))
        plt.title(f"{metric.capitalize()} Comparison by Movement Mode")
        plt.xlabel("Mode")
        plt.ylabel(metric.capitalize())
        plt.savefig(f"comparisonplots/{metric}_comparison.png")
        plt.close()

def save_heatmap(grid_series, label):
    snapshots = [grid_series[10], grid_series[len(grid_series)//2], grid_series[-1]]
    titles = ['Early', 'Mid', 'Late']
    for snap, t in zip(snapshots, titles):
        plt.figure()
        plt.imshow(snap, cmap='viridis')
        plt.title(f"{label} - {t} Phase")
        plt.colorbar()
        plt.savefig(f"heatmaps/heatmap_{label.lower()}_{t.lower()}.png")
        plt.close()

# -------------------------------
# Run Simulations
# -------------------------------
all_modes = ["none", "random", "restricted"]
metrics = {}
for mode in all_modes:
    print(f"Running simulations for mode: {mode}")
    run_results = [simulate_once(mode, save_animation=(run == 0), run_id=run, mode_label=mode) for run in range(num_runs)]
    series_only = [(s, i, r) for s, i, r, _ in run_results]
    mean_I, var_I, peak_intensity, peak_time, peak_duration, equilibrium_time = analyze_runs(series_only)
    metrics[mode] = {"peak": peak_intensity, "duration": peak_duration, "equilibrium": equilibrium_time}
    print(f"{mode.title()} Mode: Peak = {peak_intensity:.2f} at t={peak_time}, Duration = {peak_duration}, Equilibrium = {equilibrium_time}")
    plot_results(mean_I, var_I, mode)
    save_heatmap(run_results[0][3], mode)

comparative_plot(metrics)
print("All simulations completed.") 
