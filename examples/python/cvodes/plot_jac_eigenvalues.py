#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

# Usage: python plot_jac_eigenvalues.py [jac_eigenvalues.txt]
if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = "jac_eigenvalues.txt"

# Read eigenvalue data
with open(fname, "r") as f:
    lines = [line for line in f if not line.startswith("#") and line.strip()]
data = np.loadtxt(lines)
t = data[:, 0]
min_eig = data[:, 1]
max_eig = data[:, 2]
ratio = np.abs(max_eig) / (np.abs(min_eig) + 1e-16)

# Read solution data
sol_file = "cv_vdp_solution.txt"
try:
    with open(sol_file, "r") as f:
        sol_lines = [line for line in f if not line.startswith("#") and line.strip()]
    sol_data = np.loadtxt(sol_lines)
    t_sol = sol_data[:, 0]
    y1 = sol_data[:, 1]
    y2 = sol_data[:, 2]
except Exception as e:
    print(f"Could not plot solution from {sol_file}: {e}")
    t_sol = y1 = y2 = None

# Create combined figure
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Top row: eigenvalues
axs[0, 0].plot(t, max_eig, label="Largest eigenvalue")
axs[0, 0].set_xlabel("Time")
axs[0, 0].set_ylabel("Largest eigenvalue")
axs[0, 0].set_title("Largest Jacobian Eigenvalue vs Time")
axs[0, 0].grid(True)
axs[0, 0].legend()

axs[0, 1].plot(t, ratio, label="|Largest|/|Smallest|")
axs[0, 1].set_xlabel("Time")
axs[0, 1].set_ylabel("Eigenvalue Ratio")
axs[0, 1].set_title("Ratio of Largest/Smallest Eigenvalue vs Time")
axs[0, 1].grid(True)
axs[0, 1].legend()

# Bottom row: solution
if t_sol is not None:
    axs[1, 0].plot(t_sol, y1, label="y1")
    axs[1, 0].set_xlabel("Time")
    axs[1, 0].set_ylabel("y1")
    axs[1, 0].set_title("y1 vs Time")
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    axs[1, 1].plot(t_sol, y2, label="y2", color="orange")
    axs[1, 1].set_xlabel("Time")
    axs[1, 1].set_ylabel("y2")
    axs[1, 1].set_title("y2 vs Time")
    axs[1, 1].grid(True)
    axs[1, 1].legend()
else:
    axs[1, 0].text(0.5, 0.5, "No solution data", ha="center", va="center")
    axs[1, 1].text(0.5, 0.5, "No solution data", ha="center", va="center")

plt.tight_layout()
plt.savefig("jac_eigenvalues_and_solution.png")
print("Combined plot saved as jac_eigenvalues_and_solution.png")
