#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

# Usage: python plot_jac_eigenvalues.py [jac_eigenvalues.txt]
if len(sys.argv) > 1:
    fname = sys.argv[1]
else:
    fname = "jac_eigenvalues.txt"


# Read stiff_test data
with open(fname, "r") as f:
    lines = [line for line in f if not line.startswith("#") and line.strip()]
data = np.loadtxt(lines)
t = data[:, 0]
stiff_test = data[:, 1]


# Plot stiff_test vs time
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(t, stiff_test, label="stiff_test")
plt.xlabel("Time")
plt.ylabel("stiff_test")
plt.title("Stiff Test Value vs Time")
plt.grid(True)
plt.legend()

# Plot solution (y1, y2) vs time
sol_file = "ark_vdp_solution.txt"
try:
    with open(sol_file, "r") as f:
        sol_lines = [line for line in f if not line.startswith("#") and line.strip()]
    sol_data = np.loadtxt(sol_lines)
    t_sol = sol_data[:, 0]
    y1 = sol_data[:, 1]
    y2 = sol_data[:, 2]
    plt.subplot(1, 2, 2)
    plt.plot(t_sol, y1, label="y1")
    plt.plot(t_sol, y2, label="y2")
    plt.xlabel("Time")
    plt.ylabel("Solution")
    plt.title("Van der Pol Solution vs Time")
    plt.grid(True)
    plt.legend()
except Exception as e:
    plt.subplot(1, 2, 2)
    plt.text(0.5, 0.5, "No solution data", ha="center", va="center")

plt.tight_layout()
plt.savefig("stiff_test_and_solution.png")
print("Plot saved as stiff_test_and_solution.png")
