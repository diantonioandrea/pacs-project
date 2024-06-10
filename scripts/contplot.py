#!/usr/bin/env python3

"""
@file contplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-13

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import sys

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.cont.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "cont":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .cont file.")
    sys.exit(-1)

# Points.
x: list[float] = []
y: list[float] = []

# Data.
numerical: list[float] = []
exact: list[float] = []
error: list[float] = []

# Reads values.
for line in lines:
    if line:
        if line[0] == "@":
            continue

    values: list[str] = line.split(",")

    try:
        value_x: float = float(values[0])
        value_y: float = float(values[1])
        value_numerical: float = float(values[2])
        value_exact: float = float(values[3])

        x.append(value_x)
        y.append(value_y)
        numerical.append(value_numerical)
        exact.append(value_exact)
        error.append(abs(value_exact - value_numerical))

    except ValueError:
        continue

# Plot.
fig, axes = plt.subplots(1, 3)
fig.suptitle(lines[0])

data: list[list[float]] = [numerical, exact, error]
titles: list[str] = ["Numerical solution", "Exact solution", "Error"]

contours: list = [None] * 3
bars: list = [None] * 3

for j in range(3):

    # Plot.
    contours[j] = axes[j].tricontourf(x, y, data[j], cmap=cm.coolwarm, levels=100)
    axes[j].set_title(titles[j])

    # Proportions.
    axes[j].set_aspect('equal', adjustable='box')

    # Colorbars.
    bars[j] = fig.colorbar(contours[j], cmap=cm.coolwarm, orientation="horizontal")
    bars[j].set_ticks([min(data[j]), max(data[j])])

    if j == 2:
        bars[j].set_ticklabels([f"{min(data[j]):1.2e}", f"{max(data[j]):1.2e}"])

plt.show()