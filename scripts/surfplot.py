#!/usr/bin/env python3

"""
@file surfplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-13

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import sys

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.surf.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "surf":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .surf file.")
    sys.exit(-1)

x: list[float] = []
y: list[float] = []
numerical: list[float] = []
exact: list[float] = []

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

    except ValueError:
        continue

# Plot.
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_trisurf(x, y, numerical, cmap=cm.coolwarm, linewidth=0, antialiased=True)
plt.show()