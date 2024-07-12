#!/usr/bin/env python3

"""
@file polyplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-04

@copyright Copyright (c) 2024
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import sys

# Font.
matplotlib.rcParams.update({'font.size': 18})

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.poly.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "poly":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .poly file.")
    sys.exit(-1)

# Black.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]

# Create a figure and axis
fig, ax = plt.subplots()

# Normalize the colormap
degrees: list[int] = []

if "--degrees" in sys.argv:
    for line in lines:
        if line:
            if line[0] == "@":
                continue

        data: list[str] = line.split(" ")
        
        try:
            degrees.append(int(data[-1]))

        except ValueError:
            continue

if not degrees:
    degrees = [1, 2]

norm = Normalize(vmin=min(degrees), vmax=max(degrees))

for line in lines:
    if line:
        if line[0] == "@":
            continue

    x: list[float] = []
    y: list[float] = []

    data: list[str] = line.split(" ")
    data: list[float] = [float(number) for number in data if number]
    
    try:
        for j in range(0, len(data) if len(data) % 2 == 0 else len(data) - 1, 2):
            x.append(float(data[j]))
            y.append(float(data[j + 1]))

    except ValueError:
        continue

    if not (x and y):
        continue

    # Color.
    color: tuple[int] = [1, 1, 1, 1] if "--degrees" not in sys.argv else list(cm.Blues(norm(int(data[-1]))))
    color[3] = 0.75 # Reduces alpha.

    # Plot.
    ax.fill(x, y, facecolor=color, edgecolor=black, linewidth=0.25)

# Create a ScalarMappable for the colorbar
sm = plt.cm.ScalarMappable(cmap=cm.Blues, norm=norm)
sm.set_array([])

# Add colorbar
if "--degrees" in sys.argv:
    fig.colorbar(sm, ax=ax, ticks=range(min(degrees), max(degrees) + 1), alpha=0.75)

ax.set_aspect('equal', adjustable='box')

# Output.
plt.show()