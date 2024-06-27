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

plt.figure()

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
    color: tuple[int] = (1, 1, 1) if "--degrees" not in sys.argv else cm.coolwarm((int(data[-1]) - 1) * 200)

    # Plot.
    plt.fill(x, y, facecolor = color, edgecolor = (0, 0, 0), linewidth = 0.5)

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

# Output.
plt.show()