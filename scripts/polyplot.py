#!/usr/bin/env python3

"""
@file polyplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-04

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
import sys

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.poly.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "poly":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

    # Postprocessing.
    for j in range(len(lines)):
        if lines[j]:
            if lines[j][0] == "@":
                continue

        lines[j] = lines[j].replace("[", "").replace("]", "")
        lines[j] = lines[j].replace("(", "").replace(")", "")
        lines[j] = lines[j].replace(" ", "")

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

    data: list[str] = line.split(",")
    
    try:

        for j in range(0, len(data), 2):
            x.append(float(data[j]))
            y.append(float(data[j + 1]))

    except ValueError:
        continue

    x.append(x[0])
    y.append(y[0])

    plt.plot(x, y, color = (0, 0, 0), linewidth = 0.5)

plt.show()