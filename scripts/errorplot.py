#!/usr/bin/env python3

"""
@file errorplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-19

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys

elements: list[int] = []
sizes: list[float] = []
errors: list[float] = []

if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

try:

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

# Simple scraping.
for line in lines:
    try:
        data: list[str] = line.split(" ")

        if "Elements" in line:
            elements.append(int(data[-1]))

        elif "Size" in line:
            sizes.append(float(data[-1]))

        elif "L2 Error" in line:
            errors.append(float(data[-1]))

    except ValueError:
        continue

# Comparisons.
elements_linear = [1.0 / (value / elements[-1]) * errors[-1] for value in elements]
elements_quadratic_root = [1.0 / ((value / elements[-1]) ** 0.5) * errors[-1] for value in elements]

sizes_linear = [value * (errors[0] / sizes[0]) for value in sizes]
sizes_quadratic = [(value / sizes[0]) ** 2 * errors[0] for value in sizes]

# Ticks.
elements_ticks = [elements[0], elements[-1]]
elements_labels = [str(tick) for tick in elements_ticks]

sizes_ticks = [sizes[0], sizes[-1]]
sizes_labels = [f"{tick:1.2f}" for tick in sizes_ticks]

if "--single" in sys.argv: # Elements only.
    # Plot.
    fig, axes = plt.subplots()
    fig.suptitle(lines[0])

    # Elements.
    axes.plot(elements, elements_linear, linewidth=1.0, linestyle=":", color="red") # Linear comparison.
    axes.plot(elements, elements_quadratic_root, linewidth=1.0, linestyle=":", color="red") # Quadratic comparison.
    axes.plot(elements, errors, color="black") # Errors.

    # Loglog scale.
    axes.set_xscale("log")
    axes.set_yscale("log")

    # Ticks.
    axes.xaxis.set_minor_formatter(NullFormatter())
    axes.set_xticks(elements_ticks, labels=elements_labels)

    plt.show()
    sys.exit(0)

# Plot.
fig, axes = plt.subplots(1, 2)
fig.suptitle(lines[0])

# Elements.
axes[0].plot(elements, elements_linear, linewidth=1.0, linestyle=":", color="red") # Linear comparison.
axes[0].plot(elements, elements_quadratic_root, linewidth=1.0, linestyle=":", color="red") # Quadratic comparison.
axes[0].plot(elements, errors, color="black") # Errors.

# Sizes.
axes[1].plot(sizes, sizes_linear, linewidth=1.0, linestyle=":", color="red") # Linear comparison.
axes[1].plot(sizes, sizes_quadratic, linewidth=1.0, linestyle=":", color="red") # Quadratic comparison.
axes[1].plot(sizes, errors, color="black") # Errors.

# Loglog scale.
axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[1].set_xscale("log")
axes[1].set_yscale("log")

# Ticks.
axes[0].xaxis.set_minor_formatter(NullFormatter())
axes[1].xaxis.set_minor_formatter(NullFormatter())
axes[0].yaxis.set_minor_formatter(NullFormatter())
axes[1].yaxis.set_minor_formatter(NullFormatter())
axes[0].set_xticks(elements_ticks, labels=elements_labels)
axes[1].set_xticks(sizes_ticks, labels=sizes_labels)

plt.show()