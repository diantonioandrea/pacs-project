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
l2_error: list[float] = []
order: int = 0

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
            l2_error.append(float(data[-1]))

        elif "L2 Error" in line:
            l2_error.append(float(data[-1]))

        elif "Degree" in line:
            order = int(data[-1])

    except ValueError:
        continue

# Comparisons.
elements_comparison: list[list[float]] = [None] * 3
sizes_comparison: list[list[float]] = [None] * 3

elements_comparison[0] = [(value / elements[-1]) ** -(order - 1) * l2_error[-1] for value in elements]
sizes_comparison[0] = [(value / sizes[0]) ** (order - 1) * l2_error[0] for value in sizes]

elements_comparison[1] = [(value / elements[-1]) ** -order * l2_error[-1] for value in elements]
sizes_comparison[1] = [(value / sizes[0]) ** order * l2_error[0] for value in sizes]

elements_comparison[2] = [(value / elements[-1]) ** -(order + 1) * l2_error[-1] for value in elements]
sizes_comparison[2] = [(value / sizes[0]) ** (order + 1) * l2_error[0] for value in sizes]

# Ticks.
elements_ticks = [elements[0], elements[-1]]
elements_labels = [str(tick) for tick in elements_ticks]

sizes_ticks = [sizes[0], sizes[-1]]
sizes_labels = [f"{tick:1.2f}" for tick in sizes_ticks]

errors_ticks = [l2_error[0], l2_error[-1]]
errors_labels = [f"{tick:1.2e}" for tick in errors_ticks]

# Title.
title: str = lines[0]

if "--p" in sys.argv:
    title += f" Order: {order}."

if "--single" in sys.argv: # Elements only.
    # Plot.
    fig, axes = plt.subplots()
    fig.suptitle(title)

    # Elements.
    if "--p" in sys.argv:
        axes.plot(elements, elements_comparison[0], linewidth=1.0, linestyle=":", color="red", label=f"Order: {order - 1}") # Comparison.
        axes.plot(elements, elements_comparison[1], linewidth=1.0, linestyle="-", color="red", label=f"Order: {order}") # Comparison.
        axes.plot(elements, elements_comparison[2], linewidth=1.0, linestyle="--", color="red", label=f"Order: {order + 1}") # Comparison.
    axes.plot(elements, l2_error, color="black") # Errors.

    # Loglog scale.
    axes.set_xscale("log")
    axes.set_yscale("log")

    # Ticks.
    axes.xaxis.set_minor_formatter(NullFormatter())
    axes.yaxis.set_minor_formatter(NullFormatter())
    axes.set_xticks(elements_ticks, labels=elements_labels)
    axes.set_yticks(errors_ticks, labels=errors_labels)

    # Label.
    axes.set_xlabel("Elements")

    # Legend.
    axes.legend()

    plt.show()
    sys.exit(0)

# Plot.
fig, axes = plt.subplots(1, 2)
fig.suptitle(title)

# Elements.
if "--p" in sys.argv:
    axes[0].plot(elements, elements_comparison[0], linewidth=1.0, linestyle=":", color="red", label=f"Order: {order - 1}") # Comparison.
    axes[0].plot(elements, elements_comparison[1], linewidth=1.0, linestyle="-", color="red", label=f"Order: {order}") # Comparison.
    axes[0].plot(elements, elements_comparison[2], linewidth=1.0, linestyle="--", color="red", label=f"Order: {order + 1}") # Comparison.
axes[0].plot(elements, l2_error, color="black") # Errors.

# Sizes.
if "--p" in sys.argv:
    axes[1].plot(sizes, sizes_comparison[0], linewidth=1.0, linestyle=":", color="red", label=f"Order: {order - 1}") # Comparison.
    axes[1].plot(sizes, sizes_comparison[1], linewidth=1.0, linestyle="-", color="red", label=f"Order: {order}") # Comparison.
    axes[1].plot(sizes, sizes_comparison[2], linewidth=1.0, linestyle="--", color="red", label=f"Order: {order + 1}") # Comparison.
axes[1].plot(sizes, l2_error, color="black") # Errors.

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
axes[0].set_yticks(errors_ticks, labels=errors_labels)
axes[1].set_xticks(sizes_ticks, labels=sizes_labels)
axes[1].set_yticks(errors_ticks, labels=errors_labels)

# Labels.
axes[0].set_xlabel("Elements")
axes[1].set_xlabel("Size")

# Legends.
axes[0].legend()
axes[1].legend()

plt.show()