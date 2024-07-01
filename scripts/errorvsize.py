#!/usr/bin/env python3

"""
@file errorvsize.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-06-19

@copyright Copyright (c) 2024
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy
import sys

# Font.
matplotlib.rcParams.update({'font.size': 18})

# Sizes.
sizes: list[float] = []

# Degree (comparison).
degree: int = -1

# Errors.
l2_errors: list[float] = []
dg_errors: list[float] = []

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

try:
    if sys.argv[1].split(".")[-1] != "error":
        raise

    file = open(sys.argv[1], "r+")
    lines: list[str] = file.read().split("\n")
    file.close()

    # Title.
    title: str = lines[0]

except FileNotFoundError:
    print("File not found.")
    sys.exit(-1)

except:
    print("Load a .error file.")
    sys.exit(-1)

# Simple scraping.
for line in lines:
    try:
        data: list[str] = line.split(" ")

        if "Size" in line:
            sizes.append(float(data[-1]))

        elif "L2 Error" in line:
            l2_errors.append(float(data[-1]))

        elif "DG Error" in line:
            dg_errors.append(float(data[-1]))

        elif "Degree" in line:
            degree = int(data[-1])

    except ValueError:
        continue

# Colors.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]

# Comparison.
l2_comparison: list[float] = []
dg_comparison: list[float] = []

l2_interp_comparison: list[float] = []
dg_interp_comparison: list[float] = []

for size in sizes:
    l2_comparison.append((size / sizes[-1]) ** (degree + 1) * l2_errors[-1])
    dg_comparison.append((size / sizes[-1]) ** degree * dg_errors[-1])

# Ticks.
sizes_ticks = [sizes[0], sizes[-1]]

l2_ticks = [l2_errors[0], l2_errors[-1]]
dg_ticks = [dg_errors[0], dg_errors[-1]]

# Labels.
sizes_labels = [f"{tick:.3f}" for tick in sizes_ticks]

l2_labels = [f"{tick:.1e}" for tick in l2_ticks]
dg_labels = [f"{tick:.1e}" for tick in dg_ticks]

# Plot.
fig, axes = plt.subplots(1, 2)
fig.suptitle("L2 and DG errors vs. size")

# L2.
axes[0].plot(sizes, l2_errors, color=black, marker="*", linewidth=1, label="$L^2$ error") # Error.
axes[0].plot(sizes, l2_comparison, color=red, linestyle="--", linewidth=0.5, label=f"$h^{degree + 1}$") # Comparison.

# DG.
axes[1].plot(sizes, dg_errors, color=black, marker="*", linewidth=1, label="$DG$ error") # Error.
axes[1].plot(sizes, dg_comparison, color=red, linestyle="--", linewidth=0.5, label=f"$h^{degree}$") # Comparison.

# Parameters.
for j in range(2):

    # Loglog scale.
    axes[j].set_xscale("log")
    axes[j].set_yscale("log")

    # Ticks.
    axes[j].xaxis.set_minor_formatter(NullFormatter())
    axes[j].yaxis.set_minor_formatter(NullFormatter())

    # Legend.
    axes[j].legend(loc="best")

# # Title.
# axes[0].set_title("L2 error")
# axes[1].set_title("DG error")

# # Labels.
# axes[0].set_xticks(sizes_ticks, labels=sizes_labels)
# axes[1].set_xticks(sizes_ticks, labels=sizes_labels)
# axes[0].set_yticks(l2_ticks, labels=l2_labels)
# axes[1].set_yticks(dg_ticks, labels=dg_labels)

axes[1].yaxis.tick_right()
axes[1].yaxis.set_label_position("right")

# Output.
plt.show()