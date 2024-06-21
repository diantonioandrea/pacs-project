#!/usr/bin/env python3

"""
@file errorvdofs.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-06-21

@copyright Copyright (c) 2024
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy
import sys

# Font.
matplotlib.rcParams.update({'font.size': 18})

# Colors.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]

# Dofs.
dofs: list[float] = []

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

        if "Dofs" in line:
            dofs.append(int(data[-1]))

        elif "L2 Error" in line:
            l2_errors.append(float(data[-1]))

        elif "DG Error" in line:
            dg_errors.append(float(data[-1]))

    except ValueError:
        continue

# Dofs.
dofs = [dof ** 0.5 for dof in dofs]

# Comparison.
l2_comparison: list[float] = []
dg_comparison: list[float] = []

l2_interp = numpy.polyfit(numpy.log(dofs), numpy.log(l2_errors), 1)
dg_interp = numpy.polyfit(numpy.log(dofs), numpy.log(dg_errors), 1)

for dof in dofs:
    l2_comparison.append(numpy.exp(l2_interp[1]) * dof ** l2_interp[0])
    dg_comparison.append(numpy.exp(dg_interp[1]) * dof ** dg_interp[0])

# Ticks.
dofs_ticks = [dofs[0], dofs[-1]]

l2_ticks = [l2_errors[0], l2_errors[-1]]
dg_ticks = [dg_errors[0], dg_errors[-1]]

# Labels.
dofs_labels = [f"{tick:.3f}" for tick in dofs_ticks]

l2_labels = [f"{tick:.1e}" for tick in l2_ticks]
dg_labels = [f"{tick:.1e}" for tick in dg_ticks]

# Plot.
fig, axes = plt.subplots(1, 2)
fig.suptitle("L2 and DG errors vs. DOFs")

# L2.
axes[0].plot(dofs, l2_errors, color=black, marker="*", linewidth=3, label="L2 error.") # Error.
axes[0].plot(dofs, l2_comparison, color=black, linestyle="--", linewidth=1.5, alpha=0.5, label=f"Interpolant, degree: {l2_interp[0]:.1f}") # Comparison.

# DG.
axes[1].plot(dofs, dg_errors, color=black, marker="*", linewidth=3, label="DG error.") # Error.
axes[1].plot(dofs, dg_comparison, color=black, linestyle="--", linewidth=1.5, alpha=0.5, label=f"Interpolant, degree: {dg_interp[0]:.1f}") # Comparison.

# Comparison.
if len(sys.argv) == 3:

    # File.
    try:
        if sys.argv[2].split(".")[-1] != "error":
            raise

        file = open(sys.argv[2], "r+")
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

    # Reset.

    # Dofs.
    dofs: list[float] = []

    # Errors.
    l2_errors: list[float] = []
    dg_errors: list[float] = []

    for line in lines:
        try:
            data: list[str] = line.split(" ")

            if "Dofs" in line:
                dofs.append(int(data[-1]))

            elif "L2 Error" in line:
                l2_errors.append(float(data[-1]))

            elif "DG Error" in line:
                dg_errors.append(float(data[-1]))

        except ValueError:
            continue

    # Dofs.
    dofs = [dof ** 0.5 for dof in dofs]

    # Comparison.
    l2_comparison: list[float] = []
    dg_comparison: list[float] = []

    l2_interp = numpy.polyfit(numpy.log(dofs), numpy.log(l2_errors), 1)
    dg_interp = numpy.polyfit(numpy.log(dofs), numpy.log(dg_errors), 1)

    for dof in dofs:
        l2_comparison.append(numpy.exp(l2_interp[1]) * dof ** l2_interp[0])
        dg_comparison.append(numpy.exp(dg_interp[1]) * dof ** dg_interp[0])

    # Plot.

    # L2.
    axes[0].plot(dofs, l2_errors, color=red, marker="*", linewidth=3, label="L2 error (C).") # Error.
    axes[0].plot(dofs, l2_comparison, color=red, linestyle="--", linewidth=1.5, alpha=0.5, label=f"Interpolant, degree: {l2_interp[0]:.1f} (C)") # Comparison.

    # DG.
    axes[1].plot(dofs, dg_errors, color=red, marker="*", linewidth=3, label="DG error (C).") # Error.
    axes[1].plot(dofs, dg_comparison, color=red, linestyle="--", linewidth=1.5, alpha=0.5, label=f"Interpolant, degree: {dg_interp[0]:.1f} (C)") # Comparison.

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

# Title.
axes[0].set_title("L2 error")
axes[1].set_title("DG error")

# Labels.
axes[0].set_xticks(dofs_ticks, labels=dofs_labels)
axes[1].set_xticks(dofs_ticks, labels=dofs_labels)
axes[0].set_yticks(l2_ticks, labels=l2_labels)
axes[1].set_yticks(dg_ticks, labels=dg_labels)

axes[1].yaxis.tick_right()
axes[1].yaxis.set_label_position("right")

# Output.
plt.show()