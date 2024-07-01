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
import sys

# Font.
matplotlib.rcParams.update({'font.size': 18})

# Colours.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]

colours: list[list[float]] = [black, red]

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

# Plot.
fig, axes = plt.subplots()
fig.suptitle("Errors vs. DOFs")

# Multiple data.
for index in range(1, len(sys.argv)):

    # Dofs.
    dofs: list[float] = []

    # Errors.
    errors: list[float] = []

    # Degree.
    degree: int = -1

    try:
        if sys.argv[index].split(".")[-1] != "error":
            raise

        file = open(sys.argv[index], "r+")
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

            elif "DG Error" in line:
                errors.append(float(data[-1]))

            elif "Degree" in line:
                degree = int(data[-1])

        except ValueError:
            continue

    # Exponent.
    exponent: float = -degree / 2

    # DG.
    axes.plot(dofs, errors, color=colours[index - 1], marker="*", linewidth=1, label="DG") # Error.

    # Comparison.
    dofs_comparison = [errors[-1] * (dof / dofs[-1]) ** exponent for dof in dofs]
    axes.plot(dofs, dofs_comparison, color=colours[index - 1], linestyle="--", linewidth=0.5, alpha=0.5, label="$DOFs^{" + str(exponent) + "}$")

# Parameters.

# Loglog scale.
axes.set_xscale("log")
axes.set_yscale("log")

# Ticks.
axes.xaxis.set_minor_formatter(NullFormatter())
axes.yaxis.set_minor_formatter(NullFormatter())

# Legend.
axes.legend(loc="best")

# Output.
plt.show()