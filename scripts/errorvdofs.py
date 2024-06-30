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

# Colors.
black: list[float] = [7 / 255, 54 / 255, 66 / 255]
red: list[float] = [220 / 255, 50 / 255, 47 / 255]

# Dofs.
dofs: list[float] = []

# Errors.
errors: list[float] = []

# Degree.
degree: int = -1

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

        elif "DG Error" in line:
            errors.append(float(data[-1]))

        elif "Degree" in line:
            degree = int(data[-1])

    except ValueError:
        continue

# Exponent.
exponent: float = -degree / 2

# Plot.
fig, axes = plt.subplots()
fig.suptitle("Errors vs. DOFs")

# DG.
axes.plot(dofs, errors, color=black, marker="*", linewidth=1, label="DG") # Error.

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
    errors: list[float] = []

    for line in lines:
        try:
            data: list[str] = line.split(" ")

            if "Dofs" in line:
                dofs.append(int(data[-1]))

            elif "DG Error" in line:
                errors.append(float(data[-1]))

        except ValueError:
            continue

    # DG.
    axes.plot(dofs, errors, color=red, marker="*", linewidth=1, label="DG (C)") # Error.

# Comparison.
dofs_comparison = [errors[0] * (dof / dofs[0]) ** exponent for dof in dofs]
axes.plot(dofs, dofs_comparison, color=black, linestyle="--", linewidth=0.5, alpha=0.5, label="$DOFs^{" + str(exponent) + "}$")

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