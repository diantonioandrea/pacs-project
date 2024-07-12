#!/usr/bin/env python3

"""
@file evd_tikz.py
@brief hpvdofs.py TikZ Python wrapper.
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-07-12

@copyright Copyright (c) 2024
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys

# Font.
matplotlib.rcParams.update({'font.size': 18})

# Colours.
colours: list[str] = ["solarized-base02", "\\accentcolor"]

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

# Template.
try:
    file = open("templates/hp_tikz.tex", "r+")
    template: str = file.read()
    file.close()

except FileNotFoundError:
    print("Templates not found.")
    sys.argv(-1)

# Plots.
plots: list[str] = []

# Multiple data.
for index in range(1, len(sys.argv)):

    # Dofs.
    dofs: list[float] = []

    # Errors.
    errors: list[float] = []

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

        except ValueError:
            continue

    # Scaling.
    dofs = [dof ** (1/3) for dof in dofs]

    # Plots.

    # Errors.
    errors_string: str = f"\n\\addplot[{colours[index - 1]}, mark=*] coordinates "
    errors_string += "{" + " ".join((f"({dofs[j]},{errors[j]})" for j in range(len(dofs)))) + "};"
    errors_string += "\n"
    errors_string += "\\addlegendentry{$DG$ Error}\n"

    # Strings.
    plots.append(errors_string)

# Text.
template = template.replace("% PLOTS", "".join(plots))

# Output.
file = open(sys.argv[1].replace(".error", "_hpvd.tex"), "w+")
file.write(template)
file.close()