#!/usr/bin/env python3

"""
@file evd_tikz.py
@brief errorvdofs TikZ Python wrapper.
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-07-03

@copyright Copyright (c) 2024
"""

import sys

# File.
if len(sys.argv) <= 1:
    print(f"Usage: {sys.argv[0]} /path/to/file.")
    sys.exit(0)

# Template.
try:
    file = open("templates/dofs_tikz.tex", "r+")
    template: str = file.read()
    file.close()

except FileNotFoundError:
    print("Templates not found.")
    sys.argv(-1)

# Plots.
plots: list[str] = []

# Colours.
colours: list[str] = ["solarized-base02", "\\accentcolor"]

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
    exponent: float = -degree/2

    # Plots.

    # Errors.
    errors_string: str = f"\n\\addplot[{colours[index - 1]}, mark=*] coordinates "
    errors_string += "{" + " ".join((f"({dofs[j]},{errors[j]})" for j in range(len(dofs)))) + "};"
    errors_string += "\n"
    errors_string += "\\addlegendentry{$DG$ Error}\n"

    # Comparison.
    dofs_comparison = [errors[-1] * (dof / dofs[-1]) ** exponent for dof in dofs]

    comparison_string: str = f"\n\\addplot[{colours[index - 1]}, dashed] coordinates "
    comparison_string += "{"+ f"({dofs[0]},{dofs_comparison[0]}) ({dofs[-1]},{dofs_comparison[-1]})" + "};\n"
    comparison_string += "\\addlegendentry{$\\mathcal{O}(DOFs^{" + str(exponent) + "})$}\n"

    # Strings.
    plots.append(errors_string)
    plots.append(comparison_string)

# Text.
template = template.replace("% PLOTS", "".join(plots))

# Output.
file = open(sys.argv[1].replace(".error", "_evd.tex"), "w+")
file.write(template)
file.close()