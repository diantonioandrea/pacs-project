#!/usr/bin/env python3

"""
@file evs_tikz.py
@brief errorvsize TikZ Python wrapper.
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
    file = open("templates/sizes_tikz.tex", "r+")
    template: str = file.read()
    file.close()

except FileNotFoundError:
    print("Templates not found.")
    sys.argv(-1)

# Dofs.
sizes: list[float] = []

# Errors.
l2_errors: list[float] = []
dg_errors: list[float] = []

# Degree.
degree: int = -1

# Plots.
l2_plots: list[str] = []
dg_plots: list[str] = []

try:
    if sys.argv[1].split(".")[-1] != "error":
        raise

    file = open(sys.argv[1], "r+")
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

# Comparison.
l2_comparison: list[float] = []
dg_comparison: list[float] = []

for size in sizes:
    l2_comparison.append((size / sizes[-1]) ** (degree + 1) * l2_errors[-1])
    dg_comparison.append((size / sizes[-1]) ** degree * dg_errors[-1])

# Plots.

# Errors.
l2_errors_string: str = "\n\\addplot[black, mark=*] coordinates "
l2_errors_string += "{" + " ".join((f"({sizes[j]},{l2_errors[j]})" for j in range(len(sizes)))) + "};"
l2_errors_string += "\n"
l2_errors_string += "\\addlegendentry{$\\LT$ Error}\n"

dg_errors_string: str = "\n\\addplot[black, mark=*] coordinates "
dg_errors_string += "{" + " ".join((f"({sizes[j]},{dg_errors[j]})" for j in range(len(sizes)))) + "};"
dg_errors_string += "\n"
dg_errors_string += "\\addlegendentry{$DG$ Error}\n"

# Comparison.
l2_comparison_string: str = f"\n\\addplot[black, dashed] coordinates "
l2_comparison_string += "{"+ f"({sizes[0]},{l2_comparison[0]}) ({sizes[-1]},{l2_comparison[-1]})" + "};\n"
l2_comparison_string += "\\addlegendentry{$\\mathcal{O}(h^{" + str(degree + 1) + "})$}\n"

dg_comparison_string: str = f"\n\\addplot[black, dashed] coordinates "
dg_comparison_string += "{"+ f"({sizes[0]},{dg_comparison[0]}) ({sizes[-1]},{dg_comparison[-1]})" + "};\n"
dg_comparison_string += "\\addlegendentry{$\\mathcal{O}(h^{" + str(degree) + "})$}\n"

# Strings.
l2_plots.append(l2_errors_string)
l2_plots.append(l2_comparison_string)

dg_plots.append(dg_errors_string)
dg_plots.append(dg_comparison_string)

# Text.
template = template.replace("% PLOTS_L2", "".join(l2_plots))
template = template.replace("% PLOTS_DG", "".join(dg_plots))

# Output.
file = open(sys.argv[1].replace(".error", "_evs.tex"), "w+")
file.write(template)
file.close()