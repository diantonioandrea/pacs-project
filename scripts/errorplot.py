#!/usr/bin/env python3

"""
@file errorplot.py
@author Andrea Di Antonio (github.com/diantonioandrea)
@date 2024-05-19

@copyright Copyright (c) 2024
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy
import sys

elements: list[int] = []
elements_second: list[int] = []
sizes: list[float] = []

dofs: list[float] = []
dofs_second: list[float] = []

l2_error: list[float] = []
dg_error: list[float] = []

l2_error_second: list[float] = [] # Second file, if present.
dg_error_second: list[float] = [] # Second file, if present.

second: bool = False

degrees: list[int] = []
degree: int = 0

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

        if "Elements" in line:
            elements.append(int(data[-1]))
        
        elif "Dofs" in line:
            dofs.append(int(data[-1]))

        elif "Size" in line:
            sizes.append(float(data[-1]))

        elif "L2 Error" in line:
            l2_error.append(float(data[-1]))

        elif "DG Error" in line:
            dg_error.append(float(data[-1]))

        elif "Degree" in line:
            degrees.append(int(data[-1]))

    except ValueError:
        continue

degree = max(degrees)

for index in range(len(sys.argv)):
    if sys.argv[index] == "-c": # Compare.
        try:
            if sys.argv[index + 1].split(".")[-1] != "error":
                raise

            file = open(sys.argv[index + 1], "r+")
            lines: list[str] = file.read().split("\n")
            file.close()

        except FileNotFoundError:
            print("File not found.")
            sys.exit(-1)

        except:
            print("Load a .error file.")
            sys.exit(-1)

        for line in lines:
            try:
                data: list[str] = line.split(" ")
                
                if "Elements" in line:
                    elements_second.append(int(data[-1]))

                if "Dofs" in line:
                    dofs_second.append(int(data[-1]))

                elif "L2 Error" in line:
                    l2_error_second.append(float(data[-1]))

                elif "DG Error" in line:
                    dg_error_second.append(float(data[-1]))

            except ValueError:
                continue
        
        second = True
        break

# Comparisons.

# L2.
elements_comparison_l2: list[list[float]] = [None] * 3
sizes_comparison_l2: list[list[float]] = [None] * 3

elements_comparison_l2[0] = [(value / elements[-1]) ** -((degree + 1) / 2 - 0.5) * l2_error[-1] for value in elements]
sizes_comparison_l2[0] = [(value / sizes[0]) ** degree * l2_error[0] for value in sizes]

elements_comparison_l2[1] = [(value / elements[-1]) ** -((degree + 1) / 2) * l2_error[-1] for value in elements]
sizes_comparison_l2[1] = [(value / sizes[0]) ** (degree + 1) * l2_error[0] for value in sizes]

elements_comparison_l2[2] = [(value / elements[-1]) ** -((degree + 1) / 2 + 0.5) * l2_error[-1] for value in elements]
sizes_comparison_l2[2] = [(value / sizes[0]) ** (degree + 2) * l2_error[0] for value in sizes]

# DG.
elements_comparison_dg: list[list[float]] = [None] * 3
sizes_comparison_dg: list[list[float]] = [None] * 3

elements_comparison_dg[0] = [(value / elements[-1]) ** -(degree / 2 - 0.5) * dg_error[-1] for value in elements]
sizes_comparison_dg[0] = [(value / sizes[0]) ** (degree - 1) * dg_error[0] for value in sizes]

elements_comparison_dg[1] = [(value / elements[-1]) ** -(degree / 2) * dg_error[-1] for value in elements]
sizes_comparison_dg[1] = [(value / sizes[0]) ** degree * dg_error[0] for value in sizes]

elements_comparison_dg[2] = [(value / elements[-1]) ** -(degree / 2 + 0.5) * dg_error[-1] for value in elements]
sizes_comparison_dg[2] = [(value / sizes[0]) ** (degree + 1) * dg_error[0] for value in sizes]

# Dofs.
dofs = [dof ** 0.5 for dof in dofs]

if second:
    dofs_second = [dof ** 0.5 for dof in dofs_second]

# Ticks.
elements_ticks = [elements[0], elements[-1]]
elements_labels = [str(tick) for tick in elements_ticks]

dofs_ticks = [dofs[0], dofs[-1]]
dofs_labels = [f"{tick:.2f}" for tick in dofs_ticks]

sizes_ticks = [sizes[0], sizes[-1]]
sizes_labels = [f"{tick:1.2f}" for tick in sizes_ticks]

errors_l2_ticks = [l2_error[0], l2_error[-1]]
errors_dg_ticks = [dg_error[0], dg_error[-1]]

if second:
    errors_l2_ticks.append(l2_error_second[-1])
    errors_dg_ticks.append(dg_error_second[-1])

    errors_l2_ticks.sort()
    errors_dg_ticks.sort()

errors_l2_labels = [f"{tick:1.2e}" for tick in errors_l2_ticks]
errors_dg_labels = [f"{tick:1.2e}" for tick in errors_dg_ticks]

if "--degrees" in sys.argv: # Degrees only.

    assert min(elements) == max(elements)
    title += f" {max(elements)} elements."

    fig, axes = plt.subplots(1, 2)
    fig.suptitle(title)

    axes[0].plot(degrees, l2_error, color="black", marker="*", linewidth=1.5) # Errors.
    axes[1].plot(degrees, dg_error, color="black", marker="*", linewidth=1.5) # Errors.

    if "--interp" in sys.argv:
        l2_interp = numpy.polyfit(degrees, numpy.log(l2_error), 1)
        dg_interp = numpy.polyfit(degrees, numpy.log(dg_error), 1)

        l2_interp_tail = numpy.polyfit(degrees[-3:], numpy.log(l2_error[-3:]), 1)
        dg_interp_tail = numpy.polyfit(degrees[-3:], numpy.log(dg_error[-3:]), 1)

        l2_comparison: list[float] = [numpy.exp(deg * l2_interp[0] + l2_interp[1]) for deg in degrees]
        dg_comparison: list[float] = [numpy.exp(deg * dg_interp[0] + dg_interp[1]) for deg in degrees]

        l2_comparison_tail: list[float] = [numpy.exp(deg * l2_interp_tail[0] + l2_interp_tail[1]) for deg in degrees]
        dg_comparison_tail: list[float] = [numpy.exp(deg * dg_interp_tail[0] + dg_interp_tail[1]) for deg in degrees]

        axes[0].plot(degrees, l2_comparison, linewidth=0.75, alpha=0.5, linestyle="-.", color="blue", label=f"Interpolant: {l2_interp[0]:.2f}")
        axes[1].plot(degrees, dg_comparison, linewidth=0.75, alpha=0.5, linestyle="-.", color="blue", label=f"Interpolant: {dg_interp[0]:.2f}")

        axes[0].plot(degrees, l2_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="-.", color="purple", label=f"Tail interpolant: {l2_interp_tail[0]:.2f}")
        axes[1].plot(degrees, dg_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="-.", color="purple", label=f"Tail interpolant: {dg_interp_tail[0]:.2f}")

    # Parameters.
    for j in range(2):

        # Loglog scale.
        axes[j].set_yscale("log")

        # Ticks.
        axes[j].xaxis.set_minor_formatter(NullFormatter())
        axes[j].yaxis.set_minor_formatter(NullFormatter())

        # Labels.
        axes[j].set_xlabel("Degree")

        # Legend.
        if "--interp" in sys.argv:
            axes[j].legend()

    # Titles.
    axes[0].set_title("L2 Error")
    axes[1].set_title("DG Error")

    # Ticks.
    axes[0].set_yticks(errors_l2_ticks, labels=errors_l2_labels)
    axes[1].set_yticks(errors_dg_ticks, labels=errors_dg_labels)

    plt.show()
    sys.exit(0)

if "--p" in sys.argv:
    title += f" Degree: {degree}."

if "--elements" in sys.argv: # Elements only.
    # Plot.
    fig, axes = plt.subplots(1, 2)
    fig.suptitle(title)

    # Elements.
    if "--p" in sys.argv:
        axes[0].plot(elements, elements_comparison_l2[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: -{(degree + 1) / 2 - 0.5}") # Comparison.
        axes[0].plot(elements, elements_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: -{(degree + 1) / 2}") # Comparison.
        axes[0].plot(elements, elements_comparison_l2[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: -{(degree + 1) / 2 + 0.5}") # Comparison.

    axes[0].plot(elements, l2_error, color="black", marker="*", linewidth=1.5) # Errors.

    if "--interp" in sys.argv:
        l2_interp = numpy.polyfit(numpy.log(elements), numpy.log(l2_error), 1)
        l2_comparison: list[float] = [numpy.exp(l2_interp[1]) * element ** l2_interp[0] for element in elements]

        l2_interp_tail = numpy.polyfit(numpy.log(elements[-3:]), numpy.log(l2_error[-3:]), 1)
        l2_comparison_tail: list[float] = [numpy.exp(l2_interp_tail[1]) * element ** l2_interp_tail[0] for element in elements]

        axes[0].plot(elements, l2_comparison, linewidth=0.75, alpha=0.5, linestyle="-.", color="blue", label=f"Interpolant: {l2_interp[0]:.2f}")
        axes[0].plot(elements, l2_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="-.", color="purple", label=f"Tail interpolant: {l2_interp_tail[0]:.2f}")

    if second:
        axes[0].plot(elements_second, l2_error_second, color="green", marker="*", linewidth=1.5) # Errors.

        if "--interp" in sys.argv:
            l2_interp = numpy.polyfit(numpy.log(elements_second), numpy.log(l2_error_second), 1)
            l2_comparison: list[float] = [numpy.exp(l2_interp[1]) * element ** l2_interp[0] for element in elements_second]

            l2_interp_tail = numpy.polyfit(numpy.log(elements_second[-3:]), numpy.log(l2_error_second[-3:]), 1)
            l2_comparison_tail: list[float] = [numpy.exp(l2_interp_tail[1]) * element ** l2_interp_tail[0] for element in elements_second]

            axes[0].plot(elements_second, l2_comparison, linewidth=0.75, alpha=0.5, linestyle="--", color="blue", label=f"Interpolant: {l2_interp[0]:.2f}")
            axes[0].plot(elements_second, l2_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="--", color="purple", label=f"Tail interpolant: {l2_interp_tail[0]:.2f}")

    if "--p" in sys.argv:
        axes[1].plot(elements, elements_comparison_dg[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: -{degree / 2 - 0.5}") # Comparison.
        axes[1].plot(elements, elements_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: -{degree / 2}") # Comparison.
        axes[1].plot(elements, elements_comparison_dg[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: -{degree / 2 + 0.5}") # Comparison.

    axes[1].plot(elements, dg_error, color="black", marker="*", linewidth=1.5) # Errors.

    if "--interp" in sys.argv:
        dg_interp = numpy.polyfit(numpy.log(elements), numpy.log(dg_error), 1)
        dg_comparison: list[float] = [numpy.exp(dg_interp[1]) * element ** dg_interp[0] for element in elements]

        dg_interp_tail = numpy.polyfit(numpy.log(elements[-3:]), numpy.log(dg_error[-3:]), 1)
        dg_comparison_tail: list[float] = [numpy.exp(dg_interp_tail[1]) * element ** dg_interp_tail[0] for element in elements]

        axes[1].plot(elements, dg_comparison, linewidth=0.75, alpha=0.5, linestyle="-.", color="blue", label=f"Interpolant: {dg_interp[0]:.2f}")
        axes[1].plot(elements, dg_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="-.", color="purple", label=f"Tail interpolant: {dg_interp_tail[0]:.2f}")

    if second:
        axes[1].plot(elements_second, dg_error_second, color="green", marker="*", linewidth=1.5) # Errors.

        if "--interp" in sys.argv:
            dg_interp = numpy.polyfit(numpy.log(elements_second), numpy.log(dg_error_second), 1)
            dg_comparison: list[float] = [numpy.exp(dg_interp[1]) * element ** dg_interp[0] for element in elements_second]

            dg_interp_tail = numpy.polyfit(numpy.log(elements_second[-3:]), numpy.log(dg_error_second[-3:]), 1)
            dg_comparison_tail: list[float] = [numpy.exp(dg_interp_tail[1]) * element ** dg_interp_tail[0] for element in elements_second]

            axes[1].plot(elements_second, dg_comparison, linewidth=0.75, alpha=0.5, linestyle="--", color="blue", label=f"Interpolant: {dg_interp[0]:.2f}")
            axes[1].plot(elements_second, dg_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="--", color="purple", label=f"Tail interpolant: {dg_interp_tail[0]:.2f}")

    # Parameters.
    for j in range(2):

        # Loglog scale.
        axes[j].set_xscale("log")
        axes[j].set_yscale("log")

        # Ticks.
        axes[j].xaxis.set_minor_formatter(NullFormatter())
        axes[j].yaxis.set_minor_formatter(NullFormatter())

        # Labels.
        axes[j].set_xlabel("Elements")

        # Legend.
        if "--p" in sys.argv or "--interp" in sys.argv:
            axes[j].legend()
    
    # Ticks.
    axes[0].set_xticks(elements_ticks, labels=elements_labels)
    axes[0].set_yticks(errors_l2_ticks, labels=errors_l2_labels)
    axes[1].set_xticks(elements_ticks, labels=elements_labels)
    axes[1].set_yticks(errors_dg_ticks, labels=errors_dg_labels)

    # Titles.
    axes[0].set_title("L2 Error")
    axes[1].set_title("DG Error")

    plt.show()
    sys.exit(0)

if "--dofs" in sys.argv: # Dofs comparison.
    # Plot.
    fig, axes = plt.subplots(1, 2)
    fig.suptitle(title)

    axes[0].plot(dofs, l2_error, color="black", marker="*", linewidth=1.5) # Errors.
    axes[1].plot(dofs, dg_error, color="black", marker="*", linewidth=1.5) # Errors.

    if "--interp" in sys.argv:
        l2_interp = numpy.polyfit(numpy.log(dofs), numpy.log(l2_error), 1)
        l2_comparison: list[float] = [numpy.exp(l2_interp[1]) * dof ** l2_interp[0] for dof in dofs]

        l2_interp_tail = numpy.polyfit(numpy.log(dofs[-3:]), numpy.log(l2_error[-3:]), 1)
        l2_comparison_tail: list[float] = [numpy.exp(l2_interp_tail[1]) * dof ** l2_interp_tail[0] for dof in dofs]

        axes[0].plot(dofs, l2_comparison, linewidth=0.75, alpha=0.5, linestyle="-.", color="blue", label=f"Interpolant: {l2_interp[0]:.2f}")
        axes[0].plot(dofs, l2_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="-.", color="purple", label=f"Tail interpolant: {l2_interp_tail[0]:.2f}")

    if second:
        axes[0].plot(dofs_second, l2_error_second, color="green", marker="*", linewidth=1.5) # Errors.
        axes[1].plot(dofs_second, dg_error_second, color="green", marker="*", linewidth=1.5) # Errors.

        if "--interp" in sys.argv:
            l2_interp = numpy.polyfit(numpy.log(dofs_second), numpy.log(l2_error_second), 1)
            l2_comparison: list[float] = [numpy.exp(l2_interp[1]) * dof ** l2_interp[0] for dof in dofs_second]

            l2_interp_tail = numpy.polyfit(numpy.log(dofs_second[-3:]), numpy.log(l2_error_second[-3:]), 1)
            l2_comparison_tail: list[float] = [numpy.exp(l2_interp_tail[1]) * dof ** l2_interp_tail[0] for dof in dofs_second]

            axes[0].plot(dofs_second, l2_comparison, linewidth=0.75, alpha=0.5, linestyle="--", color="blue", label=f"Interpolant: {l2_interp[0]:.2f}")
            axes[0].plot(dofs_second, l2_comparison_tail, linewidth=0.75, alpha=0.5, linestyle="--", color="purple", label=f"Tail interpolant: {l2_interp_tail[0]:.2f}")

    # Parameters.
    for j in range(2):

        # Loglog scale.
        axes[j].set_xscale("log")
        axes[j].set_yscale("log")

        # Ticks.
        axes[j].xaxis.set_minor_formatter(NullFormatter())
        axes[j].yaxis.set_minor_formatter(NullFormatter())

        if "--interp" in sys.argv:
            axes[j].legend()

        # Labels.
        axes[j].set_xlabel("DOFs ^ 1/2")
    
    # Ticks.
    axes[0].set_xticks(dofs_ticks, labels=dofs_labels)
    axes[0].set_yticks(errors_l2_ticks, labels=errors_l2_labels)
    axes[1].set_xticks(dofs_ticks, labels=dofs_labels)
    axes[1].set_yticks(errors_dg_ticks, labels=errors_dg_labels)

    # Titles.
    axes[0].set_title("L2 Error")
    axes[1].set_title("DG Error")

    plt.show()
    sys.exit(0)

# Plot.
fig, axes = plt.subplots(2, 2)
fig.suptitle(title)

# Elements.
if "--p" in sys.argv:
    axes[0, 0].plot(elements, elements_comparison_l2[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: -{(degree + 1) / 2 -0.5}") # Comparison.
    axes[0, 0].plot(elements, elements_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: -{(degree + 1) / 2}") # Comparison.
    axes[0, 0].plot(elements, elements_comparison_l2[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: -{(degree + 1) / 2 + 0.5}") # Comparison.
axes[0, 0].plot(elements, l2_error, color="black", marker="*", linewidth=1.5) # Errors.

if "--p" in sys.argv:
    axes[1, 0].plot(elements, elements_comparison_dg[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: -{degree / 2 - 0.5}") # Comparison.
    axes[1, 0].plot(elements, elements_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: -{degree / 2}") # Comparison.
    axes[1, 0].plot(elements, elements_comparison_dg[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: -{degree / 2 + 0.5}") # Comparison.
axes[1, 0].plot(elements, dg_error, color="black", marker="*", linewidth=1.5) # Errors.

# Sizes.
if "--p" in sys.argv:
    axes[0, 1].plot(sizes, sizes_comparison_l2[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: {degree}") # Comparison.
    axes[0, 1].plot(sizes, sizes_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: {degree + 1}") # Comparison.
    axes[0, 1].plot(sizes, sizes_comparison_l2[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: {degree + 2}") # Comparison.
axes[0, 1].plot(sizes, l2_error, color="black", marker="*", linewidth=1.5) # Errors.

if "--p" in sys.argv:
    axes[1, 1].plot(sizes, sizes_comparison_dg[0], linewidth=0.75, alpha=0.5, linestyle=":", color="red", label=f"Degree: {degree - 1}") # Comparison.
    axes[1, 1].plot(sizes, sizes_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Degree: {degree}") # Comparison.
    axes[1, 1].plot(sizes, sizes_comparison_dg[2], linewidth=0.75, alpha=0.5, linestyle="--", color="red", label=f"Degree: {degree + 1}") # Comparison.
axes[1, 1].plot(sizes, dg_error, color="black", marker="*", linewidth=1.5) # Errors.

# Parameters.
for j in range(2):
    for k in range(2):

        # Loglog scale.
        axes[j, k].set_xscale("log")
        axes[j, k].set_yscale("log")

        # Formatter.
        axes[j, k].xaxis.set_minor_formatter(NullFormatter())
        axes[j, k].yaxis.set_minor_formatter(NullFormatter())

        # Legend.
        if "--p" in sys.argv:
            axes[j, k].legend()

    # Ticks.
    axes[j, 0].set_xticks(elements_ticks, labels=elements_labels)
    axes[j, 1].set_xticks(sizes_ticks, labels=sizes_labels)
    axes[0, j].set_yticks(errors_l2_ticks, labels=errors_l2_labels)
    axes[1, j].set_yticks(errors_dg_ticks, labels=errors_dg_labels)
    
    # Labels.
    axes[j, 0].set_xlabel("Elements")
    axes[j, 1].set_xlabel("Size")

    # Titles.
    axes[0, j].set_title("L2 Error")
    axes[1, j].set_title("DG Error")

plt.show()