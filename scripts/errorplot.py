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
elements_second: list[int] = []
sizes: list[float] = []

l2_error: list[float] = []
dg_error: list[float] = []

l2_error_second: list[float] = [] # Second file, if present.
dg_error_second: list[float] = [] # Second file, if present.

second: bool = False
order: int = 0

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

        elif "Size" in line:
            sizes.append(float(data[-1]))

        elif "L2 Error" in line:
            l2_error.append(float(data[-1]))

        elif "DG Error" in line:
            dg_error.append(float(data[-1]))

        elif "Degree" in line:
            order = int(data[-1])

    except ValueError:
        continue

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

elements_comparison_l2[0] = [(value / elements[-1]) ** -((order + 1) / 2 - 0.5) * l2_error[-1] for value in elements]
sizes_comparison_l2[0] = [(value / sizes[0]) ** order * l2_error[0] for value in sizes]

elements_comparison_l2[1] = [(value / elements[-1]) ** -((order + 1) / 2) * l2_error[-1] for value in elements]
sizes_comparison_l2[1] = [(value / sizes[0]) ** (order + 1) * l2_error[0] for value in sizes]

elements_comparison_l2[2] = [(value / elements[-1]) ** -((order + 1) / 2 + 0.5) * l2_error[-1] for value in elements]
sizes_comparison_l2[2] = [(value / sizes[0]) ** (order + 2) * l2_error[0] for value in sizes]

# DG.
elements_comparison_dg: list[list[float]] = [None] * 3
sizes_comparison_dg: list[list[float]] = [None] * 3

elements_comparison_dg[0] = [(value / elements[-1]) ** -(order / 2 - 0.5) * dg_error[-1] for value in elements]
sizes_comparison_dg[0] = [(value / sizes[0]) ** (order - 1) * dg_error[0] for value in sizes]

elements_comparison_dg[1] = [(value / elements[-1]) ** -(order / 2) * dg_error[-1] for value in elements]
sizes_comparison_dg[1] = [(value / sizes[0]) ** order * dg_error[0] for value in sizes]

elements_comparison_dg[2] = [(value / elements[-1]) ** -(order / 2 + 0.5) * dg_error[-1] for value in elements]
sizes_comparison_dg[2] = [(value / sizes[0]) ** (order + 1) * dg_error[0] for value in sizes]

# Ticks.
elements_ticks = [elements[0], elements[-1]]
elements_labels = [str(tick) for tick in elements_ticks]

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

if "--p" in sys.argv:
    title += f" Order: {order}."

if "--elements" in sys.argv: # Elements only.
    # Plot.
    fig, axes = plt.subplots(1, 2)
    fig.suptitle(title)

    # Elements.
    if "--p" in sys.argv:
        axes[0].plot(elements, elements_comparison_l2[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: -{(order + 1) / 2 - 0.5}") # Comparison.
        axes[0].plot(elements, elements_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: -{(order + 1) / 2}") # Comparison.
        axes[0].plot(elements, elements_comparison_l2[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: -{(order + 1) / 2 + 0.5}") # Comparison.
    axes[0].plot(elements, l2_error, color="black", marker="*", linewidth=1.0) # Errors.

    if second:
        axes[0].plot(elements_second, l2_error_second, color="green", marker="*", linewidth=1.0) # Errors.

    if "--p" in sys.argv:
        axes[1].plot(elements, elements_comparison_dg[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: -{order / 2 - 0.5}") # Comparison.
        axes[1].plot(elements, elements_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: -{order / 2}") # Comparison.
        axes[1].plot(elements, elements_comparison_dg[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: -{order / 2 + 0.5}") # Comparison.
    axes[1].plot(elements_second, dg_error, color="black", marker="*", linewidth=1.0) # Errors.

    if second:
        axes[1].plot(elements_second, dg_error_second, color="green", marker="*", linewidth=1.0) # Errors.

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

# Plot.
fig, axes = plt.subplots(2, 2)
fig.suptitle(title)

# Elements.
if "--p" in sys.argv:
    axes[0, 0].plot(elements, elements_comparison_l2[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: -{(order + 1) / 2 -0.5}") # Comparison.
    axes[0, 0].plot(elements, elements_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: -{(order + 1) / 2}") # Comparison.
    axes[0, 0].plot(elements, elements_comparison_l2[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: -{(order + 1) / 2 + 0.5}") # Comparison.
axes[0, 0].plot(elements, l2_error, color="black", marker="*", linewidth=1.0) # Errors.

if "--p" in sys.argv:
    axes[1, 0].plot(elements, elements_comparison_dg[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: -{order / 2 - 0.5}") # Comparison.
    axes[1, 0].plot(elements, elements_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: -{order / 2}") # Comparison.
    axes[1, 0].plot(elements, elements_comparison_dg[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: -{order / 2 + 0.5}") # Comparison.
axes[1, 0].plot(elements, dg_error, color="black", marker="*", linewidth=1.0) # Errors.

# Sizes.
if "--p" in sys.argv:
    axes[0, 1].plot(sizes, sizes_comparison_l2[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: {order}") # Comparison.
    axes[0, 1].plot(sizes, sizes_comparison_l2[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: {order + 1}") # Comparison.
    axes[0, 1].plot(sizes, sizes_comparison_l2[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: {order + 2}") # Comparison.
axes[0, 1].plot(sizes, l2_error, color="black", marker="*", linewidth=1.0) # Errors.

if "--p" in sys.argv:
    axes[1, 1].plot(sizes, sizes_comparison_dg[0], linewidth=1.0, alpha=0.5, linestyle=":", color="red", label=f"Order: {order - 1}") # Comparison.
    axes[1, 1].plot(sizes, sizes_comparison_dg[1], linewidth=1.0, alpha=1, linestyle="-", color="red", label=f"Order: {order}") # Comparison.
    axes[1, 1].plot(sizes, sizes_comparison_dg[2], linewidth=1.0, alpha=0.5, linestyle="--", color="red", label=f"Order: {order + 1}") # Comparison.
axes[1, 1].plot(sizes, dg_error, color="black", marker="*", linewidth=1.0) # Errors.

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