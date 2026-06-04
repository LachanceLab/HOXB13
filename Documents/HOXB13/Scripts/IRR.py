# -*- coding: utf-8 -*-
# Plot IRR with CIs by HOXB13 carrier status and RFX6 genotype (CC/TC/TT)

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- Input ----------
FILE_PATH = "IRR_AllofUsandMADCaPGhanaUganda.csv"  # Change this to your filename

# ---------- Load ----------
p = Path(FILE_PATH)
if not p.exists():
    raise FileNotFoundError(f"Couldn't find {FILE_PATH}. Update FILE_PATH to your data file.")

sep = "\t" if p.suffix.lower() in [".tsv", ".txt"] else ","
df = pd.read_csv(p, sep=sep)

# Split genotype info
parts = df["Genotype"].str.split("_", n=1, expand=True)
df["carrier_raw"] = parts[0]  # Non or Car
df["dosage"] = parts[1].astype(int)

carrier_map = {"Non": "Non-carriers", "Car": "X285K Carriers"}
df["carrier"] = df["carrier_raw"].map(carrier_map)

colors = {"Non-carriers": "tab:blue", "X285K Carriers": "tab:red"}

# Convert numeric
for col in ["IRR", "IRR_low", "IRR_high"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ---------- X positions ----------
order_carriers = ["Non-carriers", "X285K Carriers"]
order_dosage = [0, 1, 2]
genotype_labels = {0: "CC", 1: "CT", 2: "TT"}

x_positions = {}
x = 0
gap = 1.0  # Reduced gap from 1.5 to 1.0 for less separation
for c in order_carriers:
    for d in order_dosage:
        x_positions[(c, d)] = x
        x += 1
    x += gap

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(6.5, 5))  # Reduced width from 7.5 to 6.5

# Reference line (Non-carriers CC)
ax.axhline(1.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

handles = {}
for c in order_carriers:
    sub = df[df["carrier"] == c].set_index("dosage")
    for d in order_dosage:
        if d not in sub.index:
            continue
        row = sub.loc[d]
        xval = x_positions[(c, d)]
        y = row["IRR"]
        ylow, yhigh = row["IRR_low"], row["IRR_high"]
        if pd.notna(ylow) and pd.notna(yhigh):
            yerr = [[y - ylow], [yhigh - y]]
        else:
            yerr = None
        h = ax.errorbar(
            xval,
            y,
            yerr=yerr,
            fmt="o",
            capsize=4,
            elinewidth=1.5,
            markersize=6,
            color=colors[c],
            ecolor=colors[c],
        )
        handles[c] = h

# ---------- Axis formatting ----------
xticks = []
xticklabels = []
for c in order_carriers:
    for d in order_dosage:
        xticks.append(x_positions[(c, d)])
        xticklabels.append(genotype_labels[d])
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, fontsize=14) 

# Carrier group labels will be placed at the top of the plot after y-limits are set
# (See placement after setting y-axis limits)

ax.set_ylabel("Incidence Rate Ratio (IRR)", fontsize=16)
ax.set_xlabel("RFX6 genotype", fontsize=14)
# ax.set_title("IRRs by HOXB13 X285K Carrier Status and RFX6 Genotype")

# Remove legend since we're using x-axis labels instead

# Make y-axis tick labels bigger
ax.tick_params(axis='y', labelsize=14)

# Fixed y-axis range
# ax.set_ylim(0, 5)
# ax.set_ylim(0, 8)
# ax.set_ylim(0, 11)
ax.set_ylim(0, 7)
# Place carrier group labels near the top of the plot (5% below the top edge)
ymin, ymax = ax.get_ylim()
y_range = ymax - ymin
label_y_position = ymax - 0.05 * y_range
for c in order_carriers:
    x_center = sum(x_positions[(c, d)] for d in order_dosage) / len(order_dosage)
    ax.text(
        x_center,
        label_y_position,
        c,
        ha="center",
        va="top",
        fontsize=14,
        color=colors[c]
    )
ax.grid(axis="y", alpha=0.2)

plt.tight_layout()
plt.savefig("irr_plot.png", dpi=300, bbox_inches='tight')
# Also save as PDF
plt.savefig("irr_plot.pdf", dpi=300, bbox_inches='tight')