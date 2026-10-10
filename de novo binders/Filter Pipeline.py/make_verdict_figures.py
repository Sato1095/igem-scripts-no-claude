#!/usr/bin/env python3
"""Two figures for the Inhibitor Design verdict: (1) three-stage selection
funnel, (2) selection landscape (interface axes + chemistry)."""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ------------------------------------------------------------------ FIGURE 1
# Three-stage funnel
fig, ax = plt.subplots(figsize=(11, 7))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, fs=10, tc="black"):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle="round,pad=0.08", fc=fc, ec="#333", lw=1.2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc, wrap=True)

def arrow(x1, y1, x2, y2, c="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=16, lw=1.6, color=c))

# central funnel spine
box(3, 9.2, 4.6, 1.0, "2,000 candidates\nRFdiffusion × ProteinMPNN → AlphaFold Multimer", "#dfe7f0", 10)

box(3, 7.4, 4.6, 1.0, "STAGE 1 · Composite Z-score\n(AlphaFold confidence ranking)", "#cfe0f5", 10)
box(3, 5.6, 4.6, 1.0, "STAGE 2 · Confidence filter → 189 survivors\n+ Trp / anchor-register flags", "#cfe0f5", 10)
box(3, 3.4, 4.6, 1.3, "STAGE 3 · Interface + chemistry ranking\n1,275 fill 3 pockets → 409 pass gates\n→ 18 meet-or-exceed PMI", "#cfe0f5", 10)

arrow(3, 8.7, 3, 7.95); arrow(3, 6.9, 3, 6.1); arrow(3, 5.1, 3, 4.1)

# candidate branches on the right
box(8.0, 7.4, 3.4, 1.0, "42.3  (Z1)\nconfidence-optimal\nBSA/res 35.4 · 40 contacts", "#fde9c9", 9)
box(8.0, 5.6, 3.4, 1.0, "267.3  (Z2)\napparent mimic\ncharge −8", "#f6cccc", 9)
box(8.0, 3.4, 3.4, 1.0, "131.2  (Z4)\ninterface-optimal\nBSA/res 42.0 · 49 contacts", "#cdebd0", 9)

arrow(5.3, 7.4, 6.3, 7.4)
arrow(5.3, 5.6, 6.3, 5.6)
arrow(5.3, 3.4, 6.3, 3.4)

# 267 rejected marker
ax.text(8.0, 4.75, "✗ rejected at chemistry gate", ha="center", fontsize=8.5,
        color="#a11", style="italic")
ax.annotate("", xy=(8.0, 4.0), xytext=(8.0, 5.05),
            arrowprops=dict(arrowstyle="-|>", color="#a11", lw=1.3, ls=":"))

# final verdict
box(3, 1.3, 4.6, 1.0, "VERDICT · 131.2 leading candidate\nmeets PMI benchmark · MD-stable (100 ns)", "#cdebd0", 10)
arrow(3, 2.75, 3, 1.85)
arrow(6.3, 3.4, 5.3, 1.5, c="#2a7d38")

ax.text(3, 9.9, "Three-stage selection funnel", ha="center", fontsize=14, weight="bold")
fig.tight_layout()
fig.savefig("figure1_pipeline_funnel.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("figure1 saved")

# ------------------------------------------------------------------ FIGURE 2
df = pd.read_csv("/mnt/user-data/uploads/structural_ranking_full.csv")
df["key"] = df["design"].astype(str)+"."+df["n"].astype(str)
PMI_BSA, PMI_CON = 40.0, 35
picks = {"42.3":("#e6a13c","s","42.3 (Z1)"),
         "267.3":("#c0392b","D","267.3 (Z2)"),
         "131.2":("#2a7d38","*","131.2 (Z4)")}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios":[1.3,1]})

# panel 1 — interface axes
p3 = df[df["pockets_occupied"]==3]; pl = df[df["pockets_occupied"]<3]
ax1.scatter(pl["apolar_bsa_per_res"], pl["contacts"], s=7, c="#dcdcdc", alpha=.5, edgecolors="none", label="< 3 pockets")
ax1.scatter(p3["apolar_bsa_per_res"], p3["contacts"], s=9, c="#a9c9e8", alpha=.7, edgecolors="none", label="all 3 pockets")
xmax = df["apolar_bsa_per_res"].max()+1; ymax = df["contacts"].max()+3
ax1.fill_between([PMI_BSA, xmax], PMI_CON, ymax, color="#c7e9c0", alpha=.35, zorder=0)
ax1.axvline(PMI_BSA, color="#d6604d", ls="--", lw=1.5); ax1.axhline(PMI_CON, color="#d6604d", ls="--", lw=1.5)
ax1.scatter([PMI_BSA],[PMI_CON], marker="X", s=230, c="#d6604d", edgecolors="black", lw=.8, zorder=6)
ax1.annotate("PMI benchmark", (PMI_BSA, PMI_CON), textcoords="offset points", xytext=(-92,-2),
             fontsize=9.5, color="#a11", weight="bold", arrowprops=dict(arrowstyle="->", color="#a11"))
ax1.annotate("meets-or-exceeds PMI\n+ 3 pockets", (xmax-.3, ymax-1), ha="right", va="top",
             fontsize=10, color="#276221", weight="bold")
for k,(c,mk,lab) in picks.items():
    r = df[df["key"]==k].iloc[0]
    ax1.scatter(r["apolar_bsa_per_res"], r["contacts"], marker=mk, s=(320 if mk=="*" else 150),
                c=c, edgecolors="black", lw=1.1, zorder=8, label=lab)
ax1.set_xlabel("Apolar buried surface area per residue (Å²)", fontsize=12)
ax1.set_ylabel("Intermolecular contacts (≤ 5 Å)", fontsize=12)
ax1.set_title("Interface-quality landscape (n = 2000)", fontsize=13, weight="bold")
ax1.set_xlim(df["apolar_bsa_per_res"].min()-1, xmax); ax1.set_ylim(df["contacts"].min()-2, ymax)
ax1.legend(fontsize=8.5, loc="lower left"); ax1.grid(alpha=.15)

# panel 2 — chemistry gate
bins = np.arange(df["net_charge"].min()-.5, df["net_charge"].max()+1.5, 1)
ax2.hist(df["net_charge"], bins=bins, color="#cfd8e3", edgecolor="#7f8c9b", lw=.4)
ax2.axvline(-4, color="#333", ls=":", lw=1.6)
ax2.axvspan(df["net_charge"].min()-1, -4, color="#f4cccc", alpha=.35)
ax2.axvspan(-2.2, -0.8, color="#c7e9c0", alpha=.4)
ax2.annotate("chemistry gate\n(charge ≥ −4)", (-4, ax2.get_ylim()[1]*.8), fontsize=9.5, ha="right", color="#333")
ax2.annotate("PMI-like\nenvelope", (-1.5, ax2.get_ylim()[1]*.9), fontsize=9, ha="center", color="#276221", weight="bold")
for k,(c,mk,lab) in picks.items():
    r = df[df["key"]==k].iloc[0]
    ax2.axvline(r["net_charge"], color=c, lw=2.4)
    ax2.annotate(k, (r["net_charge"], ax2.get_ylim()[1]*(0.5-0.09*list(picks).index(k))),
                 ha="center", fontsize=9, color=c, weight="bold",
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=c, lw=.8))
ax2.set_xlabel("Net peptide charge (pH 7.4)", fontsize=12)
ax2.set_ylabel("Number of candidates", fontsize=12)
ax2.set_title("Chemistry gate — 267.3 rejected", fontsize=13, weight="bold")
ax2.grid(alpha=.15)

fig.suptitle("Selection landscape: where the three candidates sit", fontsize=14, weight="bold", y=1.00)
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig("figure2_selection_landscape.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("figure2 saved")
