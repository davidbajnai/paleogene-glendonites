# The following code is used to:
# Compare dual clumped isotope measurements with equilibrium values
#
# INPUT:    dual_clumped_data.csv
#
# OUTPUT:   dual_clumped_plot.png
#           dual_clumped_plot.svg


# Import libraries -------------------------------------------------------------
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from _helper_functions import *
from _plot_style import *
from pathlib import Path

# Define directories -----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
source_data_dir = Path(__file__).resolve().parent.parent / "data" / "input"

# Apply plot style -------------------------------------------------------------
apply_plot_style()
plt.rcParams["figure.figsize"] = (3.5, 3.5)

# Load data --------------------------------------------------------------------
data = pd.read_csv(source_data_dir / "dual_clumped_data.csv")

# Define functions -------------------------------------------------------------
def fiebig_2024(t_in_c: ArrayLike) -> tuple[ArrayLike, ArrayLike]:
    t = t_in_c + 273.15
    d63 = -5.897 * 1/t + -3.521 * 10**3 / t**2 + \
        2.391 * 10**7 / t**3 - 3.541 * 10**9 / t**4
    d64 = 6.002 * 1/t - 1.299 * 10**4 / t**2 + \
        8.996 * 10**6 / t**3 - 7.423 * 10**8 / t**4
    d47 = 1.038* d63 + 0.1848
    d48 = 1.038* d64 + 0.1214
    return d47, d48

# Define colors ----------------------------------------------------------------
light_blue = C["cyan5"]
dark_gray = C["cgray7"]
light_gray = C["cgray3"]

sample_styles = {
    "1567C": {"fc": light_blue, "ec": light_blue, "group": "Age"},
    "1568A": {"fc": light_blue, "ec": C["violet5"], "group": "Age"},
    "1570A": {"fc": light_blue, "ec": "black", "group": "Age"},
    "+15": {"fc": C["red"], "ec": C["red"], "group": "Other"},
    "Congo fan": {"fc": dark_gray, "ec": dark_gray, "group": "Site"},
}

legend_handles = []
legend_sections = {
    "Age": "Exp. 396 glendonites",
    "Site": "modern ikaites",
    "Other": "Fur Fm. glendonites",
}
current_legend_section = None

# Create figure ----------------------------------------------------------------
fig, ax = plt.subplots()

# Setup axes -------------------------------------------------------------------
ax.set_xlabel(r"$\Delta_{48}$ (‰ CDES90)")
ax.set_ylabel(r"$\Delta_{47}$ (‰ CDES90)")
ax.set_ylim(0.617, 0.673)
ax.set_xlim(0.235, 0.305)

# Plot data --------------------------------------------------------------------

# Plot equilibrium line
t_min = 0
t_max = 15
t_range = np.linspace(t_min, t_max, 1000)
d47_eq, d48_eq = fiebig_2024(t_range)
ax.plot(d48_eq, d47_eq, color="black", ls=":", lw=0.5)

t_markers = np.arange(t_min, t_max + 1, 5)
d47_markers, d48_markers = fiebig_2024(t_markers)
ax.scatter(
    d48_markers,
    d47_markers,
    marker=".",
    facecolors="none",
    edgecolors="black",
)
for t, d47, d48 in zip(t_markers, d47_markers, d48_markers):
    ax.annotate(
        f"{t}°C",
        xy=(d48, d47),
        xytext=(2, 0),
        textcoords="offset points",
        fontsize=6,
        va="top",
        ha="left",
        zorder = 10
    )

# Plot data points with error bars
for sample_name, style in sample_styles.items():
    sample_data = data[data["SampleName"] == sample_name]
    if sample_data.empty:
        continue

    age = sample_data["Age"].iloc[0]
    if style["group"] == "Age":
        legend_label = f"{age} ({sample_name})"
    elif style["group"] == "Site":
        legend_label = sample_name
    else:
        legend_label = f"{sample_name} ({age})"

    if style["group"] != current_legend_section:
        legend_handles.append(
            Line2D([], [], color="none", linestyle="none",
                   label=legend_sections[style["group"]])
        )
        current_legend_section = style["group"]

    ax.errorbar(
        sample_data["D48"],
        sample_data["D47"],
        yerr=sample_data["2SE_D47"],
        xerr=sample_data["2SE_D48"],
        fmt="D",
        markerfacecolor=style["fc"],
        markeredgecolor=style["ec"],
        markeredgewidth=0.8,
        ecolor=style["fc"] if style["group"] == "Age" else style["ec"],
        markersize=4,
        capsize=0,
        linestyle="none",
        label=legend_label,
    )
    legend_handles.append(
        Line2D(
            [],
            [],
            color="none",
            marker="D",
            markerfacecolor=style["fc"],
            markeredgecolor=style["ec"],
            markeredgewidth=0.8,
            markersize=4,
            linestyle="none",
            label=legend_label,
        )
    )


# Add legend -------------------------------------------------------------------
legend = ax.legend(handles=legend_handles, loc="lower right",
                   fancybox=False,
                   edgecolor="k",
                   facecolor="none",
                   markerscale=0.9,
                   fontsize=6,
                   borderpad=0.5, labelspacing=0.25,
                   handlelength=1.0, handletextpad=0.3,
                   borderaxespad=0.4
                   )
for text in legend.get_texts():
    if text.get_text() in legend_sections.values():
        text.set_weight("bold")

save_figure(fig, figures_dir,  "dual_clumped_plot")
