# Import libraries ------------------------------------------------------------
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from _helper_functions import *
from _plot_style import *
from pathlib import Path


# Define directories ----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
source_data_dir = Path(__file__).resolve().parent.parent / "source_data"


# Apply plot style ------------------------------------------------------------
apply_plot_style()
plt.rcParams["figure.figsize"] = (3.5, 3.5)

# Load data -------------------------------------------------------------------
data = pd.read_csv(source_data_dir / "dual_clumped_data.csv")


# Define functions ------------------------------------------------------------
def Fiebig_2024(T_in_C):
    T = T_in_C + 273.15
    D63 = -5.897 * 1/T + -3.521 * 10**3 / T**2 + \
        2.391 * 10**7 / T**3 - 3.541 * 10**9 / T**4
    D64 = 6.002 * 1/T - 1.299 * 10**4 / T**2 + \
        8.996 * 10**6 / T**3 - 7.423 * 10**8 / T**4
    D47 = 1.038* D63 + 0.1848
    D48 = 1.038* D64 + 0.1214
    return D47, D48

# Define colors ---------------------------------------------------------------

light_blue = "#25AAE1"
dark_gray = "#231F20"
light_gray = "#B3B4B4"

sample_styles = {
    "1567C": {"fc": light_blue, "ec": light_blue, "group": "Age"},
    "1568A": {"fc": light_blue, "ec": C["red"], "group": "Age"},
    "1570A": {"fc": light_blue, "ec": "black", "group": "Age"},
    "+15": {"fc": C["red"], "ec": C["red"], "group": "Other"},
    "Congo fan": {"fc": dark_gray, "ec": dark_gray, "group": "Site"},
}

legend_handles = []
legend_sections = {
    "Age": "Exp 396 glendonites",
    "Site": "modern ikaites",
    "Other": "Fur Fm. glendonites",
}
current_legend_section = None

# Create figure and axes ------------------------------------------------------
fig, ax = plt.subplots()

# Plot equilibrium line
Tmin = 0
Tmax = 15
T_range = np.linspace(Tmin, Tmax, 1000)
D47_eq, D48_eq = Fiebig_2024(T_range)
ax.plot(D48_eq, D47_eq, color="black", ls=":", lw=0.5)

T_markers = np.arange(Tmin, Tmax + 1, 5)
D47_markers, D48_markers = Fiebig_2024(T_markers)
ax.scatter(
    D48_markers,
    D47_markers,
    marker=".",
    facecolors="none",
    edgecolors="black",
)
for T, D47, D48 in zip(T_markers, D47_markers, D48_markers):
    ax.annotate(
        f"{T}°C",
        xy=(D48, D47),
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
        ecolor=style["ec"],
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

ax.set_xlabel(r"$\Delta_{48}$ (‰ CDES90)")
ax.set_ylabel(r"$\Delta_{47}$ (‰ CDES90)")

ax.set_ylim(0.617, 0.673)
ax.set_xlim(0.235, 0.305)

legend = ax.legend(handles=legend_handles, loc="lower right",
                   fancybox=False, markerscale=1,
                   fontsize=6, edgecolor="k", facecolor="none",
                   borderpad=0.5, labelspacing=0.5,
                   handlelength=1.0, handletextpad=0.4,
                   borderaxespad=0.4)
for text in legend.get_texts():
    if text.get_text() in legend_sections.values():
        text.set_weight("bold")

save_figure(fig, figures_dir,  "dual_clumped_plot")
