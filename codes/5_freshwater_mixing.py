# Import libraries ------------------------------------------------------------
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import gridspec
from _helper_functions import *
from _plot_style import *


# Apply plot style ------------------------------------------------------------
apply_plot_style()


# Defne directory paths ------------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
data_dir = Path(__file__).resolve().parent.parent / "export_tables" / "csv"


# Additional functions ----------------------------------------------------------------
def mix_properties(
    freshwater_fraction,
    sw_salinity=35.0,
    fw_salinity=0.0,
    sw_d18o=-1.0,
    fw_d18o=-15.0
):
    """
    Compute salinity and d18O for binary mixing between seawater and freshwater.

    Parameters
    ----------
    freshwater_fraction : np.ndarray
        Fraction of freshwater in the mixture, from 0 to 1.
    sw_salinity : float
        Salinity of seawater.
    fw_salinity : float
        Salinity of freshwater.
    sw_d18o : float
        d18O value of seawater in per mil.
    fw_d18o : float
        d18O value of freshwater in per mil.

    Returns
    -------
    salinity : np.ndarray
        Mixed salinity values.
    d18o : np.ndarray
        Mixed d18O values.
    """
    seawater_fraction = 1.0 - freshwater_fraction

    salinity = seawater_fraction * sw_salinity + freshwater_fraction * fw_salinity
    d18o = seawater_fraction * sw_d18o + freshwater_fraction * fw_d18o

    return salinity, d18o


def d18o_to_salinity(d18o, sw_salinity=35.0, fw_salinity=0.0, sw_d18o=-1.0, fw_d18o=-15.0):
    """Convert d18O to salinity using a linear two-endmember mixing model."""
    freshwater_fraction = (d18o - sw_d18o) / (fw_d18o - sw_d18o)
    return (1.0 - freshwater_fraction) * sw_salinity + freshwater_fraction * fw_salinity


def salinity_to_d18o(salinity, sw_salinity=35.0, fw_salinity=0.0, sw_d18o=-1.0, fw_d18o=-15.0):
    """Inverse transform for secondary y-axis: salinity to d18O."""
    freshwater_fraction = (salinity - sw_salinity) / (fw_salinity - sw_salinity)
    return (1.0 - freshwater_fraction) * sw_d18o + freshwater_fraction * fw_d18o


sw_salinity = 35.0
sw_d18o = -1.1

# Freshwater endmember is itself a mix of ice melt and meteoric water.
fw_salinity = 0.0
d18o_ice = -32.25
d18o_meteoric = -16.84
ice_fraction_in_fw = np.arange(0.0, 1.0 + 0.1, 0.1)
fw_d18o = (
    ice_fraction_in_fw * d18o_ice
    + (1.0 - ice_fraction_in_fw) * d18o_meteoric
)

fig = plt.figure(figsize=(4.8, 3.6))
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[20, 1], wspace=0.35)
ax = fig.add_subplot(gs[0, 0])
cax = fig.add_subplot(gs[0, 1])

# Model a continuous range of mixed-water d18O values.
modeled_d18o_values = np.arange(-1, -10.0, -0.01)
cmap = plt.get_cmap("cividis_r")
norm = plt.Normalize(vmin=modeled_d18o_values.min(), vmax=modeled_d18o_values.max())

for d18o_value in modeled_d18o_values:
    freshwater_fraction_in_seawater = (d18o_value - sw_d18o) / (fw_d18o - sw_d18o)
    freshening_pct = freshwater_fraction_in_seawater * 100.0

    # Keep only physically meaningful freshening values.
    mask = (freshening_pct >= 0.0) & (freshening_pct <= 100.0)
    x_curve = ice_fraction_in_fw[mask] * 100.0
    salinity_curve = sw_salinity * (1.0 - freshwater_fraction_in_seawater[mask])

    if np.any(mask):
        ax.plot(
            x_curve,
            salinity_curve,
            color=cmap(norm(d18o_value)),
            lw=1.8,
            linestyle="-",
        )

# Add the modeld glendonite parent woter points for reference.
modeled_points = {
    "early Eocene": -4.6,
    "PETM": -7.2,
    "pre-PETM": -4.9,
}

colors = {
    "early Eocene": "k",
    "PETM": C["red"],
    "pre-PETM": "#25AAE1",
}

for label, d18o_value in modeled_points.items():
    freshwater_fraction_in_seawater = (d18o_value - sw_d18o) / (fw_d18o - sw_d18o)
    freshening_pct = freshwater_fraction_in_seawater * 100.0

    # Keep only physically meaningful freshening values.
    mask = (freshening_pct >= 0.0) & (freshening_pct <= 100.0)
    if np.any(mask):
        x_curve = ice_fraction_in_fw[mask] * 100.0
        salinity_value = sw_salinity * (1.0 - freshwater_fraction_in_seawater[mask])
        ax.plot(
            x_curve,
            salinity_value,
            color=colors[label],
            label=f"{label}",
            lw = 1.5,
            zorder=2,
        )

ax.text(50, 15, "tolerance limit for euryhaline taxa",  va="bottom", ha="center")
ax.hlines(15, 0, 100, color="black", lw=1.5, linestyle="--", zorder=5)

ax.text(50, 25, "tolerance limit for stenohaline taxa",  va="bottom", ha="center")
ax.hlines(25, 0, 100, color="black", lw=1.5, linestyle="--", zorder=5)

add_standard_legend(ax, loc="lower right",
                    bbox_to_anchor=(1.0, 0.1))

ax.set_xlabel("Snowmelt contribution to freshwater endmember (%)")
ax.set_ylabel("Salinity (g/kg)")
ax.set_xlim(0, 100)
ax.set_ylim(14, sw_salinity)
ax.set_xticks(np.arange(0, 101, 10))

# # Use a compact colorbar instead of a large legend for sequential curves.
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label(r"Modeled mixed-water $\delta^{18}O$ (‰)")

# Right axis is freshening (%), defined as freshwater contribution to seawater.
secax = ax.secondary_yaxis(
    "right",
    functions=(
        lambda sal: 100.0 * (1.0 - sal / sw_salinity),
        lambda fresh: sw_salinity * (1.0 - fresh / 100.0),
    ),
)
secax.set_ylabel("Freshwater contribution to seawater (%)")


save_figure(fig, figures_dir, "salinity", tight_layout=False)