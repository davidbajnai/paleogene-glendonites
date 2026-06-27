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
from reconstruct_water_Dp17O import reconstruct_water_Dp17O
from waters import calculate_water_endmembers


# Define directories ----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
export_tables_dir = Path(__file__).resolve().parent.parent / "export_tables" / "csv"

d18O_ice_free, Dp17O_ice_free, d18O_ice, Dp17O_ice, d18O_modern_seawater, Dp17O_modern_seawater, d18O_meteoric_water, Dp17O_meteoric_water = calculate_water_endmembers()


# Apply plot style ------------------------------------------------------------
apply_plot_style()


# Define functions ------------------------------------------------------------
def plot_calcite_equilibrium(Dp17Ow, d18Ow, Tmin, Tmax, ax, color="k"):
    d17Ow = d17O(d18Ow, Dp17Ow)

    # Plot the water point
    ax.scatter(prime(d18Ow), Dp17Ow,
               marker='*', ec=color, fc="w", s=60, zorder=10)

    # Plot the equilibrium line
    equilibrium_temperatures = np.arange(Tmin, 290, 1) + 273.15
    d18O_mineral = d18O_cc(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17O_cc(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral),
         equilibrium_temperatures]
    ).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            c=color, ls=":", zorder=0)

    # Plot the equilibrium line, marking every 5 °C
    equilibrium_temperatures = np.arange(Tmin, Tmax + 1, 5) + 273.15
    d18O_mineral = d18O_cc(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17O_cc(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral),
         equilibrium_temperatures]
    ).T
    ax.scatter(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
               s=15, marker=".", fc="white", ec=color, zorder=0)

    # Return equilibrium data as a dataframe
    equilibrium_df = pd.DataFrame(mineral_equilibrium)
    equilibrium_df[2] = equilibrium_df[2] - 273.15
    equilibrium_df = equilibrium_df.rename(
        columns={0: "d18O", 1: "Dp17O", 2: "temperature"}
    )
    return equilibrium_df


# Import data -----------------------------------------------------------------
df_carbonates = pd.read_csv(export_tables_dir / "combined_replicate_data.csv")

reconstruction_temperature_C = 4
reconstruction_temperature_error_C = 4
non_modern = df_carbonates["Age"] != "modern"

df_carbonates["temperature"] = np.nan
df_carbonates["temperature_error"] = np.nan
df_carbonates["reconstructed_water_d18O"] = np.nan
df_carbonates["reconstructed_water_Dp17O"] = np.nan

df_carbonates.loc[non_modern, "temperature"] = reconstruction_temperature_C
df_carbonates.loc[non_modern, "temperature_error"] = reconstruction_temperature_error_C

rectangle_columns = [
    "reconstructed_water_d18O_prime_min",
    "reconstructed_water_d18O_prime_max",
    "reconstructed_water_Dp17O_min",
    "reconstructed_water_Dp17O_max",
]
df_carbonates[rectangle_columns] = np.nan
if non_modern.any():
    reconstructed = reconstruct_water_Dp17O(
        df_carbonates.loc[non_modern, "d18O_AC"],
        df_carbonates.loc[non_modern, "Dp17O_AC"],
        df_carbonates.loc[non_modern, "temperature"],
        d18O_error=df_carbonates.loc[non_modern, "d18O_error"],
        Dp17O_error=df_carbonates.loc[non_modern, "Dp17O_error"],
        temperature_error_C=df_carbonates.loc[non_modern, "temperature_error"],
    )
    reconstructed_columns = [
        "reconstructed_water_d18O",
        "reconstructed_water_Dp17O",
    ] + rectangle_columns
    df_carbonates.loc[non_modern, reconstructed_columns] = reconstructed[
        reconstructed_columns
    ]

export_table(df_carbonates, "all_data")

# Plotting --------------------------------------------------------------------
fig, ax = plt.subplots()

light_blue = "#25AAE1"
dark_gray = "#231F20"
light_gray = "#B3B4B4"

sample_styles = {
    "pre-PETM": {"fc": light_blue, "ec": light_blue, "group": "Age"},
    "PETM": {"fc": light_blue, "ec": C["red"], "group": "Age"},
    "early Eocene": {"fc": light_blue, "ec": "black", "group": "Age"},
    "Congo fan": {"fc": dark_gray, "ec": dark_gray, "group": "Site"},
    "South Georgia": {"fc": light_gray, "ec": dark_gray, "group": "Site"},
    "Bransfield Strait": {"fc": dark_gray, "ec": light_gray, "group": "Site"},
    "Nankai Trough": {"fc": "white", "ec": "black", "group": "Site"},
}

# Plot modern seawater --------------------------------------------------------

ax.text(prime(d18O_modern_seawater)+2, Dp17O_modern_seawater,
        "modern\nseawater",
        va="bottom", ha="left", color="#858379", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
# Plot calcite line for modern seawater
eq = plot_calcite_equilibrium(
    Dp17O_modern_seawater, d18O_modern_seawater, 0, 15, ax, "#858379")
ax.annotate(
    str(round(eq.iloc[0, 2])) + " °C",
    (prime(eq.iloc[0, 0]), eq.iloc[0, 1]),
    xytext=(prime(eq.iloc[0, 0]) + 3.5, eq.iloc[0, 1] + 8),
    ha="left", va="center", fontsize=6,
    # bbox=dict(fc="white", ec="none", pad=1),
    arrowprops=dict(arrowstyle="->", color="k", shrinkA=0))
ax.annotate(
    str(round(eq.iloc[-1, 2])) + " °C",
    (prime(eq.iloc[-1, 0]), eq.iloc[-1, 1]),
    xytext=(prime(eq.iloc[-1, 0]) + 3.5, eq.iloc[-1, 1] + 8),
    ha="left", va="center", fontsize=6,
    # bbox=dict(fc="white", ec="none", pad=1),
    arrowprops=dict(arrowstyle="->", color="k", shrinkA=0))

# Plot ice-free seawater composition ------------------------------------------

ax.scatter(prime(d18O_ice_free), Dp17O_ice_free,
           marker='*', ec="#1455C0", fc="w", s=60, zorder = 10)
ax.text(prime(d18O_ice_free) + 2, Dp17O_ice_free,
        "ice-free\nseawater",
        va="top", ha="left", color="#1455C0", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
eq = plot_calcite_equilibrium(
    Dp17O_ice_free, d18O_ice_free, 0, 15, ax, "#1455C0")

# Plot modeled fluids -----------------

ax.text(-9, -35,
        "reconstructed\nparent waters",
        va="center", ha="right", color="k", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )

# Plot samples
order = ['pre-PETM', 'PETM', 'early Eocene', 'modern']
order_map = {key: i for i, key in enumerate(order)}
df_carbonates = df_carbonates.sort_values(
    by='Age',
    key=lambda x: x.map(lambda val: next((order_map[k] for k in order if str(val).startswith(k)), float('inf')))
)

legend_handles = []
legend_sections = {
    "Age": "Exp 396 glendonites",
    "Site": "modern ikaites",
}
current_legend_section = None

df_reconstructed = df_carbonates[df_carbonates["Age"] != "modern"].copy()
if not df_reconstructed.empty:
    for _, row in df_reconstructed.iterrows():
        x_min = row["reconstructed_water_d18O_prime_min"]
        x_max = row["reconstructed_water_d18O_prime_max"]
        y_min = row["reconstructed_water_Dp17O_min"]
        y_max = row["reconstructed_water_Dp17O_max"]
        rectangle_color = sample_styles[row["Age"]]["ec"]
        ax.add_patch(Rectangle(
            (x_min, y_min),
            x_max - x_min,
            y_max - y_min,
            fc=mcolors.to_rgba(rectangle_color, 0.3),
            ec="none",
            zorder=2,
        ))

    age_mean_handles = []
    for age in order:
        if age == "modern":
            continue
        age_dat = df_reconstructed[df_reconstructed["Age"] == age]
        if age_dat.empty:
            continue

        style = sample_styles[age]
        mean_d18O_water = age_dat["reconstructed_water_d18O"].mean()
        mean_Dp17O_water = age_dat["reconstructed_water_Dp17O"].mean()
        ax.scatter(
            prime(mean_d18O_water),
            mean_Dp17O_water,
            marker="*",
            s=60,
            fc=style["fc"],
            ec=style["ec"],
            zorder=14,
        )
        print(f"Modeled parent water for {age}: d18O = {mean_d18O_water:.1f} ‰, Dp17O = {mean_Dp17O_water:.0f} ppm")

for style_key, style in sample_styles.items():
    dat = df_carbonates[df_carbonates[style["group"]] == style_key]
    if dat.empty:
        continue

    if style["group"] != current_legend_section:
        legend_handles.append(
            Line2D([], [], color="none", linestyle="none",
                   label=legend_sections[style["group"]])
        )
        current_legend_section = style["group"]

    x_vals = prime(dat["d18O_AC"])
    y_vals = dat["Dp17O_AC"]
    marker_label = f"{style_key}"

    sample_handle = ax.scatter(
        x_vals, y_vals,
        marker="D",
        s=20,
        fc=style["fc"],
        ec=style["ec"],
        zorder=10,
        lw=0.8,
        label=marker_label
    )
    legend_handles.append(sample_handle)

    ax.errorbar(
        x_vals, y_vals,
        xerr=dat["d18O_error"], yerr=dat["Dp17O_error"],
        fmt="none",
        color=style["ec"],
        zorder=5
    )

# Add meteoric water to the mix -----------------------------------------------

df_mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=d18O_meteoric_water, D17O_B=Dp17O_meteoric_water)
ax.plot(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"],
        color="#3CB5AE", ls="-", lw=1.5)
ax.scatter(prime(df_mix["mix_d18O"].iloc[-1]), df_mix["mix_Dp17O"].iloc[[-1]],
            marker = "*", ec="#3CB5AE", fc="w", s=60, zorder=10)
ax.text(prime(df_mix["mix_d18O"].iloc[-1]) + 2, df_mix["mix_Dp17O"].iloc[-1],
        "modeled\nmeteoric water",
        va="center", ha="left", color="#3CB5AE", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )

# Add ice to the mix ----------------------------------------------------------

# Mixing line between modern seawater and Greenland ice
df_mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=d18O_ice, D17O_B=Dp17O_ice)
ax.plot(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"], color="#55B9E6", linestyle="-", lw=1.5)
ax.scatter(prime(df_mix["mix_d18O"].iloc[-1]), df_mix["mix_Dp17O"].iloc[[-1]],
            marker = "*", ec="#55B9E6", fc="w", s=60, zorder=100)
ax.text(prime(d18O_ice) + 2, Dp17O_ice,
        "modeled\nsnowmelt",
        va="center", ha="left", color="#55B9E6", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )


df_mix_meteoric_ice = mix_d17O(d18O_A=d18O_meteoric_water, D17O_A=Dp17O_meteoric_water,
                               d18O_B=d18O_ice, D17O_B=Dp17O_ice, step=10)

# Add back the original endmember line (meteoric ↔ snowmelt)
# while preserving the current gradient-colored fan.
ax.plot(
    prime(df_mix_meteoric_ice["mix_d18O"]),
    df_mix_meteoric_ice["mix_Dp17O"],
    color="#49B7C8",
    linestyle="--",
    lw=1.2,
    zorder=1,
)

# for every point in df_mix, draw mixing line between the endmember and ice-free ocean
meteoric_rgb = np.array(mcolors.to_rgb("#3CB5AE"))
ice_rgb = np.array(mcolors.to_rgb("#55B9E6"))

for i in range(len(df_mix_meteoric_ice)):
    blend_ratio = i / (len(df_mix_meteoric_ice) - 1) if len(df_mix_meteoric_ice) > 1 else 0
    line_color = mcolors.to_hex((1 - blend_ratio) * meteoric_rgb + blend_ratio * ice_rgb)

    mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=df_mix_meteoric_ice["mix_d18O"].iloc[i], D17O_B=df_mix_meteoric_ice["mix_Dp17O"].iloc[i],
                  step=50)
    ax.scatter(prime(mix["mix_d18O"]), mix["mix_Dp17O"], fc=line_color, ec ="none", zorder=0, s=1, marker=".")


# Axis and legend settings ----------------------------------------------------
legend = ax.legend(handles=legend_handles, loc='lower left',
                   fancybox=False, markerscale=0.9,
                   framealpha=1, fontsize=6, edgecolor="k", facecolor="none")
for text in legend.get_texts():
    if text.get_text() in legend_sections.values():
        text.set_weight("bold")
ax.set_xlim(-35, 45)
ax.set_ylim(-127, 52)
ax.set_ylabel(r"$\Delta\prime^{17}$O (ppm)")
ax.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW)")

# set minor ticks
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(10))

# Save figure
save_figure(fig, figures_dir, "glendonites")

# Manuscript parts:

print(f'\nModern glendonites: ∆17O = {df_carbonates[df_carbonates["Age"] == "modern"]["Dp17O_AC"].min():.0f} to {df_carbonates[df_carbonates["Age"] == "modern"]["Dp17O_AC"].max():.0f} ppm and d18O = {df_carbonates[df_carbonates["Age"] == "modern"]["d18O_AC"].min():.1f} to {df_carbonates[df_carbonates["Age"] == "modern"]["d18O_AC"].max():.1f} ‰. Count = {len(df_carbonates[df_carbonates["Age"] == "modern"])}')
print(f'pre-PETM glendonites: ∆17O = {df_carbonates[df_carbonates["Age"] == "pre-PETM"]["Dp17O_AC"].min():.0f} to {df_carbonates[df_carbonates["Age"] == "pre-PETM"]["Dp17O_AC"].max():.0f} ppm and d18O = {df_carbonates[df_carbonates["Age"] == "pre-PETM"]["d18O_AC"].min():.1f} to {df_carbonates[df_carbonates["Age"] == "pre-PETM"]["d18O_AC"].max():.1f} ‰. Count = {len(df_carbonates[df_carbonates["Age"] == "pre-PETM"])}')
print(f'PETM glendonites: ∆17O = {df_carbonates[df_carbonates["Age"] == "PETM"]["Dp17O_AC"].min():.0f} to {df_carbonates[df_carbonates["Age"] == "PETM"]["Dp17O_AC"].max():.0f} ppm and d18O = {df_carbonates[df_carbonates["Age"] == "PETM"]["d18O_AC"].min():.1f} to {df_carbonates[df_carbonates["Age"] == "PETM"]["d18O_AC"].max():.1f} ‰. Count = {len(df_carbonates[df_carbonates["Age"] == "PETM"])}')
print(f'early Eocene glendonites: ∆17O = {df_carbonates[df_carbonates["Age"] == "early Eocene"]["Dp17O_AC"].min():.0f} to {df_carbonates[df_carbonates["Age"] == "early Eocene"]["Dp17O_AC"].max():.0f} ppm and d18O = {df_carbonates[df_carbonates["Age"] == "early Eocene"]["d18O_AC"].min():.1f} to {df_carbonates[df_carbonates["Age"] == "early Eocene"]["d18O_AC"].max():.1f} ‰. Count = {len(df_carbonates[df_carbonates["Age"] == "early Eocene"])}')
