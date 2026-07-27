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
sources_dir = Path(__file__).resolve().parent.parent
figures_dir = sources_dir / "figures"
export_tables_dir = sources_dir / "export_tables" / "csv"

# Get water endmembers ------------------------------------------------------
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

reconstruction_temperature_by_age_C = {
    "pre-PETM": 5,
    "PETM": 3,
    "early Eocene": 7,
}
reconstruction_temperature_error_by_age_C = {
    "pre-PETM": 3,
    "PETM": 2,
    "early Eocene": 2,
}
non_modern = df_carbonates["age"] != "modern"

df_carbonates["temperature"] = np.nan
df_carbonates["temperature_error"] = np.nan
df_carbonates["reconstructed_water_d18O"] = np.nan
df_carbonates["reconstructed_water_Dp17O"] = np.nan

df_carbonates["temperature"] = df_carbonates["age"].map(reconstruction_temperature_by_age_C)
df_carbonates["temperature_error"] = df_carbonates["age"].map(reconstruction_temperature_error_by_age_C)

missing_temperature_ages = df_carbonates.loc[
    non_modern
    & (
        df_carbonates["temperature"].isna()
        | df_carbonates["temperature_error"].isna()
    ),
    "age",
].unique()
if len(missing_temperature_ages) > 0:
    raise ValueError(
        "Missing predefined temperature values for ages: "
        + ", ".join(map(str, missing_temperature_ages))
    )

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
    "pre-PETM": {"fc": light_blue, "ec": light_blue, "group": "age"},
    "PETM": {"fc": light_blue, "ec": C["red"], "group": "age"},
    "early Eocene": {"fc": light_blue, "ec": "black", "group": "age"},
    "Congo fan": {"fc": dark_gray, "ec": dark_gray, "group": "site"},
    "South Georgia": {"fc": light_gray, "ec": dark_gray, "group": "site"},
    "Bransfield Strait": {"fc": dark_gray, "ec": light_gray, "group": "site"},
    "Nankai Trough": {"fc": "white", "ec": "black", "group": "site"},
}

# Plot modern seawater --------------------------------------------------------

ax.text(prime(d18O_modern_seawater)+2, Dp17O_modern_seawater,
        "modern\nseawater",
        va="bottom", ha="left", color="#858379", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
# Plot calcite line for modern seawater
eq = plot_calcite_equilibrium(
    Dp17O_modern_seawater, d18O_modern_seawater, 0, 10, ax, "#858379")
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
    Dp17O_ice_free, d18O_ice_free, 0, 10, ax, "#1455C0")

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
    by='age',
    key=lambda x: x.map(lambda val: next((order_map[k] for k in order if str(val).startswith(k)), float('inf')))
)

legend_handles = []
legend_sections = {
    "age": "Exp 396 glendonites",
    "site": "modern ikaites",
}
current_legend_section = None

df_reconstructed = df_carbonates[df_carbonates["age"] != "modern"].copy()
if not df_reconstructed.empty:
    for _, row in df_reconstructed.iterrows():
        x_min = row["reconstructed_water_d18O_prime_min"]
        x_max = row["reconstructed_water_d18O_prime_max"]
        y_min = row["reconstructed_water_Dp17O_min"]
        y_max = row["reconstructed_water_Dp17O_max"]
        rectangle_color = sample_styles[row["age"]]["ec"]
        ax.add_patch(Rectangle(
            (x_min, y_min),
            x_max - x_min,
            y_max - y_min,
            fc=mcolors.to_rgba(rectangle_color, 0.3),
            ec="none",
            zorder=2,
        ))

    age_mean_handles = []
    parent_water_d18O_uncertainties = []
    parent_water_Dp17O_uncertainties = []
    for age in order:
        if age == "modern":
            continue
        age_dat = df_reconstructed[df_reconstructed["age"] == age]
        if age_dat.empty:
            continue

        style = sample_styles[age]
        mean_d18O_water = age_dat["reconstructed_water_d18O"].mean()
        mean_Dp17O_water = age_dat["reconstructed_water_Dp17O"].mean()
        sample_d18O_uncertainties = (
            unprime(age_dat["reconstructed_water_d18O_prime_max"])
            - unprime(age_dat["reconstructed_water_d18O_prime_min"])
        ) / 2.0
        sample_Dp17O_uncertainties = (
            age_dat["reconstructed_water_Dp17O_max"]
            - age_dat["reconstructed_water_Dp17O_min"]
        ) / 2.0
        mean_d18O_uncertainty = sample_d18O_uncertainties.mean()
        mean_Dp17O_uncertainty = sample_Dp17O_uncertainties.mean()
        parent_water_d18O_uncertainties.extend(sample_d18O_uncertainties)
        parent_water_Dp17O_uncertainties.extend(sample_Dp17O_uncertainties)
        ax.scatter(
            prime(mean_d18O_water),
            mean_Dp17O_water,
            marker="*",
            s=60,
            fc=style["fc"],
            ec=style["ec"],
            zorder=14,
        )
        print(
            f"Modeled parent water for {age}: "
            f"δ18O = {mean_d18O_water:.1f} ‰, "
            f"Δ’17O = {mean_Dp17O_water:.0f} ppm"
        )

    print(
        "Average parent-water reconstruction uncertainty: "
        f"δ18O = {np.nanmean(parent_water_d18O_uncertainties):.1f} ‰, "
        f"Δ’17O = {np.nanmean(parent_water_Dp17O_uncertainties):.0f} ppm"
        "\n"
    )

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
        "Eocene\nmeteoric water",
        va="center", ha="left", color="#3CB5AE", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
df_mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=d18O_meteoric_water, D17O_B=Dp17O_meteoric_water, step=10)
ax.scatter(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)
ax.annotate("50%",
            (prime(df_mix["mix_d18O"].iloc[5]), df_mix["mix_Dp17O"].iloc[5]),
            xytext=(1, 10), textcoords="offset points",
            va="center", ha="left", color="k", fontsize=4,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.5,
                    connectionstyle="arc3,rad=-0.2", shrinkA=0),
            )

# Add ice to the mix ----------------------------------------------------------

# Mixing line between modern seawater and Greenland ice
df_mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=d18O_ice, D17O_B=Dp17O_ice)
ax.plot(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"], color="#55B9E6", linestyle="-", lw=1.5)
ax.scatter(prime(df_mix["mix_d18O"].iloc[-1]), df_mix["mix_Dp17O"].iloc[[-1]],
            marker = "*", ec="#55B9E6", fc="w", s=60, zorder=100)
ax.text(prime(d18O_ice) + 2, Dp17O_ice,
        "snowmelt",
        va="center", ha="left", color="#55B9E6", fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )



df_mix_meteoric_ice = mix_d17O(d18O_A=d18O_meteoric_water, D17O_A=Dp17O_meteoric_water,
                               d18O_B=d18O_ice, D17O_B=Dp17O_ice, step=100)
ax.plot(
    prime(df_mix_meteoric_ice["mix_d18O"]),
    df_mix_meteoric_ice["mix_Dp17O"],
    color="#49B7C8",
    linestyle="--",
    lw=1.2,
    zorder=1,
)
df_mix_meteoric_ice = mix_d17O(d18O_A=d18O_meteoric_water, D17O_A=Dp17O_meteoric_water,
                               d18O_B=d18O_ice, D17O_B=Dp17O_ice, step=4)
ax.scatter(prime(df_mix_meteoric_ice["mix_d18O"]), df_mix_meteoric_ice["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)
ax.annotate("50%",
            (prime(df_mix_meteoric_ice["mix_d18O"].iloc[2]), df_mix_meteoric_ice["mix_Dp17O"].iloc[2]),
            xytext=(1, 10), textcoords="offset points",
            va="center", ha="left", color="k", fontsize=4,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.5,
                    connectionstyle="arc3,rad=-0.2", shrinkA=0),
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
    mix = mix_d17O(d18O_A=d18O_ice_free, D17O_A=Dp17O_ice_free,
                  d18O_B=df_mix_meteoric_ice["mix_d18O"].iloc[i], D17O_B=df_mix_meteoric_ice["mix_Dp17O"].iloc[i],
                  step=10)
    ax.scatter(prime(mix["mix_d18O"]), mix["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)


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

print(f'\nModern glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "modern"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "modern"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "modern"])}')
print(f'pre-PETM glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "pre-PETM"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "pre-PETM"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "pre-PETM"])}')
print(f'PETM glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "PETM"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "PETM"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "PETM"])}')
print(f'early Eocene glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "early Eocene"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "early Eocene"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "early Eocene"])}')
