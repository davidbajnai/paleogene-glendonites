# The following code is used to:
# Reconstruct parent-water compositions for Paleogene glendonites
#
# INPUT:    combined_replicate_data.csv
#           modern_seawater.csv
#           Aron_2021.csv
#           Pang_2015.csv
#           Touzeau_2016.csv
#           Landais_2008_S1.csv
#           Landais_2008_S2.csv
#           Steen_Larsen_2014.csv
#           Landais_2012.csv
#           Beershaw_2020.csv
#
# OUTPUT:   all_data.csv
#           all_data.xlsx
#           glendonites.png
#           glendonites.svg


# Import libraries -------------------------------------------------------------
from matplotlib.ticker import MultipleLocator
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from _helper_functions import *
from _plot_style import *
from pathlib import Path
from reconstruct_water_Dp17O import reconstruct_water_dp17o
from _water_functions import calculate_water_endmembers

# Define directories -----------------------------------------------------------
sources_dir = Path(__file__).resolve().parent.parent
figures_dir = sources_dir / "figures"
export_tables_dir = sources_dir / "data" / "exported" / "csv"

# Get water endmembers ---------------------------------------------------------
d18o_ice_free, dp17o_ice_free, d18o_ice, dp17o_ice, d18o_modern_seawater, dp17o_modern_seawater, d18o_meteoric_water, dp17o_meteoric_water = calculate_water_endmembers()

# Apply plot style -------------------------------------------------------------
apply_plot_style()

# Define functions -------------------------------------------------------------
def plot_calcite_equilibrium(
    dp17ow: float,
    d18ow: float,
    t_min: float,
    t_max: float,
    ax: Axes,
    color: str = "k",
) -> pd.DataFrame:
    d17ow = d17o(d18ow, dp17ow)

    # Plot the water point
    ax.scatter(prime(d18ow), dp17ow,
               marker='*', ec=color, fc="w", s=60, zorder=10)

    # Plot the equilibrium line
    equilibrium_temperatures = np.arange(t_min, 290, 1) + 273.15
    d18o_mineral = d18o_cc(equilibrium_temperatures, d18ow)
    d17o_mineral = d17o_cc(equilibrium_temperatures, d17ow)
    mineral_equilibrium = np.array(
        [d18o_mineral, dp17o(d17o_mineral, d18o_mineral),
         equilibrium_temperatures]
    ).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            c=color, ls=":", zorder=0)

    # Plot the equilibrium line, marking every 5 °C
    equilibrium_temperatures = np.arange(t_min, t_max + 1, 5) + 273.15
    d18o_mineral = d18o_cc(equilibrium_temperatures, d18ow)
    d17o_mineral = d17o_cc(equilibrium_temperatures, d17ow)
    mineral_equilibrium = np.array(
        [d18o_mineral, dp17o(d17o_mineral, d18o_mineral),
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

# Import data ------------------------------------------------------------------
df_carbonates = pd.read_csv(export_tables_dir / "combined_replicate_data.csv")

# Temperature assignment:
# Not all samples analyzed for Δ’17O have corresponding Δ47 measurements,
# while some additional samples were analyzed only for Δ47.
# Where available, sample-specific temperatures derived from Δ47 were used;
# otherwise, the mean Δ47-derived temperature for the corresponding age group
# was calculated from the manuscript's clumped isotope supplementary table.
reconstruction_temperature_by_age_c = {
    "pre-PETM": 5,
    "PETM": 3,
    "early Eocene": 7,
}
reconstruction_temperature_error_by_age_c = {
    "pre-PETM": 3,
    "PETM": 2,
    "early Eocene": 2,
}
non_modern = df_carbonates["age"] != "modern"

df_carbonates["reconstructed_water_d18O"] = np.nan
df_carbonates["reconstructed_water_Dp17O"] = np.nan

df_carbonates["clumped_T"] = df_carbonates["clumped_T"].fillna(
    df_carbonates["age"].map(reconstruction_temperature_by_age_c)
)
df_carbonates["clumped_T_err"] = df_carbonates["clumped_T_err"].fillna(
    df_carbonates["age"].map(reconstruction_temperature_error_by_age_c)
)

missing_temperature_ages = df_carbonates.loc[
    non_modern
    & (
        df_carbonates["clumped_T"].isna()
        | df_carbonates["clumped_T_err"].isna()
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
    reconstructed = reconstruct_water_dp17o(
        df_carbonates.loc[non_modern, "d18O_AC"],
        df_carbonates.loc[non_modern, "Dp17O_AC"],
        df_carbonates.loc[non_modern, "clumped_T"],
        d18o_error=df_carbonates.loc[non_modern, "d18O_error"],
        dp17o_error=df_carbonates.loc[non_modern, "Dp17O_error"],
        temperature_error_c=df_carbonates.loc[non_modern, "clumped_T_err"],
    )
    reconstructed_columns = [
        "reconstructed_water_d18O",
        "reconstructed_water_Dp17O",
    ] + rectangle_columns
    df_carbonates.loc[non_modern, reconstructed_columns] = reconstructed[
        reconstructed_columns
    ]

export_table(df_carbonates, "all_data")

# Create figure ----------------------------------------------------------------
fig, ax = plt.subplots()

# Setup axes -------------------------------------------------------------------
ax.set_xlim(-35, 45)
ax.set_ylim(-127, 52)
ax.set_ylabel(r"$\Delta\prime^{17}$O (ppm)")
ax.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW)")
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(10))

# Plot data --------------------------------------------------------------------

light_blue = C["sample_blue"]
dark_gray = C["sample_dark_gray"]
light_gray = C["sample_light_gray"]

sample_styles = {
    "pre-PETM": {"fc": light_blue, "ec": light_blue, "group": "age"},
    "PETM": {"fc": light_blue, "ec": C["red"], "group": "age"},
    "early Eocene": {"fc": light_blue, "ec": "black", "group": "age"},
    "Congo fan": {"fc": dark_gray, "ec": dark_gray, "group": "site"},
    "South Georgia": {"fc": light_gray, "ec": dark_gray, "group": "site"},
    "Bransfield Strait": {"fc": dark_gray, "ec": light_gray, "group": "site"},
    "Nankai Trough": {"fc": "white", "ec": "black", "group": "site"},
}

# Plot modern seawater ---------------------------------------------------------
ax.text(prime(d18o_modern_seawater)+2, dp17o_modern_seawater,
        "modern\nseawater",
        va="bottom", ha="left", color=C["wgray5"], fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
# Plot calcite line for modern seawater
eq = plot_calcite_equilibrium(
    dp17o_modern_seawater, d18o_modern_seawater, 0, 10, ax, C["wgray5"])
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

# Plot ice-free seawater composition -------------------------------------------
ax.scatter(prime(d18o_ice_free), dp17o_ice_free,
           marker='*', ec=C["blue"], fc="w", s=60, zorder = 10)
ax.text(prime(d18o_ice_free) + 2, dp17o_ice_free,
        "ice-free\nseawater",
        va="top", ha="left", color=C["blue"], fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
eq = plot_calcite_equilibrium(
    dp17o_ice_free, d18o_ice_free, 0, 10, ax, C["blue"])

# Plot modeled fluids ----------------------------------------------------------
ax.text(-11, -38,
        "reconstructed\nparent waters",
        va="center", ha="right", color="k", fontsize=6,
        )
ax.annotate("",
            (-11, -38),
            xytext=(10, 10),
            textcoords="offset points",
            va="center", ha="left", color="k", fontsize=4,
            arrowprops=dict(arrowstyle="<-", color="k", lw=0.5,
                    connectionstyle="arc3,rad=-0.2", shrinkA=0),
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
    parent_water_d18o_uncertainties = []
    parent_water_dp17o_uncertainties = []
    for age in order:
        if age == "modern":
            continue
        age_dat = df_reconstructed[df_reconstructed["age"] == age]
        if age_dat.empty:
            continue

        style = sample_styles[age]
        mean_d18o_water = age_dat["reconstructed_water_d18O"].mean()
        mean_dp17o_water = age_dat["reconstructed_water_Dp17O"].mean()
        sample_d18o_uncertainties = (
            unprime(age_dat["reconstructed_water_d18O_prime_max"])
            - unprime(age_dat["reconstructed_water_d18O_prime_min"])
        ) / 2.0
        sample_dp17o_uncertainties = (
            age_dat["reconstructed_water_Dp17O_max"]
            - age_dat["reconstructed_water_Dp17O_min"]
        ) / 2.0
        sample_count = len(age_dat)

        # Propagate the independent per-sample reconstruction errors into the
        # group mean and combine them with the standard error derived from the
        # between-sample 1 SD. For a single sample, the latter is undefined and
        # contributes zero additional variance.
        propagated_d18o_mean_error = (
            np.sqrt(np.sum(sample_d18o_uncertainties**2)) / sample_count
        )
        propagated_dp17o_mean_error = (
            np.sqrt(np.sum(sample_dp17o_uncertainties**2)) / sample_count
        )
        if sample_count > 1:
            d18o_standard_error = (
                age_dat["reconstructed_water_d18O"].std(ddof=1)
                / np.sqrt(sample_count)
            )
            dp17o_standard_error = (
                age_dat["reconstructed_water_Dp17O"].std(ddof=1)
                / np.sqrt(sample_count)
            )
        else:
            d18o_standard_error = 0.0
            dp17o_standard_error = 0.0

        mean_d18o_uncertainty = np.hypot(
            propagated_d18o_mean_error,
            d18o_standard_error,
        )
        mean_dp17o_uncertainty = np.hypot(
            propagated_dp17o_mean_error,
            dp17o_standard_error,
        )
        parent_water_d18o_uncertainties.append(mean_d18o_uncertainty)
        parent_water_dp17o_uncertainties.append(mean_dp17o_uncertainty)
        ax.scatter(
            prime(mean_d18o_water),
            mean_dp17o_water,
            marker="*",
            s=60,
            fc=style["fc"],
            ec=style["ec"],
            zorder=14,
        )
        print(
            f"Modeled parent water for {age}: "
            f"δ18O = {mean_d18o_water:.1f} ± {mean_d18o_uncertainty:.1f} ‰, "
            f"Δ’17O = {mean_dp17o_water:.0f} ± {mean_dp17o_uncertainty:.0f} ppm "
            "(propagated 1 SD)"
        )

    print(
        "Average parent-water reconstruction uncertainty: "
        f"δ18O = {np.nanmean(parent_water_d18o_uncertainties):.1f} ‰, "
        f"Δ’17O = {np.nanmean(parent_water_dp17o_uncertainties):.0f} ppm "
        "(mean propagated 1 SD)"
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

# Add meteoric water to the mix ------------------------------------------------
df_mix = mix_d17o(d18o_a=d18o_ice_free, dp17o_a=dp17o_ice_free,
                  d18o_b=d18o_meteoric_water, dp17o_b=dp17o_meteoric_water)
ax.plot(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"],
        color=C["meteoric"], ls="-", lw=1.5)
ax.scatter(prime(df_mix["mix_d18O"].iloc[-1]), df_mix["mix_Dp17O"].iloc[[-1]],
            marker = "*", ec=C["meteoric"], fc="w", s=60, zorder=10)
ax.text(prime(df_mix["mix_d18O"].iloc[-1]) + 2, df_mix["mix_Dp17O"].iloc[-1],
        "Eocene\nmeteoric water",
        va="center", ha="left", color=C["meteoric"], fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )
df_mix = mix_d17o(d18o_a=d18o_ice_free, dp17o_a=dp17o_ice_free,
                  d18o_b=d18o_meteoric_water, dp17o_b=dp17o_meteoric_water, step=10)
ax.scatter(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)
ax.annotate("50%",
            (prime(df_mix["mix_d18O"].iloc[5]), df_mix["mix_Dp17O"].iloc[5]),
            xytext=(1, 10), textcoords="offset points",
            va="center", ha="left", color="k", fontsize=4,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.5,
                    connectionstyle="arc3,rad=-0.2", shrinkA=0),
            )

# Add ice to the mix -----------------------------------------------------------
# Mixing line between modern seawater and Greenland ice
df_mix = mix_d17o(d18o_a=d18o_ice_free, dp17o_a=dp17o_ice_free,
                  d18o_b=d18o_ice, dp17o_b=dp17o_ice)
ax.plot(prime(df_mix["mix_d18O"]), df_mix["mix_Dp17O"], color=C["snowmelt"], linestyle="-", lw=1.5)
ax.scatter(prime(df_mix["mix_d18O"].iloc[-1]), df_mix["mix_Dp17O"].iloc[[-1]],
            marker = "*", ec=C["snowmelt"], fc="w", s=60, zorder=100)
ax.text(prime(d18o_ice) + 2, dp17o_ice,
        "snowmelt",
        va="center", ha="left", color=C["snowmelt"], fontsize=6,
        # bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8)
        )


df_mix_meteoric_ice = mix_d17o(d18o_a=d18o_meteoric_water, dp17o_a=dp17o_meteoric_water,
                               d18o_b=d18o_ice, dp17o_b=dp17o_ice, step=100)
ax.plot(
    prime(df_mix_meteoric_ice["mix_d18O"]),
    df_mix_meteoric_ice["mix_Dp17O"],
    color=C["mixed_water"],
    linestyle="--",
    lw=1.2,
    zorder=1,
)
df_mix_meteoric_ice = mix_d17o(d18o_a=d18o_meteoric_water, dp17o_a=dp17o_meteoric_water,
                               d18o_b=d18o_ice, dp17o_b=dp17o_ice, step=4)
ax.scatter(prime(df_mix_meteoric_ice["mix_d18O"]), df_mix_meteoric_ice["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)
ax.annotate("50%",
            (prime(df_mix_meteoric_ice["mix_d18O"].iloc[2]), df_mix_meteoric_ice["mix_Dp17O"].iloc[2]),
            xytext=(1, 10), textcoords="offset points",
            va="center", ha="left", color="k", fontsize=4,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.5,
                    connectionstyle="arc3,rad=-0.2", shrinkA=0),
            )


# for every point in df_mix, draw mixing line between the endmember and ice-free ocean
meteoric_rgb = np.array(mcolors.to_rgb(C["meteoric"]))
ice_rgb = np.array(mcolors.to_rgb(C["snowmelt"]))

for i in range(len(df_mix_meteoric_ice)):
    blend_ratio = i / (len(df_mix_meteoric_ice) - 1) if len(df_mix_meteoric_ice) > 1 else 0
    line_color = mcolors.to_hex((1 - blend_ratio) * meteoric_rgb + blend_ratio * ice_rgb)

    mix = mix_d17o(d18o_a=d18o_ice_free, dp17o_a=dp17o_ice_free,
                  d18o_b=df_mix_meteoric_ice["mix_d18O"].iloc[i], dp17o_b=df_mix_meteoric_ice["mix_Dp17O"].iloc[i],
                  step=50)
    ax.scatter(prime(mix["mix_d18O"]), mix["mix_Dp17O"], fc=line_color, ec ="none", zorder=0, s=1, marker=".")
    mix = mix_d17o(d18o_a=d18o_ice_free, dp17o_a=dp17o_ice_free,
                  d18o_b=df_mix_meteoric_ice["mix_d18O"].iloc[i], dp17o_b=df_mix_meteoric_ice["mix_Dp17O"].iloc[i],
                  step=10)
    ax.scatter(prime(mix["mix_d18O"]), mix["mix_Dp17O"], marker="|", c="k", s=1, zorder = 9)

# Add legend -------------------------------------------------------------------
legend = ax.legend(handles=legend_handles, loc='lower left',
                   fancybox=False, markerscale=0.9,
                   framealpha=1, fontsize=6, edgecolor="k", facecolor="none")
for text in legend.get_texts():
    if text.get_text() in legend_sections.values():
        text.set_weight("bold")


# Save figure
save_figure(fig, figures_dir, "glendonites")

# Manuscript parts:

print(f'\nModern glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "modern"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "modern"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "modern"])}')
print(f'pre-PETM glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "pre-PETM"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "pre-PETM"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "pre-PETM"])}')
print(f'PETM glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "PETM"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "PETM"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "PETM"])}')
print(f'early Eocene glendonites: ∆17O = {df_carbonates[df_carbonates["age"] == "early Eocene"]["Dp17O_AC"].mean():.0f} ppm and d18O = {df_carbonates[df_carbonates["age"] == "early Eocene"]["d18O_AC"].mean():.1f} ‰. Count = {len(df_carbonates[df_carbonates["age"] == "early Eocene"])}')
