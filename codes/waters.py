# Import libraries ------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from _helper_functions import *
from _plot_style import *
from pathlib import Path


# Define directories ----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
source_data_dir = Path(__file__).resolve().parent.parent / "source_data"
export_tables_dir = Path(__file__).resolve().parent.parent / "export_tables" / "csv"


def _resolve_data_dirs(source_data_dir=None, export_tables_dir=None):
    if source_data_dir is None:
        source_data_dir = Path(__file__).resolve().parent.parent / "source_data"
    if export_tables_dir is None:
        export_tables_dir = Path(__file__).resolve().parent.parent / "export_tables" / "csv"
    return source_data_dir, export_tables_dir


def load_water_data(source_data_dir=None, export_tables_dir=None):
    source_data_dir, export_tables_dir = _resolve_data_dirs(source_data_dir, export_tables_dir)

    df_sw = pd.read_csv(source_data_dir / "water_data" / "modern_seawater.csv")
    df_Aron = pd.read_csv(source_data_dir / "water_data" / "Aron_2021.csv")
    df_Pang = pd.read_csv(source_data_dir / "water_data" / "Pang_2015.csv")
    df_Touzeau = pd.read_csv(source_data_dir / "water_data" / "Touzeau_2016.csv")
    df_Landais2008S1 = pd.read_csv(source_data_dir / "water_data" / "Landais_2008_S1.csv")
    df_Landais2008S2 = pd.read_csv(source_data_dir / "water_data" / "Landais_2008_S2.csv")
    df_Steen_Larsen_2014 = pd.read_csv(source_data_dir / "water_data" / "Steen_Larsen_2014.csv")
    df_Landais_2012b = pd.read_csv(source_data_dir / "water_data" / "Landais_2012.csv")
    df_Beershaw = pd.read_csv(source_data_dir / "water_data" / "Beershaw_2020.csv")

    df_carbonates = pd.read_csv(export_tables_dir / "combined_replicate_data.csv")

    return {
        "df_sw": df_sw,
        "df_Aron": df_Aron,
        "df_Pang": df_Pang,
        "df_Touzeau": df_Touzeau,
        "df_Landais2008S1": df_Landais2008S1,
        "df_Landais2008S2": df_Landais2008S2,
        "df_Steen_Larsen_2014": df_Steen_Larsen_2014,
        "df_Landais_2012b": df_Landais_2012b,
        "df_Beershaw": df_Beershaw,
        "df_carbonates": df_carbonates,
    }


# Calcualte ice-free ocean composition from mass balance ----------------------

def calculate_ice_free_ocean(
    d18O_ocean_now,
    Dp17O_ocean_now,
    d18O_antarctic,
    Dp17O_antarctic,
    d18O_Greenland,
    Dp17O_Greenland,
):
    """Return ice-free ocean d18O and Dp17O from three-reservoir mass balance."""
    # Water mass volumes (×10^6 km^3)
    V_ocean = 1360.0
    V_antarctic = 24.7  # Lythe et al. (2001)
    V_Greenland = 2.9  # Bamber et al. (2001)
    V_total = V_ocean + V_antarctic + V_Greenland

    d18O_num = (
        d18O_ocean_now * V_ocean
        + d18O_antarctic * V_antarctic
        + d18O_Greenland * V_Greenland
    )
    d18O_ice_free = d18O_num / V_total

    d17O_ocean_now = d17O(d18O_ocean_now, Dp17O_ocean_now)
    d17O_antarctic = d17O(d18O_antarctic, Dp17O_antarctic)
    d17O_Greenland = d17O(d18O_Greenland, Dp17O_Greenland)

    d17O_num = (
        d17O_ocean_now * V_ocean
        + d17O_antarctic * V_antarctic
        + d17O_Greenland * V_Greenland
    )
    d17O_ice_free = d17O_num / V_total
    Dp17O_ice_free = Dp17O(d17O_ice_free, d18O_ice_free)

    return d18O_ice_free, Dp17O_ice_free


def _build_water_state(source_data_dir=None, export_tables_dir=None, latitude=60, altitude=1500, Dp17O_meteoric_water=30):
    data = load_water_data(source_data_dir, export_tables_dir)

    df_sw = data["df_sw"]
    df_Aron = data["df_Aron"]
    df_Pang = data["df_Pang"]
    df_Touzeau = data["df_Touzeau"]
    df_Landais2008S1 = data["df_Landais2008S1"]
    df_Landais2008S2 = data["df_Landais2008S2"]
    df_Steen_Larsen_2014 = data["df_Steen_Larsen_2014"]
    df_Landais_2012b = data["df_Landais_2012b"]
    df_Beershaw = data["df_Beershaw"]

    d18O_modern_seawater = df_sw["d18O"].median()
    d18O_modern_seawater_SD = df_sw["d18O"].std()
    Dp17O_modern_seawater = df_sw["Dp17O"].median()
    Dp17O_modern_seawater_SD = df_sw["Dp17O"].std()

    df_Landais2008S1 = df_Landais2008S1.loc[:, ["Elevation", "δ18O", "Dp17O"]]
    df_Landais2008S1_mean = df_Landais2008S1.groupby("Elevation").mean().reset_index()
    df_Landais2008S2 = df_Landais2008S2.loc[:, ["Depth", "δ18O", "Dp17O"]]
    df_Landais2008S2_mean = df_Landais2008S2.groupby("Depth").mean().reset_index()

    df_antarctic_combined = pd.concat([df_Pang, df_Touzeau, df_Landais2008S1_mean, df_Landais2008S2_mean])
    df_antarctic_combined = df_antarctic_combined.loc[:, ["δ18O", "Dp17O"]]
    d18O_antarctic = df_antarctic_combined["δ18O"].median()
    Dp17O_antarctic = df_antarctic_combined["Dp17O"].median()

    df_greenland_combined = pd.concat([df_Steen_Larsen_2014, df_Landais_2012b])
    df_greenland_combined = df_greenland_combined.loc[:, ["δ18O", "Dp17O"]]
    d18O_greenland = df_greenland_combined["δ18O"].median()
    Dp17O_greenland = df_greenland_combined["Dp17O"].median()

    df_meteoric_combined = pd.concat([df_Aron, df_Beershaw])
    df_meteoric_combined = df_meteoric_combined.loc[:, ["δ18O", "Dp17O"]]

    d18O_ice_free, Dp17O_ice_free = calculate_ice_free_ocean(
        d18O_modern_seawater,
        Dp17O_modern_seawater,
        d18O_antarctic,
        Dp17O_antarctic,
        d18O_greenland,
        Dp17O_greenland,
    )

    Dd18O_ice_free_modern = d18O_modern_seawater - d18O_ice_free
    Dd17O_ice_free_modern = Dp17O_modern_seawater - Dp17O_ice_free

    d18O_greenland_shifted = d18O_greenland - Dd18O_ice_free_modern
    Dp17O_greenland_shifted = Dp17O_greenland - Dd17O_ice_free_modern

    d18O_meteoric_water = -0.0051 * latitude**2 + 0.1805 * latitude - 0.002 * altitude - 5.247
    d18O_meteoric_water_shifted = d18O_meteoric_water - Dd18O_ice_free_modern
    Dp17O_meteoric_water_shifted = Dp17O_meteoric_water - Dd17O_ice_free_modern

    return {
        "df_sw": df_sw,
        "df_antarctic_combined": df_antarctic_combined,
        "df_greenland_combined": df_greenland_combined,
        "df_meteoric_combined": df_meteoric_combined,
        "d18O_modern_seawater": d18O_modern_seawater,
        "d18O_modern_seawater_SD": d18O_modern_seawater_SD,
        "Dp17O_modern_seawater": Dp17O_modern_seawater,
        "Dp17O_modern_seawater_SD": Dp17O_modern_seawater_SD,
        "d18O_antarctic": d18O_antarctic,
        "Dp17O_antarctic": Dp17O_antarctic,
        "d18O_greenland": d18O_greenland,
        "Dp17O_greenland": Dp17O_greenland,
        "d18O_ice_free": d18O_ice_free,
        "Dp17O_ice_free": Dp17O_ice_free,
        "Dd18O_ice_free_modern": Dd18O_ice_free_modern,
        "Dd17O_ice_free_modern": Dd17O_ice_free_modern,
        "d18O_greenland_shifted": d18O_greenland_shifted,
        "Dp17O_greenland_shifted": Dp17O_greenland_shifted,
        "d18O_meteoric_water": d18O_meteoric_water,
        "Dp17O_meteoric_water": Dp17O_meteoric_water,
        "d18O_meteoric_water_shifted": d18O_meteoric_water_shifted,
        "Dp17O_meteoric_water_shifted": Dp17O_meteoric_water_shifted,
        "latitude": latitude,
        "altitude": altitude,
    }


def calculate_water_endmembers(source_data_dir=None, export_tables_dir=None, latitude=60, altitude=1500, Dp17O_meteoric_water=30):
    state = _build_water_state(source_data_dir, export_tables_dir, latitude, altitude, Dp17O_meteoric_water)
    return (
        state["d18O_ice_free"],
        state["Dp17O_ice_free"],
        state["d18O_greenland_shifted"],
        state["Dp17O_greenland_shifted"],
        state["d18O_modern_seawater"],
        state["Dp17O_modern_seawater"],
        state["d18O_meteoric_water_shifted"],
        state["Dp17O_meteoric_water_shifted"],
    )


def print_water_summary(state):
    print(
        f"Modern seawater (median): d18O = {state['d18O_modern_seawater']:.2f}(±{state['d18O_modern_seawater_SD']:.2f}; 1σ)‰, "
        f"Dp17O = {state['Dp17O_modern_seawater']:.2f}(±{state['Dp17O_modern_seawater_SD']:.2f}; 1σ)ppm"
    )
    print(f"Ice-free ocean: Δ'17O = {state['Dp17O_ice_free']:.2f} ppm, δ18O = {state['d18O_ice_free']:.2f} ‰")
    print(
        f"Difference between modern seawater and ice-free ocean: Δd18O = {state['Dd18O_ice_free_modern']:.2f}‰, "
        f"ΔDp17O = {state['Dd17O_ice_free_modern']:.2f} ppm"
    )
    print(f"Antarctic ice average: d18O = {state['d18O_antarctic']:.2f}‰, Dp17O = {state['Dp17O_antarctic']:.2f} ppm")
    print(f"Greenland snow average: d18O = {state['d18O_greenland']:.2f}‰, Dp17O = {state['Dp17O_greenland']:.2f} ppm")
    print(
        f"Greenland snow average shifted to ice-free ocean: d18O = {state['d18O_greenland_shifted']:.2f}‰, "
        f"Dp17O = {state['Dp17O_greenland_shifted']:.2f} ppm"
    )
    print(
        f"Estimated meteoric water: d18O = {state['d18O_meteoric_water_shifted']:.2f}‰, "
        f"Dp17O = {state['Dp17O_meteoric_water_shifted']:.2f} ppm"
    )


def plot_water_compositions(state, figures_dir):
    fig, ax = plt.subplots()

    ax.scatter(prime(state["df_sw"]["d18O"]), state["df_sw"]["Dp17O"],
               marker='o', fc="#858379", alpha=0.5, ec="none", label="modern seawater")
    ax.scatter(prime(state["d18O_modern_seawater"]), state["Dp17O_modern_seawater"],
               marker="*", fc="#858379", ec="k", s=80, zorder=10)

    ax.scatter(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"],
               marker='*', ec="k", fc="#1455C0", s=60, zorder=10)
    ax.text(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"] - 5,
            "ice-free\nseawater",
            va="top", ha="center", color="#1455C0", fontsize=6,
            bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8))

    ax.annotate("", xy=(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"]),
                xytext=(prime(state["d18O_modern_seawater"]), state["Dp17O_modern_seawater"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter(prime(state["df_antarctic_combined"]["δ18O"]), state["df_antarctic_combined"]["Dp17O"],
               marker='o', fc="#004B6D", ec="none", alpha=0.5, zorder=2, label="Antarctic ice and snow")
    ax.scatter(prime(state["d18O_antarctic"]), state["Dp17O_antarctic"],
               marker="*", fc="#004B6D", ec="k", s=80, zorder=10)

    ax.scatter(prime(state["df_greenland_combined"]["δ18O"]), state["df_greenland_combined"]["Dp17O"],
               marker='o', fc="#309FD1", ec="none", alpha=0.5, zorder=2, label="Greenland snow")

    ax.scatter(prime(state["d18O_greenland"]), state["Dp17O_greenland"],
               marker="*", fc="#309FD1", ec="k", zorder=10)

    ax.scatter(prime(state["d18O_greenland_shifted"]), state["Dp17O_greenland_shifted"],
               marker="*", fc="#309FD1", ec="k", s=80, zorder=10)
    ax.annotate("", xy=(prime(state["d18O_greenland_shifted"]), state["Dp17O_greenland_shifted"]),
                xytext=(prime(state["d18O_greenland"]), state["Dp17O_greenland"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter(prime(state["df_meteoric_combined"]["δ18O"]), state["df_meteoric_combined"]["Dp17O"],
               marker='o', fc="#3CB5AE", ec="none", zorder=2, alpha=0.5,
               label="meteoric waters")

    ax.scatter(prime(state["d18O_meteoric_water"]), state["Dp17O_meteoric_water"],
               marker="*", fc="#3CB5AE", ec="k", zorder=10)
    ax.text(prime(state["d18O_meteoric_water"]) + 2, state["Dp17O_meteoric_water"],
            f"modeled modern\nmeteoric water\n({state['latitude']}° N, {state['altitude']} m)",
            va="center", ha="left", color="#3CB5AE", fontsize=6,
            bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8))
    ax.scatter(prime(state["d18O_meteoric_water_shifted"]), state["Dp17O_meteoric_water_shifted"],
               marker="*", fc="#3CB5AE", ec="k", s=80, zorder=10)
    ax.annotate("", xy=(prime(state["d18O_meteoric_water_shifted"]), state["Dp17O_meteoric_water_shifted"]),
                xytext=(prime(state["d18O_meteoric_water"]), state["Dp17O_meteoric_water"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter([], [], marker="*", fc="none", ec="k", s=80, zorder=10, label="median compositions")

    ax.legend(loc='lower left', fancybox=False, markerscale=0.9,
              framealpha=1, fontsize=6, edgecolor="k")
    ax.set_ylabel(r"$\Delta\prime^{17}$O (ppm)")
    ax.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW)")

    save_figure(fig, figures_dir, "waters")


def main():
    apply_plot_style()
    figures_dir = Path(__file__).resolve().parent.parent / "figures"
    state = _build_water_state()
    print_water_summary(state)
    plot_water_compositions(state, figures_dir)


if __name__ == "__main__":
    main()