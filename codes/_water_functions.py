# Import libraries -------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from _helper_functions import *
from _plot_style import *
from pathlib import Path
from typing import Any

# Define directories -----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
literature_data_dir = Path(__file__).resolve().parent.parent / "data" / "literature"
export_tables_dir = Path(__file__).resolve().parent.parent / "data" / "exported" / "csv"


def _resolve_data_dirs(
    literature_data_dir: Path | None = None,
    export_tables_dir: Path | None = None,
) -> tuple[Path, Path]:
    if literature_data_dir is None:
        literature_data_dir = Path(__file__).resolve().parent.parent / "data" / "literature"
    if export_tables_dir is None:
        export_tables_dir = Path(__file__).resolve().parent.parent / "data" / "exported" / "csv"
    return literature_data_dir, export_tables_dir


def load_water_data(
    literature_data_dir: Path | None = None,
    export_tables_dir: Path | None = None,
) -> dict[str, pd.DataFrame]:
    literature_data_dir, export_tables_dir = _resolve_data_dirs(literature_data_dir, export_tables_dir)

    df_sw = pd.read_csv(literature_data_dir / "modern_seawater.csv")
    df_aron = pd.read_csv(literature_data_dir / "Aron_2021.csv")
    df_pang = pd.read_csv(literature_data_dir / "Pang_2015.csv")
    df_touzeau = pd.read_csv(literature_data_dir / "Touzeau_2016.csv")
    df_landais2008s1 = pd.read_csv(literature_data_dir / "Landais_2008_S1.csv")
    df_landais2008s2 = pd.read_csv(literature_data_dir / "Landais_2008_S2.csv")
    df_steen_larsen_2014 = pd.read_csv(literature_data_dir / "Steen_Larsen_2014.csv")
    df_landais_2012b = pd.read_csv(literature_data_dir / "Landais_2012.csv")
    df_beershaw = pd.read_csv(literature_data_dir / "Beershaw_2020.csv")

    df_carbonates = pd.read_csv(export_tables_dir / "combined_replicate_data.csv")

    return {
        "df_sw": df_sw,
        "df_Aron": df_aron,
        "df_Pang": df_pang,
        "df_Touzeau": df_touzeau,
        "df_Landais2008S1": df_landais2008s1,
        "df_Landais2008S2": df_landais2008s2,
        "df_Steen_Larsen_2014": df_steen_larsen_2014,
        "df_Landais_2012b": df_landais_2012b,
        "df_Beershaw": df_beershaw,
        "df_carbonates": df_carbonates,
    }

# Calculate ice-free ocean composition from mass balance -----------------------
def calculate_ice_free_ocean(
    d18o_ocean_now: float,
    dp17o_ocean_now: float,
    d18o_antarctic: float,
    dp17o_antarctic: float,
    d18o_greenland: float,
    dp17o_greenland: float,
) -> tuple[float, float]:
    """Return ice-free ocean d18O and Dp17O from three-reservoir mass balance."""
    # Water mass volumes (×10^6 km^3)
    v_ocean = 1360.0
    v_antarctic = 24.7  # Lythe et al. (2001)
    v_greenland = 2.9  # Bamber et al. (2001)
    v_total = v_ocean + v_antarctic + v_greenland

    d18o_num = (
        d18o_ocean_now * v_ocean
        + d18o_antarctic * v_antarctic
        + d18o_greenland * v_greenland
    )
    d18o_ice_free = d18o_num / v_total

    d17o_ocean_now = d17o(d18o_ocean_now, dp17o_ocean_now)
    d17o_antarctic = d17o(d18o_antarctic, dp17o_antarctic)
    d17o_greenland = d17o(d18o_greenland, dp17o_greenland)

    d17o_num = (
        d17o_ocean_now * v_ocean
        + d17o_antarctic * v_antarctic
        + d17o_greenland * v_greenland
    )
    d17o_ice_free = d17o_num / v_total
    dp17o_ice_free = dp17o(d17o_ice_free, d18o_ice_free)

    return d18o_ice_free, dp17o_ice_free


def _build_water_state(
    literature_data_dir: Path | None = None,
    export_tables_dir: Path | None = None,
    latitude: float = 60,
    altitude: float = 1500,
    dp17o_meteoric_water: float = 30,
) -> dict[str, Any]:
    data = load_water_data(literature_data_dir, export_tables_dir)

    df_sw = data["df_sw"]
    df_aron = data["df_Aron"]
    df_pang = data["df_Pang"]
    df_touzeau = data["df_Touzeau"]
    df_landais2008s1 = data["df_Landais2008S1"]
    df_landais2008s2 = data["df_Landais2008S2"]
    df_steen_larsen_2014 = data["df_Steen_Larsen_2014"]
    df_landais_2012b = data["df_Landais_2012b"]
    df_beershaw = data["df_Beershaw"]

    d18o_modern_seawater = df_sw["d18O"].median()
    d18o_modern_seawater_sd = df_sw["d18O"].std()
    dp17o_modern_seawater = df_sw["Dp17O"].median()
    dp17o_modern_seawater_sd = df_sw["Dp17O"].std()

    df_landais2008s1 = df_landais2008s1.loc[:, ["Elevation", "δ18O", "Dp17O"]]
    df_landais2008s1_mean = df_landais2008s1.groupby("Elevation").mean().reset_index()
    df_landais2008s2 = df_landais2008s2.loc[:, ["Depth", "δ18O", "Dp17O"]]
    df_landais2008s2_mean = df_landais2008s2.groupby("Depth").mean().reset_index()

    df_antarctic_combined = pd.concat([df_pang, df_touzeau, df_landais2008s1_mean, df_landais2008s2_mean])
    df_antarctic_combined = df_antarctic_combined.loc[:, ["δ18O", "Dp17O"]]
    d18o_antarctic = df_antarctic_combined["δ18O"].median()
    dp17o_antarctic = df_antarctic_combined["Dp17O"].median()

    df_greenland_combined = pd.concat([df_steen_larsen_2014, df_landais_2012b])
    df_greenland_combined = df_greenland_combined.loc[:, ["δ18O", "Dp17O"]]
    d18o_greenland = df_greenland_combined["δ18O"].median()
    dp17o_greenland = df_greenland_combined["Dp17O"].median()

    df_meteoric_combined = pd.concat([df_aron, df_beershaw])
    df_meteoric_combined = df_meteoric_combined.loc[:, ["δ18O", "Dp17O"]]

    d18o_ice_free, dp17o_ice_free = calculate_ice_free_ocean(
        d18o_modern_seawater,
        dp17o_modern_seawater,
        d18o_antarctic,
        dp17o_antarctic,
        d18o_greenland,
        dp17o_greenland,
    )

    dd18o_ice_free_modern = d18o_modern_seawater - d18o_ice_free
    dd17o_ice_free_modern = dp17o_modern_seawater - dp17o_ice_free

    d18o_greenland_shifted = d18o_greenland - dd18o_ice_free_modern
    dp17o_greenland_shifted = dp17o_greenland - dd17o_ice_free_modern

    d18o_meteoric_water = -0.0051 * latitude**2 + 0.1805 * latitude - 0.002 * altitude - 5.247
    d18o_meteoric_water_shifted = d18o_meteoric_water - dd18o_ice_free_modern
    dp17o_meteoric_water_shifted = dp17o_meteoric_water - dd17o_ice_free_modern

    return {
        "df_sw": df_sw,
        "df_antarctic_combined": df_antarctic_combined,
        "df_greenland_combined": df_greenland_combined,
        "df_meteoric_combined": df_meteoric_combined,
        "d18O_modern_seawater": d18o_modern_seawater,
        "d18O_modern_seawater_SD": d18o_modern_seawater_sd,
        "Dp17O_modern_seawater": dp17o_modern_seawater,
        "Dp17O_modern_seawater_SD": dp17o_modern_seawater_sd,
        "d18O_antarctic": d18o_antarctic,
        "Dp17O_antarctic": dp17o_antarctic,
        "d18O_greenland": d18o_greenland,
        "Dp17O_greenland": dp17o_greenland,
        "d18O_ice_free": d18o_ice_free,
        "Dp17O_ice_free": dp17o_ice_free,
        "Dd18O_ice_free_modern": dd18o_ice_free_modern,
        "Dd17O_ice_free_modern": dd17o_ice_free_modern,
        "d18O_greenland_shifted": d18o_greenland_shifted,
        "Dp17O_greenland_shifted": dp17o_greenland_shifted,
        "d18O_meteoric_water": d18o_meteoric_water,
        "Dp17O_meteoric_water": dp17o_meteoric_water,
        "d18O_meteoric_water_shifted": d18o_meteoric_water_shifted,
        "Dp17O_meteoric_water_shifted": dp17o_meteoric_water_shifted,
        "latitude": latitude,
        "altitude": altitude,
    }


def calculate_water_endmembers(
    literature_data_dir: Path | None = None,
    export_tables_dir: Path | None = None,
    latitude: float = 60,
    altitude: float = 1500,
    dp17o_meteoric_water: float = 30,
) -> tuple[float, ...]:
    state = _build_water_state(literature_data_dir, export_tables_dir, latitude, altitude, dp17o_meteoric_water)
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


def print_water_summary(state: dict[str, Any]) -> None:
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


def plot_water_compositions(state: dict[str, Any], figures_dir: Path) -> None:
    # Create figure ------------------------------------------------------------
    fig, ax = plt.subplots()

    # Setup axes ---------------------------------------------------------------
    ax.set_ylabel(r"$\Delta\prime^{17}$O (ppm)")
    ax.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW)")

    # Plot data ----------------------------------------------------------------

    ax.scatter(prime(state["df_sw"]["d18O"]), state["df_sw"]["Dp17O"],
               marker='o', fc=C["wgray5"], alpha=0.5, ec="none", label="modern seawater")
    ax.scatter(prime(state["d18O_modern_seawater"]), state["Dp17O_modern_seawater"],
               marker="*", fc=C["wgray5"], ec="k", s=80, zorder=10)

    ax.scatter(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"],
               marker='*', ec="k", fc=C["blue"], s=60, zorder=10)
    ax.text(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"] - 5,
            "ice-free\nseawater",
            va="top", ha="center", color=C["blue"], fontsize=6,
            bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8))

    ax.annotate("", xy=(prime(state["d18O_ice_free"]), state["Dp17O_ice_free"]),
                xytext=(prime(state["d18O_modern_seawater"]), state["Dp17O_modern_seawater"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter(prime(state["df_antarctic_combined"]["δ18O"]), state["df_antarctic_combined"]["Dp17O"],
               marker='o', fc=C["cyan8"], ec="none", alpha=0.5, zorder=2, label="Antarctic ice and snow")
    ax.scatter(prime(state["d18O_antarctic"]), state["Dp17O_antarctic"],
               marker="*", fc=C["cyan8"], ec="k", s=80, zorder=10)

    ax.scatter(prime(state["df_greenland_combined"]["δ18O"]), state["df_greenland_combined"]["Dp17O"],
               marker='o', fc=C["cyan5"], ec="none", alpha=0.5, zorder=2, label="Greenland snow")

    ax.scatter(prime(state["d18O_greenland"]), state["Dp17O_greenland"],
               marker="*", fc=C["cyan5"], ec="k", zorder=10)

    ax.scatter(prime(state["d18O_greenland_shifted"]), state["Dp17O_greenland_shifted"],
               marker="*", fc=C["cyan5"], ec="k", s=80, zorder=10)
    ax.annotate("", xy=(prime(state["d18O_greenland_shifted"]), state["Dp17O_greenland_shifted"]),
                xytext=(prime(state["d18O_greenland"]), state["Dp17O_greenland"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter(prime(state["df_meteoric_combined"]["δ18O"]), state["df_meteoric_combined"]["Dp17O"],
               marker='o', fc=C["turquoise4"], ec="none", zorder=2, alpha=0.5,
               label="meteoric waters")

    ax.scatter(prime(state["d18O_meteoric_water"]), state["Dp17O_meteoric_water"],
               marker="*", fc=C["turquoise4"], ec="k", zorder=10)
    ax.text(prime(state["d18O_meteoric_water"]) + 2, state["Dp17O_meteoric_water"],
            f"modeled modern\nmeteoric water\n({state['latitude']}° N, {state['altitude']} m)",
            va="center", ha="left", color=C["turquoise4"], fontsize=6,
            bbox=dict(fc="white", ec="none", pad=0.2, alpha=0.8))
    ax.scatter(prime(state["d18O_meteoric_water_shifted"]), state["Dp17O_meteoric_water_shifted"],
               marker="*", fc=C["turquoise4"], ec="k", s=80, zorder=10)
    ax.annotate("", xy=(prime(state["d18O_meteoric_water_shifted"]), state["Dp17O_meteoric_water_shifted"]),
                xytext=(prime(state["d18O_meteoric_water"]), state["Dp17O_meteoric_water"]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=0.5, connectionstyle="arc3,rad=3"),
                )

    ax.scatter([], [], marker="*", fc="none", ec="k", s=80, zorder=10, label="median compositions")

    # Add legend ---------------------------------------------------------------
    ax.legend(loc='lower left', fancybox=False, markerscale=0.9,
              framealpha=1, fontsize=6, edgecolor="k")

    save_figure(fig, figures_dir, "waters")
