# Import libraries ------------------------------------------------------------
from pathlib import Path
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from _helper_functions import *
from _plot_style import *


# Define directories ----------------------------------------------------------
sources_dir = Path(__file__).resolve().parent.parent
source_data_dir = sources_dir / "source_data"
export_tables_dir = sources_dir / "export_tables" / "csv"
figures_dir = sources_dir / "figures"


# Apply plot style ------------------------------------------------------------
apply_plot_style()
plt.rcParams["figure.figsize"] = (6.5, 3.5)


# Additional functions --------------------------------------------------------

def retrieve_seawater_parameters(samples, model_path):
    """Retrieve monthly mean seawater d18O and temperature by site and depth."""

    d18O_SW_MODEL_SIGMA = 0.2
    TEMPERATURE_MODEL_SIGMA = 0.5

    with nc.Dataset(model_path) as dataset:
        latitudes = np.asarray(dataset["lat_1deg_center"][:])
        longitudes = np.asarray(dataset["lon_1deg_center"][:])

        # Add a surface level by duplicating the shallowest model layer, as in
        # SK_2_get_sw_parameters.py. Model depths are stored as negative values.
        depths = -np.concatenate(
            (
                np.asarray(dataset["depth_edge"][:1]),
                np.asarray(dataset["depth_center"][:]),
            )
        )
        d18O = np.ma.concatenate(
            (
                dataset["D18O_1deg"][:, :1, :, :],
                dataset["D18O_1deg"][:, :, :, :],
            ),
            axis=1,
        ).filled(np.nan)
        temperature = np.ma.concatenate(
            (
                dataset["THETA_1deg"][:, :1, :, :],
                dataset["THETA_1deg"][:, :, :, :],
            ),
            axis=1,
        ).filled(np.nan)

    grid_x = np.cos(np.deg2rad(longitudes)) * np.cos(np.deg2rad(latitudes))
    grid_y = np.sin(np.deg2rad(longitudes)) * np.cos(np.deg2rad(latitudes))
    grid_z = np.sin(np.deg2rad(latitudes))

    retrieved = []
    for sample in samples.itertuples(index=False):
        if pd.isna(sample.lat) or pd.isna(sample.long) or pd.isna(sample.water_depth):
            raise ValueError(
                f"Missing latitude, longitude, or water depth for {sample.SampleName}."
            )

        sample_x = np.cos(np.deg2rad(sample.long)) * np.cos(np.deg2rad(sample.lat))
        sample_y = np.sin(np.deg2rad(sample.long)) * np.cos(np.deg2rad(sample.lat))
        sample_z = np.sin(np.deg2rad(sample.lat))
        squared_distance = (
            (grid_x - sample_x) ** 2
            + (grid_y - sample_y) ** 2
            + (grid_z - sample_z) ** 2
        )

        # Select the nearest horizontal grid cell containing water at the first
        # model level at or below the sample depth. For samples deeper than the
        # model, use the deepest model level rather than extrapolating.
        depth_index = min(
            np.searchsorted(depths, sample.water_depth, side="left"),
            len(depths) - 1,
        )
        squared_distance[~np.isfinite(d18O[0, depth_index, :, :])] = np.inf
        if not np.isfinite(squared_distance).any():
            raise ValueError(
                f"No valid seawater d18O grid cell found for {sample.SampleName}."
            )
        lat_index, lon_index = np.unravel_index(
            np.argmin(squared_distance), squared_distance.shape
        )

        monthly_d18O_values = []
        monthly_temperature_values = []
        interpolation_depths = []
        for monthly_d18O_profile, monthly_temperature_profile in zip(
            d18O[:, :, lat_index, lon_index],
            temperature[:, :, lat_index, lon_index],
        ):
            valid = np.isfinite(monthly_d18O_profile) & np.isfinite(
                monthly_temperature_profile
            )
            if not valid.any():
                continue

            valid_depths = depths[valid]
            interpolation_depth = np.clip(
                sample.water_depth, valid_depths[0], valid_depths[-1]
            )
            monthly_d18O_values.append(
                np.interp(
                    interpolation_depth,
                    valid_depths,
                    monthly_d18O_profile[valid],
                )
            )
            monthly_temperature_values.append(
                np.interp(
                    interpolation_depth,
                    valid_depths,
                    monthly_temperature_profile[valid],
                )
            )
            interpolation_depths.append(interpolation_depth)

        if not monthly_d18O_values:
            raise ValueError(
                "No valid monthly seawater d18O and temperature values found for "
                f"{sample.SampleName}."
            )

        monthly_d18O_values = np.asarray(monthly_d18O_values)
        monthly_temperature_values = np.asarray(monthly_temperature_values)
        monthly_d18O_sigma = (
            monthly_d18O_values.std(ddof=1)
            if len(monthly_d18O_values) > 1
            else 0.0
        )
        monthly_temperature_sigma = (
            monthly_temperature_values.std(ddof=1)
            if len(monthly_temperature_values) > 1
            else 0.0
        )
        model_depth = min(interpolation_depths)
        if not np.isclose(model_depth, sample.water_depth):
            print(
                f"{sample.SampleName} has a water depth of {sample.water_depth:g} m, but the deepest available model value is at {model_depth:g} m.\n"
                "Using the deepest available seawater d18O and temperature values.\n"
            )

        retrieved.append(
            {
                "SampleName": sample.SampleName,
                "water_d18O": round(monthly_d18O_values.mean(), 3),
                "water_d18O_error": round(
                    np.hypot(d18O_SW_MODEL_SIGMA, monthly_d18O_sigma), 3
                ),
                "water_temperature": monthly_temperature_values.mean(),
                "water_temperature_error": np.hypot(
                    TEMPERATURE_MODEL_SIGMA,
                    monthly_temperature_sigma,
                ),
            }
        )

    return pd.DataFrame(retrieved)


def calculate_fractionation_factor(
    mineral_value, water_value, mineral_error, water_error
):
    """Calculate mineral-water alpha and propagate independent errors."""
    mineral_ratio = 1000.0 + mineral_value
    water_ratio = 1000.0 + water_value
    alpha = mineral_ratio / water_ratio
    alpha_error = np.sqrt(
        (mineral_error / water_ratio) ** 2
        + (mineral_ratio * water_error / water_ratio**2) ** 2
    )
    return alpha, alpha_error


def calculate_d17O_error(d18O_value, Dp17O_value, d18O_error, Dp17O_error):
    """Propagate independent d18O and Dp17O errors through d17O()."""
    d18O_contribution = (
        d17O(d18O_value + d18O_error, Dp17O_value)
        - d17O(d18O_value - d18O_error, Dp17O_value)
    ) / 2.0
    Dp17O_contribution = (
        d17O(d18O_value, Dp17O_value + Dp17O_error)
        - d17O(d18O_value, Dp17O_value - Dp17O_error)
    ) / 2.0
    return np.hypot(d18O_contribution, Dp17O_contribution)


def calculate_theta_from_inputs(
    carbonate_d18O,
    carbonate_Dp17O,
    water_d18O,
    water_Dp17O,
    carbonate_d17O_offset,
):
    """Calculate theta while retaining the shared d18O dependence."""
    carbonate_d17O = (
        d17O(carbonate_d18O, carbonate_Dp17O) + carbonate_d17O_offset
    )
    water_d17O = d17O(water_d18O, water_Dp17O)
    alpha18 = (1000.0 + carbonate_d18O) / (1000.0 + water_d18O)
    alpha17 = (1000.0 + carbonate_d17O) / (1000.0 + water_d17O)
    return np.log(alpha17) / np.log(alpha18)


def a18_KO97(T):
    # Kim et O'Neil (1997) - calcite
    return np.exp((18.04 * 1000 / T - 32.18) / 1000)

def a18_Daeron2019(T):
    # Daeron et al. (2019) - calcite
    return np.exp((17.57 * 1000 / T - 29.13) / 1000)


def theta_GZ19_cc(T):
    # Guo and Zhou (2019) - calcite
    return 59.1047/T**2 + -1.4089/T + 0.5297


# def theta_W20(T):
#     # Wostbrock et al. (2020) - calcite
#     return -1.39 / T + 0.5305


# def theta_W20(T):
#     # Wostbrock et al. (2020) – aragonite
#     return -1.53 / T + 0.5305


# def theta_GZ19_ar(T):
#     # Guo and Zhou (2019) – aragonite
#     return 78.1173/T**2 + -1.5152/T + 0.5299


def theta_H18(T):

    # Hayles et al. (2018) - calcite
    K_calcite = 1.019124E+09 / T**5 + -2.117501E+07 / T**4 + 1.686453E+05 / T**3 + -5.784679E+02 / T**2 + 1.489666E-01 / T + 0.5304852
    B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    K_water = 7.625734E+06 / T**5 + 1.216102E+06 / T**4 + -2.135774E+04 / T**3 + 1.323782E+02 / T**2 + -4.931630E-01 / T + 0.5306551
    B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    a18 = np.exp(B_calcite) / np.exp(B_water)
    return K_calcite + (K_calcite-K_water) * (B_water / np.log(a18))


# Get modern data -------------------------------------------------------------
data = pd.read_csv(export_tables_dir / "combined_replicate_data.csv")
data = data[data["age"] == "modern"].copy()


# Retrieve seawater isotope values -------------------------------------------
model_values = retrieve_seawater_parameters(
    data,
    source_data_dir / "Breitkreuz_et_al_2018.nc",
)
data = data.merge(model_values, on="SampleName", validate="one_to_one")

# Assign modern seawater Dp17O and its 1-sigma uncertainty (ppm).
data["water_Dp17O"] = -5.0
data["water_Dp17O_error"] = 1.0


# Calculate clumped-based water d18O for comparison --------------------------
temperature_K = data["clumped_T"] + 273.15
temperature_K_error = data["clumped_T_err"]
daeron_alpha18 = a18_Daeron2019(temperature_K)
data["clumped_water_d18O"] = (
    (1000.0 + data["d18O_AC"]) / daeron_alpha18 - 1000.0
)

clumped_water_d18O_temperature_plus = (
    (1000.0 + data["d18O_AC"])
    / a18_Daeron2019(temperature_K + temperature_K_error)
    - 1000.0
)
clumped_water_d18O_temperature_minus = (
    (1000.0 + data["d18O_AC"])
    / a18_Daeron2019(temperature_K - temperature_K_error)
    - 1000.0
)
data["clumped_water_d18O_error"] = np.hypot(
    data["d18O_error"] / daeron_alpha18,
    (
        clumped_water_d18O_temperature_plus
        - clumped_water_d18O_temperature_minus
    )
    / 2.0,
)

# Calculate ikaite-water fractionation factors using database seawater d18O --
data["water_d17O"] = d17O(data["water_d18O"], data["water_Dp17O"])

# Propagate water d18O and Dp17O errors through the d17O calculation.
data["water_d17O_error"] = calculate_d17O_error(
    data["water_d18O"],
    data["water_Dp17O"],
    data["water_d18O_error"],
    data["water_Dp17O_error"],
)

# Calculate carbonate d17O error from d18O and Dp17O, retaining their
# correlation instead of using the independently reported d17O error.
data["d17O_AC_error_propagated"] = calculate_d17O_error(
    data["d18O_AC"],
    data["Dp17O_AC"],
    data["d18O_error"],
    data["Dp17O_error"],
)

data["18O_alpha"], data["18O_alpha_error"] = calculate_fractionation_factor(
    data["d18O_AC"],
    data["water_d18O"],
    data["d18O_error"],
    data["water_d18O_error"],
)
data["17O_alpha"], data["17O_alpha_error"] = calculate_fractionation_factor(
    data["d17O_AC"],
    data["water_d17O"],
    data["d17O_AC_error_propagated"],
    data["water_d17O_error"],
)

# Calculate temperature variables and theta ----------------------------------
data["1000_over_T"] = 1000.0 / (data["clumped_T"] + 273.15)
data["1000_over_T_error"] = (
    1000.0
    * data["clumped_T_err"]
    / (data["clumped_T"] + 273.15) ** 2
)
data["theta"] = np.log(data["17O_alpha"]) / np.log(data["18O_alpha"])

# Propagate theta directly from the independent carbonate and water d18O and
# Dp17O inputs. This retains the correlated d18O contribution to both alphas.
carbonate_d17O_offset = data["d17O_AC"] - d17O(
    data["d18O_AC"], data["Dp17O_AC"]
)
theta_inputs = {
    "carbonate_d18O": data["d18O_AC"],
    "carbonate_Dp17O": data["Dp17O_AC"],
    "water_d18O": data["water_d18O"],
    "water_Dp17O": data["water_Dp17O"],
}
theta_input_errors = {
    "carbonate_d18O": data["d18O_error"],
    "carbonate_Dp17O": data["Dp17O_error"],
    "water_d18O": data["water_d18O_error"],
    "water_Dp17O": data["water_Dp17O_error"],
}
theta_variance = np.zeros(len(data))
for input_name, input_error in theta_input_errors.items():
    inputs_plus = theta_inputs.copy()
    inputs_minus = theta_inputs.copy()
    inputs_plus[input_name] = theta_inputs[input_name] + input_error
    inputs_minus[input_name] = theta_inputs[input_name] - input_error
    theta_plus = calculate_theta_from_inputs(
        **inputs_plus,
        carbonate_d17O_offset=carbonate_d17O_offset,
    )
    theta_minus = calculate_theta_from_inputs(
        **inputs_minus,
        carbonate_d17O_offset=carbonate_d17O_offset,
    )
    theta_contribution = (theta_plus - theta_minus) / 2.0
    theta_variance += np.asarray(theta_contribution) ** 2
data["theta_error"] = np.sqrt(theta_variance)

# Fit theta with the theoretical infinite-temperature intercept fixed at 0.5305.
theta_intercept = 0.5305
inverse_temperature = 1.0 / (data["clumped_T"] + 273.15)
theta_ikaite_water_slope = np.sum(
    inverse_temperature * (data["theta"] - theta_intercept)
) / np.sum(inverse_temperature**2)
theta_residuals = (
    data["theta"]
    - theta_intercept
    - theta_ikaite_water_slope * inverse_temperature
)
theta_ikaite_water_slope_error = np.sqrt(
    np.sum(theta_residuals**2)
    / (len(data) - 1)
    / np.sum(inverse_temperature**2)
)
data["theta_ikaite_water_fit"] = (
    theta_ikaite_water_slope / (data["clumped_T"] + 273.15)
    + theta_intercept
)
print(
    "theta_ikaite_water = "
    f"{theta_ikaite_water_slope:.2f}/T + {theta_intercept:.4f}"
)
theta_slope_mathml_sign = "-" if theta_ikaite_water_slope < 0 else "+"
theta_mathml = (
    '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">'
    "<mrow>"
    "<msub><mi>&#x03B8;</mi><mtext>ikaite/water</mtext></msub>"
    "<mo>=</mo>"
    f"<mo>{theta_slope_mathml_sign}</mo>"
    f"<mn>{abs(theta_ikaite_water_slope):.2f}</mn>"
    "<mo>/</mo>"
    "<mi>T</mi>"
    "<mo>+</mo>"
    f"<mn>{theta_intercept:.4f}</mn>"
    "</mrow>"
    "</math>"
)
print(theta_mathml,"\n")

# Fit the ikaite 18O fractionation intercept while fixing the slope to the
# Daeron et al. (2019) value.
a18_ikaite_water_slope = 17.57
data["1000_ln_18O_alpha"] = 1000.0 * np.log(data["18O_alpha"])
data["1000_ln_18O_alpha_error"] = (
    1000.0 * data["18O_alpha_error"] / data["18O_alpha"]
)
data["1000_ln_18O_alpha_total_fit_error"] = np.hypot(
    data["1000_ln_18O_alpha_error"],
    a18_ikaite_water_slope * data["1000_over_T_error"],
)
a18_fit_weights = 1.0 / data["1000_ln_18O_alpha_total_fit_error"] ** 2
a18_ikaite_water_intercept = np.average(
    data["1000_ln_18O_alpha"]
    - a18_ikaite_water_slope * data["1000_over_T"],
    weights=a18_fit_weights,
)
a18_ikaite_water_intercept_error = np.sqrt(1.0 / a18_fit_weights.sum())
data["18O_alpha_ikaite_water_fit"] = np.exp(
    (
        a18_ikaite_water_slope
        * 1000.0
        / (data["clumped_T"] + 273.15)
        + a18_ikaite_water_intercept
    )
    / 1000.0
)
data["1000_ln_18O_alpha_ikaite_water_fit"] = 1000.0 * np.log(
    data["18O_alpha_ikaite_water_fit"]
)
a18_ikaite_water_intercept_sign = (
    "+" if a18_ikaite_water_intercept >= 0 else "-"
)
print(
    "1000 ln 18alpha = "
    f"{a18_ikaite_water_slope:.2f} * 10**3 / T "
    f"{a18_ikaite_water_intercept_sign} "
    f"{abs(a18_ikaite_water_intercept):.2f}"
)
a18_mathml = (
    '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">'
    "<mrow>"
    "<msup><mn>10</mn><mn>3</mn></msup>"
    "<mi>ln</mi>"
    "<msub><mmultiscripts><mi>&#x03B1;</mi>"
    "<mprescripts/><none/>"
    "<mrow><mn>18</mn><mo>/</mo><mn>16</mn></mrow>"
    "</mmultiscripts>"
    "<mtext>ikaite&#x005C;water</mtext></msub>"
    "<mo>=</mo>"
    f"<mn>{a18_ikaite_water_slope:.2f}</mn>"
    "<mo>&#x00D7;</mo>"
    "<msup><mn>10</mn><mn>3</mn></msup>"
    "<mo>/</mo>"
    "<mi>T</mi>"
    f"<mo>{a18_ikaite_water_intercept_sign}</mo>"
    f"<mn>{abs(a18_ikaite_water_intercept):.2f}</mn>"
    "</mrow>"
    "</math>"
)
print(a18_mathml,"\n")

# save to csv
data.to_csv(export_tables_dir / "modern_data_with_sw_values.csv", index=False)


# Plot ikaite-water fractionation --------------------------------------------
def inverse_temperature_to_celsius(inverse_temperature):
    """Convert 1000/T in K^-1 to temperature in degrees Celsius."""
    return 1000.0 / np.asarray(inverse_temperature) - 273.15


def celsius_to_inverse_temperature(temperature_C):
    """Convert temperature in degrees Celsius to 1000/T in K^-1."""
    return 1000.0 / (np.asarray(temperature_C) + 273.15)


fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True)

for ax in fig.get_axes():
    i = fig.get_axes().index(ax)
    ax.text(0.02, 0.98, chr(65 + i),
            size=12, weight="bold", ha="left", va="top",
            transform=ax.transAxes)

dark_gray = "#231F20"
light_gray = "#B3B4B4"
sample_styles = {
    "Congo fan": {"fc": dark_gray, "ec": dark_gray},
    "South Georgia": {"fc": light_gray, "ec": dark_gray},
    "Bransfield Strait": {"fc": dark_gray, "ec": light_gray},
    "Nankai Trough": {"fc": "white", "ec": "black"},
}

for _, sample in data.iterrows():
    style = sample_styles[sample["site"]]

    ax1.errorbar(
        sample["1000_over_T"],
        sample["1000_ln_18O_alpha"],
        xerr=sample["1000_over_T_error"],
        yerr=sample["1000_ln_18O_alpha_error"],
        fmt="none",
        ecolor="k",
        capsize=0,
        zorder=0,
    )
    ax1.scatter(
        sample["1000_over_T"],
        sample["1000_ln_18O_alpha"],
        fc=style["fc"],
        ec=style["ec"],
        marker="D",
        s=30,
        linewidths=1,
        label=sample["site"],
        zorder=1,
    )

    ax2.errorbar(
        sample["1000_over_T"],
        sample["theta"],
        xerr=sample["1000_over_T_error"],
        yerr=sample["theta_error"],
        fmt="none",
        ecolor="k",
        capsize=0,
        zorder=0,
    )
    ax2.scatter(
        sample["1000_over_T"],
        sample["theta"],
        fc=style["fc"],
        ec=style["ec"],
        marker="D",
        s=30,
        linewidths=1,
        zorder=1,
    )

site_handles, site_labels = ax1.get_legend_handles_labels()

x_min = (data["1000_over_T"] - data["1000_over_T_error"]).min()
x_max = (data["1000_over_T"] + data["1000_over_T_error"]).max()
x_padding = 0.05 * (x_max - x_min)
xlim = (x_min - x_padding, x_max + x_padding)
ax1.set_xlim(xlim)
ax2.set_xlim(xlim)

for ax in (ax1, ax2):
    temperature_axis = ax.secondary_xaxis(
        "top",
        functions=(
            inverse_temperature_to_celsius,
            celsius_to_inverse_temperature,
        ),
    )
    temperature_axis.set_xlabel("Temperature (°C)")

# Plot published fractionation relationships over the measured temperature
# range, with a small extension beyond the outermost data points.
model_1000_over_T = np.linspace(
    xlim[0] - 0.005,
    xlim[1] + 0.005,
    200,
)
model_temperature_K = 1000.0 / model_1000_over_T

# Extend the theta models to the infinite-temperature limit at 1000/T = 0.
theta_model_1000_over_T = np.linspace(
    xlim[0],
    xlim[1] + 0.005,
    500,
)
theta_model_temperature_K = np.divide(
    1000.0,
    theta_model_1000_over_T,
    out=np.full_like(theta_model_1000_over_T, np.inf),
    where=theta_model_1000_over_T != 0,
)

line_a18_KO97, = ax1.plot(
    model_1000_over_T,
    1000.0 * np.log(a18_KO97(model_temperature_K)),
    color=C["red"],
    ls = ":",
    label="Kim and O'Neil (1997)",
    zorder=-1,
)
line_a18_Daeron2019, = ax1.plot(
    model_1000_over_T,
    1000.0 * np.log(a18_Daeron2019(model_temperature_K)),
    color=C["blue"],
    ls="--",
    label="Daëron et al. (2019)",
    zorder=-1,
)
line_a18_ikaite_water, = ax1.plot(
    model_1000_over_T,
    a18_ikaite_water_slope * 1000.0 / model_temperature_K
    + a18_ikaite_water_intercept,
    color=C["black"],
    ls="-",
    lw = 1,
    label="this study",
    zorder=-1,
)
ax1.fill_between(
    model_1000_over_T,
    (
        a18_ikaite_water_slope * 1000.0 / model_temperature_K
        + a18_ikaite_water_intercept
        - a18_ikaite_water_intercept_error
    ),
    (
        a18_ikaite_water_slope * 1000.0 / model_temperature_K
        + a18_ikaite_water_intercept
        + a18_ikaite_water_intercept_error
    ),
    color=C["black"],
    alpha=0.15,
    linewidth=0,
    zorder=-2,
)

# Annotate the offset between Daeron et al. (2019) and this study at 1000/T = 3.66.
difference_x = 3.62
difference_temperature_K = 1000.0 / difference_x
difference_daeron = 1000.0 * np.log(a18_Daeron2019(difference_temperature_K))
difference_this_study = (
    a18_ikaite_water_slope * difference_x + a18_ikaite_water_intercept
)
difference_a18 = abs(difference_daeron - difference_this_study)
ax1.annotate(
    "",
    xy=(difference_x, difference_daeron),
    xytext=(difference_x, difference_this_study),
    arrowprops=dict(arrowstyle="<->", color="black", lw=0.7),
    zorder=3,
)
ax1.text(
    difference_x + 0.001,
    (difference_daeron + difference_this_study) / 1.99,
    f"∆ = {difference_a18:.1f}",
    rotation=90,
    ha="left",
    va="center",
    fontsize=6,
)

line_theta_H18, = ax2.plot(
    theta_model_1000_over_T,
    theta_H18(theta_model_temperature_K),
    color=C["red"],
    ls="--",
    label="Hayles et al. (2018)",
    zorder=-1,
)
line_theta_GZ19, = ax2.plot(
    theta_model_1000_over_T,
    theta_GZ19_cc(theta_model_temperature_K),
    color=C["yellow6"],
    ls="-.",
    label="Guo & Zhou (2019)",
    zorder=-1,
)
line_theta_ikaite_water, = ax2.plot(
    theta_model_1000_over_T,
    theta_ikaite_water_slope / theta_model_temperature_K + theta_intercept,
    color=C["black"],
    ls="-",
    lw = 1,
    label="this study",
    zorder=-1,
)
theta_ikaite_water_fit = (
    theta_ikaite_water_slope / theta_model_temperature_K + theta_intercept
)
theta_ikaite_water_confidence = (
    theta_ikaite_water_slope_error / theta_model_temperature_K
)
ax2.fill_between(
    theta_model_1000_over_T,
    theta_ikaite_water_fit - theta_ikaite_water_confidence,
    theta_ikaite_water_fit + theta_ikaite_water_confidence,
    color=C["black"],
    alpha=0.15,
    linewidth=0,
    zorder=-2,
)

ax1.set_xlabel(r"Temperature ($10^{3}\,K^{-1}$)")
ax1.set_ylabel(r"$10^{3}\,\ln\,{}^{18}\alpha_{\mathrm{ikaite \backslash water}}$")

add_standard_legend(
    ax1,
    handles=[line_a18_KO97, line_a18_Daeron2019, line_a18_ikaite_water],
    loc="lower right",
    fontsize=6
)

ax2.set_xlabel(r"Temperature ($10^{3}\,K^{-1}$)")
ax2.set_ylabel(r"$\theta_{\mathrm{ikaite \backslash water}}$")

add_standard_legend(
    ax2,
    handles=[
        line_theta_H18,
        line_theta_GZ19,
        line_theta_ikaite_water,
        # high_temperature_limit,
    ],
    loc="lower right",
    fontsize=6,
)

sample_legend = fig.legend(
    site_handles,
    site_labels,
    title="modern ikaites",
    alignment="left",
    loc="upper center",
    bbox_to_anchor=(0.5, -0.0),
    ncol=4,
    fontsize=6,
    frameon=True,
    fancybox=False,
    edgecolor="black",
    facecolor="white",
    framealpha=1.0,
)
sample_legend.get_title().set_fontweight("bold")
sample_legend.get_title().set_fontsize(6)

save_figure(fig, figures_dir, "ikaite_water_equilibrium")