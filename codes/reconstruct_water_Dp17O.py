# The following code is used to:
# Calculate parent-water compositions and reconstruction uncertainties
#
# INPUT:    Carbonate isotope values, temperatures, and uncertainties
#
# OUTPUT:   Reconstructed water values and uncertainty bounds (returned DataFrame)

import numpy as np
from numpy.typing import ArrayLike
import pandas as pd

from _helper_functions import dp17o, a17_cc, a18_cc, d17o, prime, theta_cc


def reconstruct_water_dp17o(
    d18o_carbonate: ArrayLike,
    dp17o_carbonate: ArrayLike,
    temperature_c: ArrayLike,
    d18o_error: ArrayLike,
    dp17o_error: ArrayLike,
    temperature_error_c: ArrayLike,
) -> pd.DataFrame:
    """Reconstruct parent-water values and corner-propagated rectangle bounds."""
    error_args = [d18o_error, dp17o_error, temperature_error_c]
    if any(value is None for value in error_args):
        raise ValueError(
            "d18O_error, Dp17O_error, and temperature_error_C must be provided together."
        )

    index = None
    for value in [d18o_carbonate, dp17o_carbonate, d18o_error, dp17o_error]:
        if isinstance(value, pd.Series):
            index = value.index
            break

    d18o_carbonate, dp17o_carbonate, temperature_c, d18o_error, dp17o_error, temperature_error_c = np.broadcast_arrays(
        d18o_carbonate,
        dp17o_carbonate,
        temperature_c,
        d18o_error,
        dp17o_error,
        temperature_error_c,
    )

    temperature_k = temperature_c + 273.15
    a18 = a18_cc(temperature_k)
    theta = theta_cc(temperature_k)
    a17 = a17_cc(temperature_k)
    if not np.allclose(a17, a18**theta):
        raise ValueError("a17_cc(T) does not match a18_cc(T) ** theta_cc(T).")

    d18o_water = (d18o_carbonate + 1000) / a18 - 1000
    d17o_carbonate = d17o(d18o_carbonate, dp17o_carbonate)
    d17o_water = (d17o_carbonate + 1000) / a17 - 1000
    dp17o_water = dp17o(d17o_water, d18o_water)

    rectangle_rows = []
    for d18o_c, dp17o_c, temp_c, d18o_err, dp17o_err, temp_err in zip(
        d18o_carbonate.flat,
        dp17o_carbonate.flat,
        temperature_c.flat,
        d18o_error.flat,
        dp17o_error.flat,
        temperature_error_c.flat,
    ):
        corner_d18o_water = []
        corner_dp17o_water = []
        for d18o_value in [d18o_c - d18o_err, d18o_c + d18o_err]:
            for dp17o_value in [dp17o_c - dp17o_err, dp17o_c + dp17o_err]:
                for temperature_value in [temp_c - temp_err, temp_c + temp_err]:
                    corner_temperature_k = temperature_value + 273.15
                    corner_a18 = a18_cc(corner_temperature_k)
                    corner_a17 = a17_cc(corner_temperature_k)
                    corner_d18o = (d18o_value + 1000) / corner_a18 - 1000
                    corner_d17o_carbonate = d17o(d18o_value, dp17o_value)
                    corner_d17o = (corner_d17o_carbonate + 1000) / corner_a17 - 1000
                    corner_dp17o = dp17o(corner_d17o, corner_d18o)
                    corner_d18o_water.append(corner_d18o)
                    corner_dp17o_water.append(corner_dp17o)

        x_values = prime(np.array(corner_d18o_water))
        y_values = np.array(corner_dp17o_water)
        rectangle_rows.append({
            "reconstructed_water_d18O_prime_min": x_values.min(),
            "reconstructed_water_d18O_prime_max": x_values.max(),
            "reconstructed_water_Dp17O_min": y_values.min(),
            "reconstructed_water_Dp17O_max": y_values.max(),
        })

    df = pd.DataFrame(rectangle_rows, index=index)
    df["reconstructed_water_d18O"] = np.asarray(d18o_water).ravel()
    df["reconstructed_water_Dp17O"] = np.asarray(dp17o_water).ravel()

    return df
