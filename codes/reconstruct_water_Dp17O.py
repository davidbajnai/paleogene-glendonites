import numpy as np
import pandas as pd

from _helper_functions import Dp17O, a17_cc, a18_cc, d17O, prime, theta_cc


def reconstruct_water_Dp17O(
    d18O_carbonate,
    Dp17O_carbonate,
    temperature_C,
    d18O_error,
    Dp17O_error,
    temperature_error_C,
):
    """Reconstruct parent-water values and corner-propagated rectangle bounds."""
    error_args = [d18O_error, Dp17O_error, temperature_error_C]
    if any(value is None for value in error_args):
        raise ValueError(
            "d18O_error, Dp17O_error, and temperature_error_C must be provided together."
        )

    index = None
    for value in [d18O_carbonate, Dp17O_carbonate, d18O_error, Dp17O_error]:
        if isinstance(value, pd.Series):
            index = value.index
            break

    d18O_carbonate, Dp17O_carbonate, temperature_C, d18O_error, Dp17O_error, temperature_error_C = np.broadcast_arrays(
        d18O_carbonate,
        Dp17O_carbonate,
        temperature_C,
        d18O_error,
        Dp17O_error,
        temperature_error_C,
    )

    temperature_K = temperature_C + 273.15
    a18 = a18_cc(temperature_K)
    theta = theta_cc(temperature_K)
    a17 = a17_cc(temperature_K)
    if not np.allclose(a17, a18**theta):
        raise ValueError("a17_cc(T) does not match a18_cc(T) ** theta_cc(T).")

    d18O_water = (d18O_carbonate + 1000) / a18 - 1000
    d17O_carbonate = d17O(d18O_carbonate, Dp17O_carbonate)
    d17O_water = (d17O_carbonate + 1000) / a17 - 1000
    Dp17O_water = Dp17O(d17O_water, d18O_water)

    rectangle_rows = []
    for d18O_c, Dp17O_c, temp_c, d18O_err, Dp17O_err, temp_err in zip(
        d18O_carbonate.flat,
        Dp17O_carbonate.flat,
        temperature_C.flat,
        d18O_error.flat,
        Dp17O_error.flat,
        temperature_error_C.flat,
    ):
        corner_d18O_water = []
        corner_Dp17O_water = []
        for d18O_value in [d18O_c - d18O_err, d18O_c + d18O_err]:
            for Dp17O_value in [Dp17O_c - Dp17O_err, Dp17O_c + Dp17O_err]:
                for temperature_value in [temp_c - temp_err, temp_c + temp_err]:
                    corner_temperature_K = temperature_value + 273.15
                    corner_a18 = a18_cc(corner_temperature_K)
                    corner_a17 = a17_cc(corner_temperature_K)
                    corner_d18O = (d18O_value + 1000) / corner_a18 - 1000
                    corner_d17O_carbonate = d17O(d18O_value, Dp17O_value)
                    corner_d17O = (corner_d17O_carbonate + 1000) / corner_a17 - 1000
                    corner_Dp17O = Dp17O(corner_d17O, corner_d18O)
                    corner_d18O_water.append(corner_d18O)
                    corner_Dp17O_water.append(corner_Dp17O)

        x_values = prime(np.array(corner_d18O_water))
        y_values = np.array(corner_Dp17O_water)
        rectangle_rows.append({
            "reconstructed_water_d18O_prime_min": x_values.min(),
            "reconstructed_water_d18O_prime_max": x_values.max(),
            "reconstructed_water_Dp17O_min": y_values.min(),
            "reconstructed_water_Dp17O_max": y_values.max(),
        })

    df = pd.DataFrame(rectangle_rows, index=index)
    df["reconstructed_water_d18O"] = np.asarray(d18O_water).ravel()
    df["reconstructed_water_Dp17O"] = np.asarray(Dp17O_water).ravel()

    return df
