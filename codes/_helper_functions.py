from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


def to_VSMOW(VPDB):
    return (VPDB*1.03092)+30.92


def prime(x):
    return 1000 * np.log(x / 1000 + 1)


def unprime(x):
    return (np.exp(x / 1000) - 1) * 1000


def Dp17O(d17O, d18O):
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


def d17O(d18O, Dp17O):
    return unprime(Dp17O / 1000 + 0.528 * prime(d18O))


def mix_d17O(d18O_A, d17O_A=None, D17O_A=None, d18O_B=None, d17O_B=None, D17O_B=None, step=100):
    ratio_B = np.arange(0, 1+1/step, 1/step)

    if d17O_A is None:
        d17O_A = unprime(D17O_A/1000 + 0.528 * prime(d18O_A))

    if d17O_B is None:
        d17O_B = unprime(D17O_B/1000 + 0.528 * prime(d18O_B))

    mix_d18O = ratio_B * float(d18O_B) + (1 - ratio_B) * float(d18O_A)
    mix_d17O = ratio_B * float(d17O_B) + (1 - ratio_B) * float(d17O_A)
    mix_D17O = Dp17O(mix_d17O, mix_d18O)
    xB = ratio_B * 100

    df = pd.DataFrame(
        {'mix_d17O': mix_d17O, 'mix_d18O': mix_d18O, 'mix_Dp17O': mix_D17O, 'xB': xB})
    return df


def a18_cc(T):

    # Kim et O'Neil (1997) - calcite
    return np.exp((18.04 * 1000 / T - 32.18) / 1000)

    # Alternative equations:

    # Daeron et al. (2019) - calcite
    # return np.exp((17.57 * 1000 / T - 29.13) / 1000)

    # Guo and Zhou (2019) – aragonite
    # return 0.0201 * (1000 / T) + 0.9642

    # Hayles et al. (2018) - calcite
    # B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    # B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    # return np.exp(B_calcite) / np.exp(B_water)


def theta_cc(T):

    # Wostbrock et al. (2020) - aragonite
    # return -1.53 / T + 0.5305

    # Alternative equations:

    # Wostbrock et al. (2020) - calcite
    # return -1.39 / T + 0.5305

    # Guo and Zhou (2019) – aragonite
    return 59.1047/T**2 + -1.4089/T + 0.5297

    # Hayles et al. (2018) - calcite
    # K_calcite = 1.019124E+09 / T**5 + -2.117501E+07 / T**4 + 1.686453E+05 / T**3 + -5.784679E+02 / T**2 + 1.489666E-01 / T + 0.5304852
    # B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    # K_water = 7.625734E+06 / T**5 + 1.216102E+06 / T**4 + -2.135774E+04 / T**3 + 1.323782E+02 / T**2 + -4.931630E-01 / T + 0.5306551
    # B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    # a18 = np.exp(B_calcite) / np.exp(B_water)
    # return K_calcite + (K_calcite-K_water) * (B_water / np.log(a18))


def a17_cc(T):
    return a18_cc(T)**theta_cc(T)


def d18O_cc(equilibrium_temperatures, d18Ow):
    return a18_cc(equilibrium_temperatures) * (d18Ow+1000) - 1000


def d17O_cc(equilibrium_temperatures, d18Ow):
    return a17_cc(equilibrium_temperatures) * (d18Ow+1000) - 1000



def applyAFF(d18O_CO2, d17O_CO2, mineral, temperature_C):

    mineral = mineral.lower()

    # Acid fractionation factors at 90 C
    if mineral == "calcite":
        a18_90 = 1.008146
        a18_25 = 1.01025
        a18_70 = 1.008709

    elif mineral == "aragonite":
        a18_25 = 1.010596
        a18_90 = 1.008541
        a18_70 = 1.00909
    
    elif mineral == "dolomite":
        # Rosenbaum and Sheppard (1986)
        a18_90 = np.exp((6.65*10**5*(1/(90+273.15)**2)+4.23)/1000)
        a18_25 = np.mean([1.01178, 1.01186])
    else:
        raise ValueError(f"Unsupported mineral: {mineral}")

    if temperature_C == 90:
        a18 = a18_90
        theta_acid = 1.2019/(90+273.15) + 0.5192
    elif temperature_C == 25:
        a18 = a18_25
        theta_acid = 1.2019/(25+273.15) + 0.5192
    elif temperature_C == 70:
        a18 = a18_70
        theta_acid = 1.2019/(70+273.15) + 0.5192

    else:
        raise ValueError("temperature_C must be 25, 70, or 90.")

    d18O_AC = (d18O_CO2 + 1000) / a18 - 1000
    d17O_AC = (d17O_CO2 + 1000) / (a18 ** theta_acid) - 1000
    Dp17O_AC = Dp17O(d17O_AC, d18O_AC)

    return d18O_AC, d17O_AC, Dp17O_AC


def confidence_ellipse(x, y, ax, n_std=3, **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns:
    --------
    ell_radius_x, ell_radius_y : float
        Radii of the ellipse before scaling (raw covariance radii).
    scale_x, scale_y : float
        Scaled semi-major and semi-minor axes (after applying standard deviations).
    mean_x, mean_y : float
        Mean of x and y (center of the ellipse).
    transform : matplotlib.transforms.Affine2D
        Full affine transformation applied to the ellipse.
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])

    # Ellipse radii before scaling (raw covariance ellipse)
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)

    # Create the ellipse
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      **kwargs)

    # Standard deviations and scaling
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    # Apply the transformation (rotation, scaling, translation)
    rotation = np.deg2rad(45)  # Rotation in radians
    transf = transforms.Affine2D() \
        .rotate(rotation) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)

    return ell_radius_x, ell_radius_y, scale_x, scale_y, mean_x, mean_y, transf


def export_table(df, file_name, print_columns=None, include_index=False):
    """Export a table DataFrame as both CSV and XLSX, and optionally print columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Table content to export.
    file_name : str
        Output file name stem (without extension).
    print_columns : list[str] | bool | None
        Optional columns to print with tabulate.
        If True, all columns are printed.
        If None or False, no table is printed.
    include_index : bool
        Whether to include DataFrame index in exported files.

    Returns
    -------
    tuple[str, str]
        Paths to the exported CSV and XLSX files.
    """
    import os
    from tabulate import tabulate

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame")

    tables_dir = data_dir = Path(__file__).resolve().parent.parent / "export_tables"

    os.makedirs(os.path.join(data_dir, "csv"), exist_ok=True)
    os.makedirs(os.path.join(data_dir, "xlsx"), exist_ok=True)

    csv_path = os.path.join(tables_dir, "csv", f"{file_name}.csv")
    xlsx_path = os.path.join(tables_dir, "xlsx", f"{file_name}.xlsx")

    df.to_csv(csv_path, index=include_index)
    df.to_excel(xlsx_path, index=include_index)

    if print_columns is True:
        printable_df = df
    elif isinstance(print_columns, list):
        missing_cols = [col for col in print_columns if col not in df.columns]
        if missing_cols:
            raise KeyError(f"Columns not found in DataFrame: {missing_cols}")
        printable_df = df[print_columns]
    elif print_columns in (None, False):
        printable_df = None
    else:
        raise TypeError("print_columns must be a list of column names, True, False, or None")

    if printable_df is not None:
        print(tabulate(printable_df, headers="keys", tablefmt="simple_outline", showindex=False))
        print(f"Number of rows: {len(printable_df)}")

    return csv_path, xlsx_path