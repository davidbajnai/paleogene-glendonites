# The following code is used to:
# Scale the reference gases using accepted IAEA-603 and NBS-18 isotope values
#
# INPUT:    Bajnai_2024_ChemGeol.csv
#
# OUTPUT:   Unscaled, measured, accepted, and scaled values in the terminal only


# Import libraries -------------------------------------------------------------
from pathlib import Path

import pandas as pd

from _helper_functions import d17o, dp17o, to_vsmow

# Define directories -----------------------------------------------------------
data_dir = Path(__file__).resolve().parent.parent / "data" / "input"

# Define additional functions --------------------------------------------------
def print_info(
    df: pd.DataFrame, d18o_col: str, dp17o_col: str, sample_name: str
) -> None:
    gas_subset = df[df["SampleName"].str.contains(sample_name)].copy()
    d18o_mean = gas_subset[d18o_col].mean()
    d18o_std = gas_subset[d18o_col].std()
    dp17o_mean = gas_subset[dp17o_col].mean()
    dp17o_std = gas_subset[dp17o_col].std()
    n_gas = len(gas_subset)
    print(
        f"{sample_name}, N = {n_gas}, "
        f"δ18O = {d18o_mean:.3f} (±{d18o_std:.3f}) ‰, "
        f"Δ’17O = {dp17o_mean:.0f} (±{dp17o_std:.0f}) ppm",
        end="",
    )


# Read data --------------------------------------------------------------------
df = pd.read_csv(data_dir / "Bajnai_2024_ChemGeol.csv")
df = df[df["measurementPeriod"].str.contains("2023-09-20 to 2023-10-08")]
df = df[df["SampleName"].str.contains("light|heavy|NBS18|IAEA")].copy()

print("\nNumber of replicates and standard deviations for the references:")
print_info(df, "d18O", "Dp17O", "light")
print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "heavy")
print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "IAEA603")
print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "NBS18")
print("\t<--- unscaled")

# Scale reference values -------------------------------------------------------
# Measured CO2 values
iaea603_d18o_measured = df[df["SampleName"].str.contains("IAEA603")]["d18O"].mean()
iaea603_d17o_measured = df[df["SampleName"].str.contains("IAEA603")]["d17O"].mean()
print("\nMeasured values:")
print(f"IAEA-603: δ18O = {iaea603_d18o_measured:.3f} ‰, δ17O = {iaea603_d17o_measured:.3f} ‰")

nbs18_d18o_measured = df[df["SampleName"].str.contains("NBS18")]["d18O"].mean()
nbs18_d17o_measured = df[df["SampleName"].str.contains("NBS18")]["d17O"].mean()
print(f"NBS-18: δ18O = {nbs18_d18o_measured:.3f} ‰, δ17O = {nbs18_d17o_measured:.3f} ‰")

# Accepted CO2 values from the original Wostbrock et al. (2020) calibration.
# These retain the 25 °C acid-digestion basis of the input measurements.
iaea603_d18o_accepted = (to_vsmow(-2.37) + 1000) * 1.01025 - 1000 # 25 °C acid digestion
# iaea603_dp17o_accepted = -147 # from Wostbrock et al. (2020)
iaea603_dp17o_accepted = -127.3-4.1 # from Chaillot et al. (2022) # from 90 °C acid digestion, corrected for a 25 °C reaction
iaea603_d17o_accepted = d17o(iaea603_d18o_accepted, iaea603_dp17o_accepted)
print("\nAccepted values:")
print(f"IAEA-603: δ18O = {iaea603_d18o_accepted:.3f} ‰, Δ’17O = {iaea603_dp17o_accepted:.0f} ppm")

nbs18_d18o_accepted = (to_vsmow(-23.01) + 1000) * 1.01025 - 1000 # 25 °C acid digestion
# nbs18_dp17o_accepted = -100 # from Wostbrock et al. (2020)
nbs18_dp17o_accepted = -101.3-4.1 # from Chaillot et al. (2022) # from 90 °C acid digestion, corrected for a 25 °C reaction
nbs18_d17o_accepted = d17o(nbs18_d18o_accepted, nbs18_dp17o_accepted)
print(f"NBS-18: δ18O = {nbs18_d18o_accepted:.3f} ‰, Δ’17O = {nbs18_dp17o_accepted:.0f} ppm")

# Calculate the scaling factors
slope_d18o = (nbs18_d18o_accepted - iaea603_d18o_accepted) / (nbs18_d18o_measured - iaea603_d18o_measured)
intercept_d18o = iaea603_d18o_accepted - slope_d18o * iaea603_d18o_measured
slope_d17o = (nbs18_d17o_accepted - iaea603_d17o_accepted) / (nbs18_d17o_measured - iaea603_d17o_measured)
intercept_d17o = iaea603_d17o_accepted - slope_d17o * iaea603_d17o_measured

# Scale the measured values
df["d18O_scaled"] = slope_d18o * df["d18O"] + intercept_d18o
df["d17O_scaled"] = slope_d17o * df["d17O"] + intercept_d17o
df["Dp17O_scaled"] = dp17o(df["d17O_scaled"], df["d18O_scaled"])

# Average reference gas replicates ---------------------------------------------
df_gases = df.loc[:, ["SampleName", "d18O_scaled", "Dp17O_scaled"]]
df_gases = df_gases[df_gases["SampleName"].str.contains("light|heavy")]
df_gases_mean = df_gases.groupby('SampleName').mean().reset_index()
df_gases_mean = df_gases_mean.rename(columns={'d18O_scaled': 'd18O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})
df_gases_std = df_gases.groupby('SampleName').std().reset_index()
df_gases_std = df_gases_std.rename(columns={'d18O_scaled': 'd18O_SD', 'Dp17O_scaled': 'Dp17O_SD'})

print("\nScaled reference gas values (replicate mean ± 1 SD):")
df_gases_summary = df_gases_mean.merge(df_gases_std, on="SampleName")
print(df_gases_summary.to_string(
    index=False,
    formatters={
        "d18O_CO2": "{:.3f}".format,
        "Dp17O_CO2": "{:.0f}".format,
        "d18O_SD": "{:.3f}".format,
        "Dp17O_SD": "{:.0f}".format,
    },
))
print("δ18O and its SD are in ‰; Δ’17O and its SD are in ppm.")
