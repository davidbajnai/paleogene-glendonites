# Import libraries ------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
from pathlib import Path
from _helper_functions import *
from _plot_style import *

# Define directories ----------------------------------------------------------
figures_dir = Path(__file__).resolve().parent.parent / "figures"
data_dir = Path(__file__).resolve().parent.parent / "source_data"


# Apply plot style ------------------------------------------------------------
apply_plot_style()
plt.rcParams["figure.figsize"] = (9, 4)
plt.rcParams.update({'font.size': 7})


# Define additional functions -------------------------------------------------

# Function to print info for scaled samples
def print_info(df, d18O_col, Dp17O_col, sample_name):
    gas_subset = df[df["SampleName"].str.contains(sample_name)].copy()
    d18O_mean = gas_subset[d18O_col].mean()
    d18O_std = gas_subset[d18O_col].std()
    Dp17O_mean = gas_subset[Dp17O_col].mean()
    Dp17O_std = gas_subset[Dp17O_col].std()
    N_gas = len(gas_subset)
    print(f"{sample_name}, N = {N_gas}, d18O = {d18O_mean:.3f}(±{d18O_std:.3f})‰, ∆'17O = {Dp17O_mean:.0f}(±{Dp17O_std:.0f}) ppm", end="")

# This function scales the data based on the accepted values of the light and heavy reference gases
# The scaling is done for each measurement period separately
def scaleData(df):

    df["dateTimeMeasured"] = pd.to_datetime(df["dateTimeMeasured"])

    # Perform the scaling for each meausrement period separately
    df_samples = pd.DataFrame()
    grouped = df.groupby("measurementPeriod")
    for FigNum, (period, group) in enumerate(grouped):

        print(f"\nMeasurement period {period}:")
        print_info(group, "d18O", "Dp17O", "light"); print("\t<--- unscaled")
        print_info(group, "d18O", "Dp17O", "heavy"); print("\t<--- unscaled")

        # Do the scaling here, based on the accepted values of the light and heavy reference gases

    # Measured CO2 values
        heavy_d18O_measured = group[group["SampleName"].str.contains("heavy")]["d18O"].mean()
        heavy_d17O_measured = group[group["SampleName"].str.contains("heavy")]["d17O"].mean()

        light_d18O_measured = group[group["SampleName"].str.contains("light")]["d18O"].mean()
        light_d17O_measured = group[group["SampleName"].str.contains("light")]["d17O"].mean()

        # Accepted CO2 values
        # Linked to 90 °C acid reaction of
        # NBS-18 and IAEA-603 in Chaillot et al. (2025)
        heavy_d18O_accepted = 76.475
        heavy_Dp17O_accepted = -165
        heavy_d17O_accepted = d17O(heavy_d18O_accepted, heavy_Dp17O_accepted)

        light_d18O_accepted = -1.139
        light_Dp17O_accepted = -167
        light_d17O_accepted = d17O(light_d18O_accepted, light_Dp17O_accepted)

        # Calculate the scaling factors
        slope_d18O = (light_d18O_accepted - heavy_d18O_accepted) / (light_d18O_measured - heavy_d18O_measured)
        intercept_d18O = heavy_d18O_accepted - slope_d18O * heavy_d18O_measured

        slope_d17O = (light_d17O_accepted - heavy_d17O_accepted) / (light_d17O_measured - heavy_d17O_measured)
        intercept_d17O = heavy_d17O_accepted - slope_d17O * heavy_d17O_measured

        # Scale the measured values
        group["d18O_scaled"] = slope_d18O*group['d18O']+intercept_d18O
        group["d17O_scaled"] = slope_d17O*group['d17O']+intercept_d17O
        group["Dp17O_scaled"] = Dp17O(group["d17O_scaled"], group["d18O_scaled"])

        # Assign colors and markers to samples
        categories = group["SampleName"].unique()
        markers = dict(zip(categories, ["o", "s", "D", "v", "^",
                                        "<", ">", "p", "P", "*"]*10))
        colors = dict(zip(categories, plt.cm.tab20(np.linspace(0, 1, 20))))

        # Figure: unscaled Dp17O vs time

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 5))

        for cat in categories:
            data = group[group["SampleName"] == cat]
            ax1.scatter(data["dateTimeMeasured"], data["Dp17O"],
                       marker=markers[cat], fc=colors[cat], label=cat)
            if np.isnan(data["Dp17OError"]).any() == False:
                ax1.errorbar(data["dateTimeMeasured"], data["Dp17O"],
                             yerr=data["Dp17OError"],
                             fmt="none", color="#cacaca", zorder=0)
            ax2.scatter(data["dateTimeMeasured"], data["pCO2Mismatch"],
                        marker=markers[cat], fc=colors[cat])
            if np.isnan(data["pCO2Mismatch_error"]).any() == False:
                ax2.errorbar(data["dateTimeMeasured"], data["pCO2Mismatch"],
                             yerr=data["pCO2Mismatch_error"],
                             fmt="none", color="#cacaca", zorder=0)
        
        ax1.legend(loc='upper right', bbox_to_anchor=(1.15, 1))

        ax1.set_title(f"Measurement period: {period}")
        ax1.set_ylabel(r"$\Delta\prime^{17}$O (ppm, unscaled CO$_2$)")
        ax2.set_ylabel(r'$\chi\prime^{smp}_{626}$ - $\chi\prime^{ref}_{626}$ (µmol mol$^{-1}$)')
        ax2.set_xlabel("Measurement date")

        setup_axes(fig)
        save_figure(fig, figures_dir, f"OFFLINE_replicates_{FigNum+1}")

        # Exclude the standards from the exported dataframe
        group = group[~group["SampleName"].str.contains("NBS18|IAEA603|light|heavy|DH11|NBS")].copy()
        df_samples = pd.concat([df_samples, group])

    return df_samples


# This function averages the scaled data from multiple measurement periods
def average_data(df):
    # Calculate the mean values from the replicate measurements
    df = df.loc[:, ["SampleName", "d18O_scaled", "d17O_scaled", "Dp17O_scaled"]]
    df_mean = df.groupby('SampleName').mean().reset_index()
    df_mean = df_mean.rename(columns={
                             'd18O_scaled': 'd18O_CO2', 'd17O_scaled': 'd17O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})

    # Calculate reproducibility from replicate measurements
    df_std = df.groupby('SampleName').std().reset_index()
    df_std = df_std.rename(columns={'d18O_scaled': 'd18O_error',
                           'd17O_scaled': 'd17O_error', 'Dp17O_scaled': 'Dp17O_error'})
    df_std['Replicates'] = df.groupby(
        'SampleName').size().reset_index(name='counts')['counts']

    confidence_level = 0.68
    alpha = 1 - confidence_level
    dof = df_std['Replicates'] - 1
    t_crit = pd.Series(np.nan, index=df_std.index, dtype=float)
    valid = dof > 0
    t_crit.loc[valid] = t.ppf(1 - alpha / 2, dof.loc[valid])

    for col in ['d18O_error', 'd17O_error', 'Dp17O_error']:
        standard_error = df_std[col] / np.sqrt(df_std['Replicates'])
        df_std[col] = standard_error * t_crit

    dfMerged = df_mean.merge(df_std, on='SampleName')
    df = dfMerged

    return df


# Data processing -------------------------------------------------------------

# Scale the data
df = scaleData(pd.read_csv(data_dir /  "replicate_data_OFFLINE.csv"))

# Export the scaled data on a replicate level
export_table(df, "replicate_data_OFFLINE_scaled")

# Average the data
df_avg = average_data(df)

print(df_avg)

# Apply acid fractionation factor
df_avg[["d18O_AC", "d17O_AC", "Dp17O_AC"]] = df_avg.apply(lambda x: applyAFF(
    x["d18O_CO2"], x["d17O_CO2"], "calcite", 25), axis=1, result_type="expand")

print(df_avg)

# Merge data with sample info
info = pd.read_csv(data_dir / "sample_info.csv")
df_avg = df_avg.merge(info, on='SampleName')


# Export the averaged data
export_table(df_avg,
             "replicate_data_OFFLINE_averaged",
             print_columns=["SampleName", "d18O_AC",
                            "Dp17O_AC",  "Dp17O_error", "Replicates"]
             )


# Figures ---------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2)

# Assign colors and markers to samples
df.sort_values(by="SampleName", inplace=True)
categories = df["SampleName"].unique()
markers = dict(zip(categories, ["o", "s", "D", "^", "v", "X", "P", "*"]*10))
colors = dict(zip(categories, plt.cm.tab20(
    np.linspace(0, 1, len(categories)))))

# Subplot a: scaled sample replicates
for cat in categories:
    data = df[df["SampleName"] == cat]
    ax1.scatter(prime(data["d18O_scaled"]), data["Dp17O_scaled"],
                marker=markers[cat], fc=colors[cat], label=cat)
    ax1.errorbar(prime(df["d18O_scaled"]), df["Dp17O_scaled"],
                 yerr=df["Dp17OError"],
                 xerr=df["d18OError"],
                 fmt="none", color="#cacaca", zorder=0)

# Axis properties
ylim = ax1.get_ylim()
xlim = ax1.get_xlim()
ax1.set_ylabel(r"$\Delta\prime^{17}$O (ppm, CO$_2$)")
ax1.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW, CO$_2$)")
ax1.text(0.02, 0.98, "e", size=14, ha="left", va="top",
         transform=ax1.transAxes, fontweight="bold")

# Subplot B: scaled sample averages
for cat in categories:
    data = df_avg[df_avg["SampleName"] == cat]
    ax2.scatter(prime(data["d18O_CO2"]), data["Dp17O_CO2"],
                marker=markers[cat], fc=colors[cat], label=cat)
    ax2.errorbar(prime(df_avg["d18O_CO2"]), df_avg["Dp17O_CO2"],
                 yerr=df_avg["Dp17O_error"],
                 xerr=df_avg["d18O_error"],
                 fmt="none", color="#cacaca", zorder=0)

ax2.legend(loc='upper right', bbox_to_anchor=(1.36, 1))

# Axis properties
ax2.set_ylim(ylim)
ax2.set_xlim(xlim)
ax2.set_ylabel(r"$\Delta\prime^{17}$O (ppm, CO$_2$)")
ax2.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW, CO$_2$)")
ax2.text(0.02, 0.98, "f", size=14, ha="left", va="top",
         transform=ax2.transAxes, fontweight="bold")

save_figure(fig, figures_dir, f"OFFLINE_replicates_summary")


combine_images(
    "OFFLINE_replicates",
    ["OFFLINE_replicates_1",
     "OFFLINE_replicates_2",
     "OFFLINE_replicates_3",
     "OFFLINE_replicates_4"],
    remove_originals=True
)