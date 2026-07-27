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
source_data_dir = Path(__file__).resolve().parent.parent / "source_data"
export_tables_dir = Path(__file__).resolve().parent.parent / "export_tables" / "csv"


# Apply plot style ------------------------------------------------------------
apply_plot_style()
plt.rcParams["figure.figsize"] = (8, 3.5)


# This function averages the scaled data from multiple measurement periods
def average_data(df):
    # Calculate the mean values from the replicate measurements
    df = df.loc[:, ["SampleName", "d18O_AC", "d17O_AC", "Dp17O_AC"]]
    df_mean = df.groupby('SampleName').mean().reset_index()
    # df_mean = df_mean.rename(columns={
    #                          'd18O_scaled': 'd18O_CO2', 'd17O_scaled': 'd17O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})

    # Calculate reproducibility from replicate measurements
    df_std = df.groupby('SampleName').std().reset_index()
    df_std = df_std.rename(columns={'d18O_AC': 'd18O_error',
                           'd17O_AC': 'd17O_error', 'Dp17O_AC': 'Dp17O_error'})
    df_std['replicates'] = df.groupby(
        'SampleName').size().reset_index(name='counts')['counts']

    confidence_level = 0.95
    alpha = 1 - confidence_level
    dof = df_std['replicates'] - 1
    t_crit = pd.Series(np.nan, index=df_std.index, dtype=float)
    valid = dof > 0
    t_crit.loc[valid] = t.ppf(1 - alpha / 2, dof.loc[valid])

    for col in ['d18O_error', 'd17O_error', 'Dp17O_error']:
        standard_error = df_std[col] / np.sqrt(df_std['replicates'])
        df_std[col] = standard_error * t_crit

    dfMerged = df_mean.merge(df_std, on='SampleName')
    df = dfMerged

    return df


df_KIEL = pd.read_csv(export_tables_dir / "replicate_data_KIEL_scaled.csv")
df_OFFLINE = pd.read_csv(export_tables_dir / "replicate_data_OFFLINE_scaled.csv")

for dfs in [df_KIEL, df_OFFLINE]:
    dfs.drop(dfs[dfs["SampleName"] == "MLVGL75"].index, inplace=True)


# Apply acid fractionation factor
df_KIEL[["d18O_AC", "d17O_AC", "Dp17O_AC"]] = df_KIEL.apply(lambda x: applyAFF(
    x["d18O_scaled"], x["d17O_scaled"], "calcite", 90), axis=1, result_type="expand")

df_OFFLINE[["d18O_AC", "d17O_AC", "Dp17O_AC"]] = df_OFFLINE.apply(lambda x: applyAFF(
    x["d18O_scaled"], x["d17O_scaled"], "calcite", 25), axis=1, result_type="expand")


df = pd.concat([df_KIEL, df_OFFLINE], ignore_index=True)
df = df[["SampleName", "d18O_AC", "d17O_AC", "Dp17O_AC"]].copy()
df["SampleName"] = df["SampleName"].str.replace(r"[_\-.]", "", regex=True)
df_avg = average_data(df)

# Merge data with sample info
info = pd.read_csv(source_data_dir / "sample_info.csv")
info["SampleName"] = info["SampleName"].str.replace(r"[_\-.]", "", regex=True)
df_avg = df_avg.merge(info, on='SampleName')

export_table(df_avg,
             "combined_replicate_data",
             print_columns=["SampleName", "d18O_AC",
                            "Dp17O_AC",  "Dp17O_error", "replicates"]
             )

# Figures ---------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2)
setup_axes(fig)

# Assign colors and markers to samples
df.sort_values(by="SampleName", inplace=True)
categories = df["SampleName"].unique()
markers = dict(zip(categories, ["o", "s", "D", "^", "v", "X", "P", "*"]*10))
colors = dict(zip(categories, plt.cm.tab20(
    np.linspace(0, 1, len(categories)))))


# Subplot a: scaled sample replicates
for cat in categories:
    data = df[df["SampleName"] == cat]
    ax1.scatter(prime(data["d18O_AC"]), data["Dp17O_AC"],
                marker=markers[cat], fc=colors[cat], label=cat)
    # ax1.errorbar(prime(df["d18O_scaled"]), df["Dp17O_scaled"],
    #              yerr=df["Dp17O_error"],
    #              xerr=df["d18O_error"],
    #              fmt="none", color="#cacaca", zorder=0)


# Axis properties
ylim = ax1.get_ylim()
xlim = ax1.get_xlim()
ax1.set_ylabel(r"$\Delta\prime^{17}$O (ppm, CO$_2$)")
ax1.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW, CO$_2$)")

# Subplot B: scaled sample averages
for cat in categories:
    data = df_avg[df_avg["SampleName"] == cat]
    ax2.scatter(prime(data["d18O_AC"]), data["Dp17O_AC"],
                marker=markers[cat], fc=colors[cat], label=cat)
    ax2.errorbar(prime(df_avg["d18O_AC"]), df_avg["Dp17O_AC"],
                 yerr=df_avg["Dp17O_error"],
                 xerr=df_avg["d18O_error"],
                 fmt="none", color="#cacaca", zorder=0)

add_standard_legend(ax2, loc='upper right', bbox_to_anchor=(1.5, 1.02))

# Axis properties
ax2.set_ylim(ylim)
ax2.set_xlim(xlim)
ax2.set_ylabel(r"$\Delta\prime^{17}$O (ppm)")
ax2.set_xlabel(r"$\delta\prime^{18}$O (‰, VSMOW)")

save_figure(fig, figures_dir, f"combined_replicates")

# Print average reproducibilty of the measurements:
print(
    f"Average reproducibility (∆'17O): {df_avg['Dp17O_error'].mean():0.1f} ppm, d18O: {df_avg['d18O_error'].mean():0.2f} ‰")
