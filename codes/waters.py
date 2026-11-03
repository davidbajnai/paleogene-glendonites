# The following code is used to:
# Calculate water endmembers and plot literature water isotope data
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
# OUTPUT:   waters.png
#           waters.svg

from pathlib import Path

from _plot_style import apply_plot_style
from _water_functions import _build_water_state, print_water_summary, plot_water_compositions

# Process water data -----------------------------------------------------------
apply_plot_style()
figures_dir = Path(__file__).resolve().parent.parent / "figures"
state = _build_water_state()
print_water_summary(state)
plot_water_compositions(state, figures_dir)
