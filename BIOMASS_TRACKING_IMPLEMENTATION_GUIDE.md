# Biomass tracking

Total colony biomass over time.

Biomass = (sum of active segment lengths) * cross_area * hy_density, i.e.
total hyphal volume times wall density. calculate_total_biomass(mycelia,
num_total_segs, params) returns it in grams. Active segments are those with
branch_id > -1.

Functions (helper_functions.py):
- calculate_total_biomass(mycelia, num_total_segs, params)
- plot_biomass_vs_time(times, biomass, ...) - single run
- plot_biomass_comparison(...) - several runs on one axis
- plot_errorbar_biomass(...) - mean +/- std across runs

Driver use: appended to count_biomass at each output step (the lengths it needs
are already on hand), stored in the output dict as 'total_biomass', and plotted
at the end. Multi-run averaging/comparison plots are produced in the averaging
block after all runs finish.

The standalone script plot_biomass_vs_time_using_csv_one_type.py reproduces the
single-run curve from a run folder's outputdata CSV.
