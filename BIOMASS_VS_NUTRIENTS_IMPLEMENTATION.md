# Biomass vs nutrients consumed

Plots biomass against cumulative glucose consumed instead of against time, to
see growth yield (biomass produced per unit nutrient).

Nutrients consumed at a step = initial_total_nutrients - current sum of
sub_e_gluc; accumulated in count_nutrients_consumed alongside count_biomass.

Functions (helper_functions.py):
- plot_biomass_vs_nutrients_consumed(biomass, nutrients, ...) - single run
- plot_biomass_vs_nutrients_comparison(...) - several runs
- plot_errorbar_biomass_vs_nutrients(...) - mean +/- std across runs

Driver use: single-run plot is made after plot_biomass_vs_time; the averaged
and comparison versions are made in the multi-run averaging block.

A roughly straight line means a constant conversion efficiency; a flattening
curve means biomass gain per nutrient is dropping.
