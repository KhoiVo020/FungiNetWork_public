# Longest branch vs nutrients consumed

Same idea as biomass-vs-nutrients, but for the longest branch: plot longest
branch length against cumulative glucose consumed instead of against time.

Functions (helper_functions.py):
- plot_longest_branch_vs_nutrients_consumed(longest, nutrients, ...) - single run
- plot_longest_branch_vs_nutrients_comparison(...) - several runs
- plot_errorbar_longest_branch_vs_nutrients(...) - mean +/- std across runs

Driver use: uses the same count_longest_branch and count_nutrients_consumed
arrays already collected for the time-series and biomass plots. Single-run plot
at the end of a run; averaged/comparison versions in the multi-run averaging
block.

Reading it: how much extension you get per unit nutrient. Compare against the
biomass-vs-nutrients curve to separate "growing longer" from "growing denser".
