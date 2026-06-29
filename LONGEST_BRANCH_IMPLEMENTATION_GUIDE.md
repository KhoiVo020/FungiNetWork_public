# Longest branch tracking

Length of the longest single branch in the network over time.

calculate_longest_branch(mycelia, num_total_segs) returns that length (sums
segment lengths along each branch and takes the max).

Functions (helper_functions.py):
- calculate_longest_branch(mycelia, num_total_segs)
- plot_longest_branch_vs_time(times, longest, ...) - single run
- plot_longest_branch_comparison(...) - several runs
- plot_errorbar_longest_branch(...) - mean +/- std across runs

Driver use: appended to count_longest_branch at each output step and once after
the main loop, stored in the output dict as 'longest_branch_length', and
plotted at the end (single run) plus averaged/compared across runs in the
averaging block.

Tracked the same way as biomass, so the two curves can be compared directly.
