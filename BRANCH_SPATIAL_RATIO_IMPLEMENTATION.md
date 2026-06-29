# Branch spatial ratio

Where each branch sits relative to the colony's radial extent, at a given time.

Definition (origin defaults to (0,0)):
- R_max: Euclidean distance from origin to the farthest tip (the expansion
  limit at time t).
- For each branch b, r_b is the distance from origin to the branch base
  (branching point), and the spatial ratio is rho_b = r_b / R_max.
- rho_b is in [0,1]: 0 = branch starts at the inoculum, ~1 = branch out near
  the periphery.

Functions (helper_functions.py):
- calculate_branch_spatial_ratios(mycelia, num_total_segs, origin=None)
  returns a dict with origin, R_max, farthest tip index/xy, num_branches, and
  branch_data: per branch its branch_id, base_xy, r_b, rho_b, base segment
  index and parent segment index.
- save_branch_spatial_ratios_csv(result, folder_string, param_string,
  current_time, run, params) writes one CSV per call.

Driver use: called at each output step and at the final time; CSVs go to
Results/{param_string}/Run{run}/.
