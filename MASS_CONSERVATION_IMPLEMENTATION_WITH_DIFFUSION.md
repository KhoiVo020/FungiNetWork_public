# Mass conservation checks

Two independent sanity checks: a biological mass balance and a numerical
diffusion check. Both are diagnostics only; they don't change the dynamics.

## Biological mass balance (glucose equivalents)

Everything is converted back to glucose equivalents and compared to the initial
total:

  Initial_total_glucose = free_glucose
                        + cell_wall / yield_c
                        + trehalose / convert_metabolite
                        + unmodeled_metabolites

cell wall and trehalose are divided by their yield fractions to recover the
glucose they came from; the remaining consumed glucose is the unmodeled part.

Functions (helper_functions.py):
- calculate_total_mass_with_metabolites(mycelia, num_total_segs, sub_e_gluc,
  sub_e_treha, params)
- check_mass_conservation(initial_mass_dict, current_mass_dict)
- print_initial_mass_state(...), print_mass_conservation_report(...)

Driver use: initial_mass_dict is computed once after setup; at each output step
the current mass is compared to it and a short report is printed.

## Diffusion (numerical) check

Confirms the glucose diffusion step neither creates nor destroys mass on the
active grid. check_diffusion_conservation() compares the masked grid sum before
and after diffusion; print_diffusion_conservation_warning() prints if the
relative error exceeds tolerance.

This is gated by DIFFUSION_CONSERVATION_CHECK_INTERVAL in the driver:
- 0  -> off (no grid copy, no check) - default
- N>0 -> check every Nth diffusion sub-step

The driver tracks how many sub-steps were checked, how many failed, and the max
relative error, and prints a one-line summary at each output step plus a final
summary at the end of the run.
