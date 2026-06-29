# Restarting a simulation from a checkpoint

How to resume a run from a saved checkpoint instead of starting over.

## How restart works

The driver reads params['restart'] from parameters.ini:
- 0 (default): fresh start. current_time = 0, and the log is opened in 'w' mode
  so it overwrites simulation_log_run0.txt. Does not resume.
- 1: loads Results/{output_path}/Run{run}/restart.pkl (name built in code).
- 2: loads the exact file named in params['restart_file']. Use this to resume
  from a timestamped restart_t=*.pkl checkpoint.

During a run the driver periodically writes restart_t={current_time}.pkl into
the run's Results/.../Run{run}/ folder. Each one stores mycelia, num_total_segs,
dtt, sub_e_gluc, sub_e_treha, current_time.

Use restart = 2 to resume from the newest timestamped checkpoint without copying
it to restart.pkl.

## Steps

Commands below are PowerShell from the project root.

1. Find the checkpoint. List them newest-last:
```powershell
Get-ChildItem "Results\Fusion_production_PatchyEnv3Grid20\Run0\restart_t=*.pkl" |
  Sort-Object { [double]($_.Name -replace 'restart_t=|\.pkl','') } |
  Select-Object -Last 5 Name
```
Pick the latest to lose the least progress (or an earlier one to roll back).

2. Back up the log if you want it - rerunning overwrites simulation_log_run0.txt:
```powershell
Copy-Item "Results\...\Run0\simulation_log_run0.txt" "Results\...\Run0\simulation_log_run0_pre-restart.txt"
```

3. Put the full checkpoint path into parameters.ini (around line 56), e.g.:
```ini
restart_file = 'D:\...\Run0\restart_t=928000.0.pkl'
```

4. Set restart = 2 in parameters.ini. Leave output_path unchanged so new output
   continues in the same results folder.

5. Check there's work left: final_time must be greater than the checkpoint time,
   otherwise the loop exits immediately.

6. Run:
```powershell
python driver_fungalGrowth_singleNutrient.py -i parameters.ini
```
It should print "Current time: 928000.0 ..." and continue. New output is written
for times after the resume point; earlier files are left alone.

7. When done, set restart = 0 again so the next fresh run doesn't try to resume.

## Notes / cautions

- The file in restart_file must exist before running. If the open fails the
  driver falls through to pickle.load and crashes with a NameError, so don't
  skip step 3.
- Quoted paths are fine - get_configs strips surrounding quotes from
  restart_file.
- Old checkpoints don't have the grid_cell_map field; it's rebuilt
  automatically on load by growth_functions._ensure_grid_cell_map() from the
  saved xy_e_idx, so old and new checkpoints both resume fine.
