# Saving console output to a log file

Each run mirrors everything printed to the console into its own log file, so a
run's full output survives after the terminal closes.

How it works: a small DualLogger class (defined inside the driver run function)
wraps stdout and writes every line to both the console and a file. stdout/stderr
are pointed at it at the start of the run and restored in a finally block, so
the log is closed cleanly even if the run errors out.

Log file: simulation_log_run{run}.txt under the run's log directory. With
num_parallel_runs > 1 each run gets its own file, so parallel runs don't
interleave.

The file is line-buffered and flushed on each write, so it stays current and
isn't lost if the job is killed. It captures everything: setup, per-step
progress, the mass-conservation reports, warnings, tracebacks, and the final
runtime summary.

If you'd rather not touch the code, you can also just redirect at the shell:
`python driver_fungalGrowth_singleNutrient.py -i parameters.ini > run.log 2>&1`.
