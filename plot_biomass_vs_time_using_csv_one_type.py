#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 10:18:44 2026

@author: Khoi

Reproduce the biomass-vs-time plot from a driver run folder.

Given a folder containing `*_outputdata_*_t=*.csv` snapshots, the script
finds the snapshot with the largest simulation time (which contains the
full cumulative history) and plots `total_biomass` vs `array_times`.

If no folder is given, the script walks `Results/` and produces one plot
per discovered run folder, so all biomass curves end up side by side in
`Biomass_vs_time_comparison/`.

Usage as a script:
    # Auto-discover every run folder under Results/
    python plot_biomass_vs_time_using_csv.py

    # Or point at a specific run folder
    python plot_biomass_vs_time_using_csv.py <run-folder>
                                              [--units days|hours|minutes|seconds]
                                              [--out-dir Biomass_vs_time_comparison]

Usage as a module:
    from plot_biomass_vs_time_using_csv import plot_BIOMASS_VS_TIME_USING_CSV
    plot_BIOMASS_VS_TIME_USING_CSV("Results/Fusion_production_PatchyEnv3Grid20/Run0")
"""

import os
import re
import csv
import ast
import glob
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')
sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})


OUTPUT_CSV_PATTERN = "*_outputdata_*_t=*.csv"
DEFAULT_OUTPUT_DIR = "Biomass_vs_time_comparison"
DEFAULT_SEARCH_ROOT = "Results"
_T_RE = re.compile(r"_t=([0-9]+(?:\.[0-9]+)?)\.csv$", re.IGNORECASE)


def discover_run_folders(search_root=DEFAULT_SEARCH_ROOT):
    """
    Walk `search_root` and return every directory that directly contains at
    least one `*_outputdata_*_t=*.csv` file. Sorted for stable output.
    """
    if not os.path.isdir(search_root):
        return []
    found = set()
    for dirpath, _, filenames in os.walk(search_root):
        for name in filenames:
            if _T_RE.search(name) and "_outputdata_" in name:
                found.add(dirpath)
                break
    return sorted(found)


def find_longest_time_csv(csv_folder):
    """
    Scan `csv_folder` for `*_outputdata_*_t=*.csv` files and return the path
    of the one with the largest simulation time encoded in its filename.

    Raises FileNotFoundError if no matching CSV is found.
    """
    candidates = glob.glob(os.path.join(csv_folder, OUTPUT_CSV_PATTERN))
    if not candidates:
        raise FileNotFoundError(
            "No '{}' files found in folder: {}".format(OUTPUT_CSV_PATTERN, csv_folder))

    best_path = None
    best_t = -float("inf")
    for path in candidates:
        match = _T_RE.search(os.path.basename(path))
        if not match:
            continue
        t = float(match.group(1))
        if t > best_t:
            best_t = t
            best_path = path

    if best_path is None:
        raise FileNotFoundError(
            "Found CSV files in {} but none matched the expected '_t=<number>.csv' pattern.".format(csv_folder))

    return best_path, best_t


def _read_biomass_series(csv_file_path):
    """Return (times, biomass) numpy arrays from a driver outputdata CSV."""
    times = None
    biomass = None
    with open(csv_file_path, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            key = row[0].strip()
            if key == 'array_times' and len(row) >= 2:
                times = np.array(ast.literal_eval(row[1]), dtype=float)
            elif key == 'total_biomass' and len(row) >= 2:
                biomass = np.array(ast.literal_eval(row[1]), dtype=float)
            if times is not None and biomass is not None:
                break

    if times is None or biomass is None:
        raise ValueError(
            "Could not find both 'array_times' and 'total_biomass' rows in {}".format(csv_file_path))
    return times, biomass


def plot_BIOMASS_VS_TIME_USING_CSV(csv_folder, params=None, output_dir=None):
    """
    Find the longest-time outputdata CSV in `csv_folder`, then plot
    `total_biomass` vs `array_times` from it.

    Parameters
    ----------
    csv_folder : str
        Folder containing one or more `*_outputdata_*_t=*.csv` snapshots.
    params : dict, optional
        Only `plot_units_time` is consulted ('days', 'hours', 'minutes', else
        seconds). Defaults to days.
    output_dir : str, optional
        Folder where the PNG is written. Created if missing.
        Defaults to `Biomass_vs_time_comparison/` in the current working directory.

    Returns
    -------
    str
        Absolute path to the saved PNG.
    """
    csv_file_path, max_t = find_longest_time_csv(csv_folder)
    print("Using CSV with largest time (t={}): {}".format(max_t, csv_file_path))

    times, biomass = _read_biomass_series(csv_file_path)

    plot_units_time = (params or {}).get('plot_units_time', 'days')
    if plot_units_time == 'days':
        plot_times = times / (60 * 60 * 24)
        time_label = 'Time (days)'
    elif plot_units_time == 'hours':
        plot_times = times / (60 * 60)
        time_label = 'Time (hours)'
    elif plot_units_time == 'minutes':
        plot_times = times / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = times
        time_label = 'Time (seconds)'

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.plot(plot_times, biomass, linewidth=2, marker='o', markersize=4, color='green')
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel('Total Biomass (g)', fontsize=12, fontweight='bold')
    ax.set_title('Total Fungal Biomass vs Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    final_biomass = biomass[-1]
    ax.text(0.95, 0.05, f'Final: {final_biomass:.6f} g',
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    sns.despine()
    plt.tight_layout()

    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    csv_stem = os.path.splitext(os.path.basename(csv_file_path))[0]
    save_path = os.path.abspath(os.path.join(output_dir, csv_stem + "_biomass_vs_time.png"))
    fig.savefig(save_path)
    plt.close(fig)

    print(f"Biomass plot saved to: {save_path}")
    print(f"Final biomass: {final_biomass:.6f} g")
    return save_path


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Plot total biomass vs time from the longest-time outputdata CSV in a run folder. "
                    "If no folder is given, walks Results/ and plots every run found.")
    parser.add_argument("csv_folder", nargs="?", default=None,
                        help="Folder containing *_outputdata_*_t=*.csv files (e.g. Results/.../Run0). "
                             "If omitted, every run folder under Results/ is processed.")
    parser.add_argument("--search-root", default=DEFAULT_SEARCH_ROOT,
                        help="Root to scan when csv_folder is omitted (default: %(default)s)")
    parser.add_argument("--units", default="days",
                        choices=["days", "hours", "minutes", "seconds"],
                        help="X-axis time units (default: days)")
    parser.add_argument("--out-dir", default=DEFAULT_OUTPUT_DIR,
                        help="Output folder for the PNG(s) (default: %(default)s)")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    params = {'plot_units_time': args.units}

    if args.csv_folder is not None:
        folders = [args.csv_folder]
    else:
        folders = discover_run_folders(args.search_root)
        if not folders:
            raise SystemExit(
                "No run folders containing '{}' were found under '{}'. "
                "Pass a folder explicitly, or use --search-root.".format(
                    OUTPUT_CSV_PATTERN, args.search_root))
        print("Auto-discovered {} run folder(s) under '{}':".format(
            len(folders), args.search_root))
        for f in folders:
            print("  - {}".format(f))

    for folder in folders:
        plot_BIOMASS_VS_TIME_USING_CSV(
            csv_folder=folder,
            params=params,
            output_dir=args.out_dir,
        )
