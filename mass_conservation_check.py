#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mass Conservation Checking Functions for Fungal Growth Model

Created on July 31, 2025
@author: Khoi

This module provides functions to check for mass conservation in the fungal 
growth simulation without modifying the structure of existing functions.
"""

import numpy as np
import warnings

class MassConservationTracker:
    """
    A class to track and check mass conservation throughout the simulation.
    
    This tracker monitors:
    - Total glucose in external and internal domains
    - Total cell wall materials
    - Total trehalose in external and internal domains
    - Mass balance during growth, uptake, metabolism, and transport processes
    """
    
    def __init__(self, params):
        """
        Initialize the mass conservation tracker.
        
        Parameters
        ----------
        params : dict
            Simulation parameters dictionary
        """
        self.params = params
        self.initial_totals = {}
        self.tolerance = {
            'relative': 1e-10,  # Relative tolerance for mass conservation
            'absolute': 1e-15   # Absolute tolerance for very small values
        }
        self.violations = []
        self.time_step = 0
        
        # Track cumulative processes
        self.cumulative_uptake = 0.0
        self.cumulative_metabolism = 0.0
        self.cumulative_growth_cost = 0.0
        self.cumulative_release = 0.0
        
    def initialize_tracking(self, sub_e_gluc, sub_e_treha, mycelia, num_total_segs):
        """
        Initialize mass tracking at the beginning of simulation.
        
        Parameters
        ----------
        sub_e_gluc : array
            External glucose grid
        sub_e_treha : array  
            External trehalose grid
        mycelia : dict
            Mycelia structure dictionary
        num_total_segs : int
            Number of active segments
        """
        self.initial_totals = self._calculate_total_masses(
            sub_e_gluc, sub_e_treha, mycelia, num_total_segs
        )
        self.time_step = 0
        print("Mass Conservation Tracking Initialized:")
        for key, value in self.initial_totals.items():
            print(f"  Initial {key}: {value:.6e}")
    
    def _calculate_total_masses(self, sub_e_gluc, sub_e_treha, mycelia, num_total_segs):
        """
        Calculate total masses in the system.
        
        Returns
        -------
        dict
            Dictionary containing total masses
        """
        totals = {}
        
        # External domain masses
        totals['glucose_external'] = np.sum(sub_e_gluc)
        totals['trehalose_external'] = np.sum(sub_e_treha)
        
        # Internal domain masses (only for active segments)
        active_mask = mycelia['branch_id'][:num_total_segs] >= 0
        totals['glucose_internal'] = np.sum(mycelia['gluc_i'][:num_total_segs][active_mask])
        totals['cellwall_internal'] = np.sum(mycelia['cw_i'][:num_total_segs][active_mask])
        totals['trehalose_internal'] = np.sum(mycelia['treha_i'][:num_total_segs][active_mask])
        
        # Total masses
        totals['glucose_total'] = totals['glucose_external'] + totals['glucose_internal']
        totals['trehalose_total'] = totals['trehalose_external'] + totals['trehalose_internal']
        totals['total_biomass'] = totals['cellwall_internal']
        
        # Total carbon (simplified; uses approximate carbon counts per molecule)
        totals['total_carbon'] = (totals['glucose_total'] * 6 +  # 6 carbons per glucose
                                 totals['trehalose_total'] * 12 + # 12 carbons per trehalose
                                 totals['cellwall_internal'] * 6)  # Approximate carbons in cell wall
        
        return totals
    
    def check_mass_conservation(self, sub_e_gluc, sub_e_treha, mycelia, num_total_segs,
                              process_name="General", detailed_output=False):
        """
        Check mass conservation at current time step.
        
        Parameters
        ----------
        sub_e_gluc : array
            Current external glucose grid
        sub_e_treha : array
            Current external trehalose grid  
        mycelia : dict
            Current mycelia structure
        num_total_segs : int
            Current number of active segments
        process_name : str
            Name of the process being checked
        detailed_output : bool
            Whether to print detailed output
            
        Returns
        -------
        bool
            True if mass is conserved within tolerance, False otherwise
        """
        self.time_step += 1
        current_totals = self._calculate_total_masses(
            sub_e_gluc, sub_e_treha, mycelia, num_total_segs
        )
        
        violations = []
        
        # Check conservation for non-metabolized components
        conserved_components = [
            'trehalose_total',  # Trehalose should be conserved (uptake = release)
        ]
        
        for component in conserved_components:
            initial = self.initial_totals[component]
            current = current_totals[component]
            
            # Calculate relative and absolute differences
            abs_diff = abs(current - initial)
            rel_diff = abs_diff / max(abs(initial), self.tolerance['absolute'])
            
            if (abs_diff > self.tolerance['absolute'] and 
                rel_diff > self.tolerance['relative']):
                violation = {
                    'component': component,
                    'initial': initial,
                    'current': current,
                    'absolute_diff': abs_diff,
                    'relative_diff': rel_diff,
                    'process': process_name,
                    'time_step': self.time_step
                }
                violations.append(violation)
        
        # Check glucose + biomass conservation (accounting for metabolism)
        # Glucose can be converted to biomass, so check: glucose_total + biomass_equivalent
        initial_glucose_biomass = (self.initial_totals['glucose_total'] + 
                                  self.initial_totals['total_biomass'] / self.params.get('yield_c', 1.0))
        current_glucose_biomass = (current_totals['glucose_total'] + 
                                  current_totals['total_biomass'] / self.params.get('yield_c', 1.0))
        
        abs_diff = abs(current_glucose_biomass - initial_glucose_biomass)
        rel_diff = abs_diff / max(abs(initial_glucose_biomass), self.tolerance['absolute'])
        
        if (abs_diff > self.tolerance['absolute'] and 
            rel_diff > self.tolerance['relative']):
            violation = {
                'component': 'glucose_plus_biomass',
                'initial': initial_glucose_biomass,
                'current': current_glucose_biomass,
                'absolute_diff': abs_diff,
                'relative_diff': rel_diff,
                'process': process_name,
                'time_step': self.time_step
            }
            violations.append(violation)
        
        # Store violations
        self.violations.extend(violations)
        
        # Print results if requested or if violations found
        if detailed_output or violations:
            self._print_conservation_report(current_totals, violations, process_name)
        
        return len(violations) == 0
    
    def check_process_balance(self, uptake_amount=0.0, metabolism_amount=0.0, 
                            growth_cost=0.0, release_amount=0.0, process_name=""):
        """
        Check mass balance for a specific process.
        
        Parameters
        ----------
        uptake_amount : float
            Amount of material taken up
        metabolism_amount : float
            Amount of material metabolized
        growth_cost : float
            Amount of material used for growth
        release_amount : float
            Amount of material released
        process_name : str
            Name of the process
            
        Returns
        -------
        bool
            True if process is balanced
        """
        # Update cumulative totals
        self.cumulative_uptake += uptake_amount
        self.cumulative_metabolism += metabolism_amount
        self.cumulative_growth_cost += growth_cost
        self.cumulative_release += release_amount
        
        # For uptake/release processes, input should equal output
        if process_name.lower() in ['uptake', 'release']:
            net_change = uptake_amount - release_amount
            if abs(net_change) > self.tolerance['absolute']:
                print(f"Warning: {process_name} not balanced. Net change: {net_change}")
                return False
        
        return True
    
    def _print_conservation_report(self, current_totals, violations, process_name):
        """Print detailed conservation report."""
        print(f"\n=== Mass Conservation Check: {process_name} (Step {self.time_step}) ===")
        
        for key, current in current_totals.items():
            if key in self.initial_totals:
                initial = self.initial_totals[key]
                change = current - initial
                if abs(initial) > self.tolerance['absolute']:
                    rel_change = change / initial * 100
                    print(f"{key:20s}: {current:12.6e} (Δ={change:+12.6e}, {rel_change:+6.2f}%)")
                else:
                    print(f"{key:20s}: {current:12.6e} (Δ={change:+12.6e})")
        
        if violations:
            print("\n*** MASS CONSERVATION VIOLATIONS ***")
            for v in violations:
                print(f"  {v['component']}: rel_diff = {v['relative_diff']:.2e} "
                      f"(tolerance = {self.tolerance['relative']:.2e})")
        else:
            print("Mass conservation satisfied within tolerance")
    
    def get_conservation_summary(self):
        """
        Get summary of all conservation violations.
        
        Returns
        -------
        dict
            Summary of violations and statistics
        """
        summary = {
            'total_violations': len(self.violations),
            'time_steps_checked': self.time_step,
            'violations_by_component': {},
            'violations_by_process': {},
            'cumulative_processes': {
                'uptake': self.cumulative_uptake,
                'metabolism': self.cumulative_metabolism,
                'growth_cost': self.cumulative_growth_cost,
                'release': self.cumulative_release
            }
        }
        
        # Group violations by component and process
        for v in self.violations:
            comp = v['component']
            proc = v['process']
            
            if comp not in summary['violations_by_component']:
                summary['violations_by_component'][comp] = 0
            summary['violations_by_component'][comp] += 1
            
            if proc not in summary['violations_by_process']:
                summary['violations_by_process'][proc] = 0
            summary['violations_by_process'][proc] += 1
        
        return summary
    
    def print_final_report(self):
        """Print final conservation report."""
        summary = self.get_conservation_summary()
        
        print("\n" + "="*60)
        print("FINAL MASS CONSERVATION REPORT")
        print("="*60)
        print(f"Time steps checked: {summary['time_steps_checked']}")
        print(f"Total violations: {summary['total_violations']}")
        
        if summary['total_violations'] == 0:
            print("All mass conservation checks passed")
        else:
            print("Mass conservation violations detected:")
            
            print("\nViolations by component:")
            for comp, count in summary['violations_by_component'].items():
                print(f"  {comp}: {count}")
            
            print("\nViolations by process:")
            for proc, count in summary['violations_by_process'].items():
                print(f"  {proc}: {count}")
        
        print(f"\nCumulative process totals:")
        for proc, total in summary['cumulative_processes'].items():
            print(f"  {proc}: {total:.6e}")


def check_diffusion_conservation(sub_e_before, sub_e_after, process_name="Diffusion"):
    """
    Check mass conservation for diffusion processes.
    
    Parameters
    ----------
    sub_e_before : array
        Substrate grid before diffusion
    sub_e_after : array
        Substrate grid after diffusion
    process_name : str
        Name of diffusion process
        
    Returns
    -------
    bool
        True if mass is conserved
    """
    total_before = np.sum(sub_e_before)
    total_after = np.sum(sub_e_after)
    diff = abs(total_after - total_before)
    
    tolerance = 1e-12
    relative_diff = diff / max(abs(total_before), 1e-15)
    
    if relative_diff > tolerance:
        print(f"Warning: {process_name} - Mass not conserved!")
        print(f"  Before: {total_before:.6e}")
        print(f"  After:  {total_after:.6e}")
        print(f"  Diff:   {diff:.6e} (relative: {relative_diff:.2e})")
        return False
    
    return True


def check_translocation_conservation(mycelia_before, mycelia_after, num_total_segs, 
                                   component='gluc_i', process_name="Translocation"):
    """
    Check mass conservation for translocation processes.
    
    Parameters
    ----------
    mycelia_before : dict
        Mycelia state before translocation
    mycelia_after : dict
        Mycelia state after translocation
    num_total_segs : int
        Number of active segments
    component : str
        Component to check ('gluc_i', 'cw_i', 'treha_i')
    process_name : str
        Name of process
        
    Returns
    -------
    bool
        True if mass is conserved
    """
    active_mask = mycelia_after['branch_id'][:num_total_segs] >= 0
    
    total_before = np.sum(mycelia_before[component][:num_total_segs][active_mask])
    total_after = np.sum(mycelia_after[component][:num_total_segs][active_mask])
    diff = abs(total_after - total_before)
    
    tolerance = 1e-12
    relative_diff = diff / max(abs(total_before), 1e-15)
    
    if relative_diff > tolerance:
        print(f"Warning: {process_name} ({component}) - Mass not conserved!")
        print(f"  Before: {total_before:.6e}")
        print(f"  After:  {total_after:.6e}")
        print(f"  Diff:   {diff:.6e} (relative: {relative_diff:.2e})")
        return False
    
    return True


def check_growth_mass_balance(growth_costs, biomass_created, conversion_factor=1.0):
    """
    Check mass balance during growth processes.
    
    Parameters
    ----------
    growth_costs : array
        Costs of growth in terms of consumed materials
    biomass_created : array
        Amount of biomass created
    conversion_factor : float
        Conversion factor between consumed material and biomass
        
    Returns
    -------
    bool
        True if growth is mass-balanced
    """
    total_consumed = np.sum(growth_costs) * conversion_factor
    total_created = np.sum(biomass_created)
    diff = abs(total_created - total_consumed)
    
    tolerance = 1e-10
    relative_diff = diff / max(abs(total_consumed), 1e-15)
    
    if relative_diff > tolerance:
        print(f"Warning: Growth - Mass balance violated!")
        print(f"  Material consumed: {total_consumed:.6e}")
        print(f"  Biomass created:   {total_created:.6e}")
        print(f"  Diff:             {diff:.6e} (relative: {relative_diff:.2e})")
        return False
    
    return True


def verify_segment_consistency(mycelia, num_total_segs):
    """
    Verify consistency of segment properties.
    
    Parameters
    ----------
    mycelia : dict
        Mycelia structure
    num_total_segs : int
        Number of active segments
        
    Returns
    -------
    bool
        True if segments are consistent
    """
    issues = []
    
    # Check for negative masses
    active_mask = mycelia['branch_id'][:num_total_segs] >= 0
    
    for component in ['gluc_i', 'cw_i', 'treha_i']:
        if component in mycelia:
            values = mycelia[component][:num_total_segs][active_mask]
            if np.any(values < 0):
                neg_count = np.sum(values < 0)
                min_val = np.min(values)
                issues.append(f"Negative {component}: {neg_count} segments, min={min_val:.2e}")
    
    # Check for unrealistic volumes
    if 'seg_vol' in mycelia:
        volumes = mycelia['seg_vol'][:num_total_segs][active_mask]
        if np.any(volumes <= 0):
            zero_count = np.sum(volumes <= 0)
            issues.append(f"Non-positive volumes: {zero_count} segments")
    
    # Check for NaN or infinite values
    for key in ['gluc_i', 'cw_i', 'treha_i', 'seg_length', 'seg_vol']:
        if key in mycelia:
            values = mycelia[key][:num_total_segs]
            if np.any(~np.isfinite(values)):
                issues.append(f"Non-finite values in {key}")
    
    if issues:
        print("Segment consistency issues found:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    
    return True
