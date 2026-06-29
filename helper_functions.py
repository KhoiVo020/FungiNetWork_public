#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 17:17:45 2020

@author: jolenebritton
"""


import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
#from matplotlib import collections, colors, transforms, 
import seaborn as sns
import sys
from matplotlib import pyplot as plt
# import pylab as pl
from matplotlib import collections  as mc
sns.set_style('white')
# sns.set_context("talk")
sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})
import configparser
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

np.set_printoptions(threshold=sys.maxsize)

# ----------------------------------------------------------------------------
# SET UP FUNCTIONS
# ----------------------------------------------------------------------------

def get_configs(config_filename):
    """
    Parameters
    ----------
    config_filename : str
        Filepath for the .ini config file containing parameters.

    Returns
    -------
    params_dict : dict
        A dictionary containing all parameters in a usable form.
    config : TYPE
        The config file that was loaded into this function.

    """
    config = configparser.ConfigParser()
    config.read(config_filename)
    # breakpoint()
    # Extract the sections
    discrete_params = config['DISCRETE PARAMS']
    growth_params = config['GROWTH PARAMS']
    nutrient_params = config['NUTRIENT PARAMS']

    # Generate extra values
    diam =  discrete_params.getfloat('hy_diam')
    sl =  discrete_params.getfloat('hy_compartment')
    cross_area = np.pi*(0.5*diam)**2
    init_vol_seg = sl*cross_area

    dy = discrete_params.getfloat('grid_len')
    # Best to have dy = 0.5*sl:
    #dy = 0.5*sl
    dz = discrete_params.getfloat('grid_height')
    vol_grid = dy*dy*dz
    diff_e_gluc = nutrient_params.getfloat('diffusion_e_gluc')
    convert_metabolite = nutrient_params.getfloat('convert_metabolite')
    diff_i_gluc = nutrient_params.getfloat('diffusion_i_gluc')
    kg1_wall = nutrient_params.getfloat('kg1_wall')
    hy_density =  growth_params.getfloat('hy_density')
    f_dw = growth_params.getfloat('f_dw')
    f_wall =  growth_params.getfloat('f_wall')
    f_cw_cellwall =  growth_params.getfloat('f_cw_cellwall')
    mw_cw = nutrient_params.getfloat('mw_cw')

    # dt = 0.99*(dy**2)/(4*diff_e_gluc)
    #dt = 0.75*min((sl**2/diff_i_gluc),(sl/kg1_wall)) # kg1_wall is radial growth rate and advection rate
    
    # The maximum rate of active transport is set to be the product of the maximum rate of radial growth, 
    # the cross sectional area of the hyphae, the dry weight of the cell, the fraction of the dry weight that is cell wall material,
    # and divided by the formula weight of the cell wall material,
    active_trsprt_vel_cw = kg1_wall*cross_area*hy_density*1.0e+12*f_dw*f_wall*f_cw_cellwall \
        /mw_cw
    advection_constant_gluc = kg1_wall*init_vol_seg # kg1_wall is the radial rate. 
                                                    # scaling it by the volume of the hyphae gives the rate that glucose
                                                    # would spread in 2-dimensions per second.

    #dt = 0.0025*min((sl**2/diff_i_gluc),(sl/kg1_wall)) #kg1_wall should be advection rate
    dt_i = 0.01*min((sl**2/diff_i_gluc),1/(active_trsprt_vel_cw*0.02))
    dt_e = 0.01*(dy**2/diff_e_gluc)
    #dt = min(dt, dy**2/(diff_e_gluc))
    #dt = 0.01*dt

    #dt = 22.5
    #dt = 0.5*0.75*min((sl**2/diff_i_gluc),(sl/linear_growth_rate))

    up_state = nutrient_params['up_state']
    if up_state == 'repressed':
        #Ku2 = nutrient_params.getfloat('Ku2a_gluc')*init_vol_seg
        Ku2 = nutrient_params.getfloat('Ku2a_gluc')*vol_grid

    else:
        #Ku2 = nutrient_params.getfloat('Ku2b_gluc')*init_vol_seg
        Ku2 = nutrient_params.getfloat('Ku2b_gluc')*vol_grid

    # Save to a dictionary
    params_dict = {
        # SECTION 1: Discretization Parameters
        'dt_i' : dt_i,
        'dt_e' : dt_e,
        'final_time' : discrete_params.getfloat('final_time'),
        'plot_units_time' : discrete_params['plot_units_time'],

        'sl' : sl,
        'dy' : dy,
        'vol_grid': vol_grid,
        'plot_units_space' : discrete_params['plot_units_space'],
        'init_segs_count' : discrete_params.getint('init_segs_count'),
        'environ_type' : discrete_params['environ_type'],
        'cross_area' : cross_area,
        'init_vol_seg' : init_vol_seg,
        'septa_len' : 1,
        'grid_scale_val' : discrete_params.getfloat('grid_scale_val'),
        'hy_diam' : discrete_params.getfloat('hy_diam'),
        'fungal_fusion' : discrete_params.getint('fungal_fusion'),
        'isTipRelease' : discrete_params.getint('isTipRelease'),
        'restrictBranching' : discrete_params.getint('restrictBranching'),
        'chance_to_fuse' : discrete_params.getfloat('chance_to_fuse'),
        'is_pathchyEnv' : discrete_params.getint('is_pathchyEnv'),
        'diffusion_on' : discrete_params.getint('diffusion_on'),
        'num_parallel_runs' : discrete_params.getint('num_parallel_runs'),
        'output_path' : discrete_params['output_path'],
        'restart' : discrete_params.getint('restart'),
        'restart_file' : discrete_params.get('restart_file', '').strip().strip("\"'"),
        'setPatchyEnv' : discrete_params.getint('setPatchyEnv'),



        # SECTION 2: Extension & Branching for Growth Parameters
        'angle_sd' : growth_params.getfloat('angle_sd')*(np.pi/180),
        'branch_mean' : growth_params.getfloat('branch_mean')*(np.pi/180),
        'branch_sd' : growth_params.getfloat('branch_sd')*(np.pi/180),
        'branch_cost' : growth_params.getfloat('branch_cost'),
        'branch_rate' : growth_params.getfloat('branch_rate'),
        'hy_density' :  growth_params.getfloat('hy_density'),
        'f_dw' :  growth_params.getfloat('f_dw'),
        'f_wall' :  growth_params.getfloat('f_wall'),
        'f_cw_cellwall' :  growth_params.getfloat('f_cw_cellwall'),

        # SECTION 3: Internal & External Nutrient Parameters
        'init_sub_e_dist' : nutrient_params['init_sub_e_dist'],
        'init_sub_e_gluc' : nutrient_params.getfloat('init_sub_e_gluc')*vol_grid,
        'init_sub_e_treha' : nutrient_params.getfloat('init_sub_e_treha')*vol_grid,
        'var_nutrient_backgrnd' : nutrient_params.getfloat('var_nutrient_backgrnd'),
        'diffusion_e_gluc' : diff_e_gluc,
        'convert_metabolite': convert_metabolite,

        'init_sub_i_gluc' : nutrient_params.getfloat('init_sub_i_gluc'),
        'diffusion_i_gluc' : nutrient_params.getfloat('diffusion_i_gluc'),
        'vel_gluc' : nutrient_params.getfloat('vel_gluc'),
        
        # 'm_gluc' : nutrient_params.getfloat('m_gluc'),
        # 'rho' : nutrient_params.getfloat('rho')*vol_seg,

        'ku1_gluc' : nutrient_params.getfloat('ku1_gluc'),
        'Ku2_gluc' : Ku2,
        'yield_u' : nutrient_params.getfloat('yield_u'),

        #'kc1_gluc' : nutrient_params.getfloat('kc1_gluc'),
        #'Kc2_gluc' : nutrient_params.getfloat('Kc2_gluc'),
        'kc1_gluc' : nutrient_params.getfloat('kc1_gluc'),
        'Kc2_gluc' : Ku2*init_vol_seg/vol_grid,
        'yield_c' : nutrient_params.getfloat('yield_c'),

        'kg1_wall' : nutrient_params.getfloat('kg1_wall'),
        # 'Kg2_wall' : nutrient_params.getfloat('Kg2_wall')*vol_seg,
        'Kg2_wall' : nutrient_params.getfloat('Kg2_wall'),
        'mw_cw' : nutrient_params.getfloat('mw_cw'),
        'mw_glucose' : nutrient_params.getfloat('mw_glucose'),
        'active_trsprt_vel_cw' : active_trsprt_vel_cw,

        'num_v' : nutrient_params.getfloat('num_v')
    }
    # The rate of glucose uptake is determined from the amount of glucose needed to support the growth rate.
    # The rate of uptake of glucose (kc1_gluc) is the product of the rate of radial growth (kg1_wall)
    # times the cross sectinoal area of the hyphae, the hyphal density, the fraction of wet cell mass that is dry cell mass
    # the fraction of dry cell mass that is cell wall material, the fraction of cell wall material that is composed of sugars (chitin and glucan),
    # the fraction of glucose mass that is provided by metabolism for cell wall raw material (yield_c),
    # and all divided by the formula weight of the chitin and glucan cell wall material.
    # The factor of 1.0e+03 is to convert the rate from moles to millimoles.
    kc1_gluc = params_dict['kg1_wall']*params_dict['cross_area']*params_dict['hy_density']\
                *params_dict['f_dw']*params_dict['f_wall']*params_dict['f_cw_cellwall'] \
                /(params_dict['mw_cw']*params_dict['yield_c'])*1.0e+03
    params_dict['kc1_gluc'] = kc1_gluc
    params_dict['ku1_gluc'] = kc1_gluc

    #if not('yield_c_in_mmoles' in params_dict):
    #	params_dict['yield_c_in_mmoles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
    params_dict['yield_c_in_mmoles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
    
    use_original = 0

    if(use_original != 1):
        # Max rate of moles of cell wall raw material used per time step:
        max_gms_cw_per_time = params_dict['kg1_wall']* np.pi*(params_dict['hy_diam']/2.0)**2.0 \
                        * params_dict['hy_density'] * params_dict['f_dw']* params_dict['f_wall'] * params_dict['f_cw_cellwall']
        max_moles_cw_per_time = max_gms_cw_per_time / params_dict['mw_cw']

        #Max rate of conversion of glucose to cell wall material is the maxrate of
        # grams of cell wall used per step/timestep / (gms cw produced per gms glucose used)/ (MW glucose)
        params_dict['yield_in_moles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
        #params_dict['kc1_gluc'] =  max_moles_cw_per_time/yield_in_moles
        #params_dict['ku1_gluc'] = params_dict['kc1_gluc']


    return params_dict, config

def get_configs_new(config_filename):
    """
    Parameters
    ----------
    config_filename : str
        Filepath for the .ini config file containing parameters.

    Returns
    -------
    params_dict : dict
        A dictionary containing all parameters in a usable form.
    config : TYPE
        The config file that was loaded into this function.

    """
    config = configparser.ConfigParser()
    config.read(config_filename)
    # breakpoint()
    # Extract the sections
    simulation_params = config['SIMULATION']
    spatial_params = config['SPATIAL']
    environment_params = config['ENVIRONMENT']
    mycelia_params = config['MYCELIA']

    # Generate extra values
    diam =  discrete_params.getfloat('hy_diam')
    sl =  discrete_params.getfloat('hy_compartment')
    cross_area = np.pi*(0.5*diam)**2
    init_vol_seg = sl*cross_area

    dy = discrete_params.getfloat('grid_len')
    # Best to have dy = 0.5*sl:
    #dy = 0.5*sl
    dz = discrete_params.getfloat('grid_height')
    vol_grid = dy*dy*dz
    diff_e_gluc = nutrient_params.getfloat('diffusion_e_gluc')
    convert_metabolite = nutrient_params.getfloat('convert_metabolite')
    diff_i_gluc = nutrient_params.getfloat('diffusion_i_gluc')
    kg1_wall = nutrient_params.getfloat('kg1_wall')
    hy_density =  growth_params.getfloat('hy_density')
    f_dw = growth_params.getfloat('f_dw')
    f_wall =  growth_params.getfloat('f_wall')
    f_cw_cellwall =  growth_params.getfloat('f_cw_cellwall')
    mw_cw = nutrient_params.getfloat('mw_cw')

    # dt = 0.99*(dy**2)/(4*diff_e_gluc)
    #dt = 0.75*min((sl**2/diff_i_gluc),(sl/kg1_wall)) # kg1_wall is radial growth rate and advection rate
    
    # The maximum rate of active transport is set to be the product of the maximum rate of radial growth, 
    # the cross sectional area of the hyphae, the dry weight of the cell, the fraction of the dry weight that is cell wall material,
    # and divided by the formula weight of the cell wall material,
    active_trsprt_vel_cw = kg1_wall*cross_area*hy_density*1.0e+12*f_dw*f_wall*f_cw_cellwall \
        /mw_cw
    advection_constant_gluc = kg1_wall*init_vol_seg # kg1_wall is the radial rate. 
                                                    # scaling it by the volume of the hyphae gives the rate that glucose
                                                    # would spread in 2-dimensions per second.

    #dt = 0.0025*min((sl**2/diff_i_gluc),(sl/kg1_wall)) #kg1_wall should be advection rate
    dt_i = 0.01*min((sl**2/diff_i_gluc),1/(active_trsprt_vel_cw*0.02))
    dt_e = 0.01*(dy**2/diff_e_gluc)
    #dt = min(dt, dy**2/(diff_e_gluc))
    #dt = 0.01*dt

    #dt = 22.5
    #dt = 0.5*0.75*min((sl**2/diff_i_gluc),(sl/linear_growth_rate))

    up_state = nutrient_params['up_state']
    if up_state == 'repressed':
        #Ku2 = nutrient_params.getfloat('Ku2a_gluc')*init_vol_seg
        Ku2 = nutrient_params.getfloat('Ku2a_gluc')*vol_grid

    else:
        #Ku2 = nutrient_params.getfloat('Ku2b_gluc')*init_vol_seg
        Ku2 = nutrient_params.getfloat('Ku2b_gluc')*vol_grid

    # Save to a dictionary
    params_dict = {
        # SECTION 1: Discretization Parameters
        'dt_i' : dt_i,
        'dt_e' : dt_e,
        'final_time' : discrete_params.getfloat('final_time'),
        'plot_units_time' : discrete_params['plot_units_time'],

        'sl' : sl,
        'dy' : dy,
        'vol_grid': vol_grid,
        'plot_units_space' : discrete_params['plot_units_space'],
        'init_segs_count' : discrete_params.getint('init_segs_count'),
        'environ_type' : discrete_params['environ_type'],
        'cross_area' : cross_area,
        'init_vol_seg' : init_vol_seg,
        'septa_len' : 1,
        'grid_scale_val' : discrete_params.getfloat('grid_scale_val'),
        'hy_diam' : discrete_params.getfloat('hy_diam'),

        # SECTION 2: Extension & Branching for Growth Parameters
        'angle_sd' : growth_params.getfloat('angle_sd')*(np.pi/180),
        'branch_mean' : growth_params.getfloat('branch_mean')*(np.pi/180),
        'branch_sd' : growth_params.getfloat('branch_sd')*(np.pi/180),
        'branch_cost' : growth_params.getfloat('branch_cost'),
        'branch_rate' : growth_params.getfloat('branch_rate'),
        'hy_density' :  growth_params.getfloat('hy_density'),
        'f_dw' :  growth_params.getfloat('f_dw'),
        'f_wall' :  growth_params.getfloat('f_wall'),
        'f_cw_cellwall' :  growth_params.getfloat('f_cw_cellwall'),

        # SECTION 3: Internal & External Nutrient Parameters
        'init_sub_e_dist' : nutrient_params['init_sub_e_dist'],
        'init_sub_e_gluc' : nutrient_params.getfloat('init_sub_e_gluc')*vol_grid,
        'init_sub_e_treha' : nutrient_params.getfloat('init_sub_e_treha')*vol_grid,
        'diffusion_e_gluc' : diff_e_gluc,
        'convert_metabolite': convert_metabolite,

        'init_sub_i_gluc' : nutrient_params.getfloat('init_sub_i_gluc'),
        'diffusion_i_gluc' : nutrient_params.getfloat('diffusion_i_gluc'),
        'vel_gluc' : nutrient_params.getfloat('vel_gluc'),
        'var_nutrient_backgrnd' : nutrient_params.getfloat('var_nutrient_backgrnd'),
        
        
        # 'm_gluc' : nutrient_params.getfloat('m_gluc'),
        # 'rho' : nutrient_params.getfloat('rho')*vol_seg,

        'ku1_gluc' : nutrient_params.getfloat('ku1_gluc'),
        'Ku2_gluc' : Ku2,
        'yield_u' : nutrient_params.getfloat('yield_u'),

        #'kc1_gluc' : nutrient_params.getfloat('kc1_gluc'),
        #'Kc2_gluc' : nutrient_params.getfloat('Kc2_gluc'),
        'kc1_gluc' : nutrient_params.getfloat('ku1_gluc'),
        'Kc2_gluc' : Ku2*init_vol_seg/vol_grid,
        'yield_c' : nutrient_params.getfloat('yield_c'),

        'kg1_wall' : nutrient_params.getfloat('kg1_wall'),
        # 'Kg2_wall' : nutrient_params.getfloat('Kg2_wall')*vol_seg,
        'Kg2_wall' : nutrient_params.getfloat('Kg2_wall'),
        'mw_cw' : nutrient_params.getfloat('mw_cw'),
        'mw_glucose' : nutrient_params.getfloat('mw_glucose'),
        'active_trsprt_vel_cw' : active_trsprt_vel_cw,

        'num_v' : nutrient_params.getfloat('num_v')
    }
    # The rate of glucose uptake is determined from the amount of glucose needed to support the growth rate.
    # The rate of uptake of glucose (kc1_gluc) is the product of the rate of radial growth (kg1_wall)
    # times the cross sectinoal area of the hyphae, the hyphal density, the fraction of wet cell mass that is dry cell mass
    # the fraction of dry cell mass that is cell wall material, the fraction of cell wall material that is composed of sugars (chitin and glucan),
    # the fraction of glucose mass that is provided by metabolism for cell wall raw material (yield_c),
    # and all divided by the formula weight of the chitin and glucan cell wall material.
    # The factor of 1.0e+03 is to convert the rate from moles to millimoles.
    kc1_gluc = params_dict['kg1_wall']*params_dict['cross_area']*params_dict['hy_density']\
                *params_dict['f_dw']*params_dict['f_wall']*params_dict['f_cw_cellwall'] \
                /(params_dict['mw_cw']*params_dict['yield_c'])*1.0e+03
    params_dict['kc1_gluc'] = kc1_gluc
    params_dict['ku1_gluc'] = kc1_gluc

    #if not('yield_c_in_mmoles' in params_dict):
    #	params_dict['yield_c_in_mmoles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
    params_dict['yield_c_in_mmoles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
    
    use_original = 0

    if(use_original != 1):
        # Max rate of moles of cell wall raw material used per time step:
        max_gms_cw_per_time = params_dict['kg1_wall']* np.pi*(params_dict['hy_diam']/2.0)**2.0 \
                        * params_dict['hy_density'] * params_dict['f_dw']* params_dict['f_wall'] * params_dict['f_cw_cellwall']
        max_moles_cw_per_time = max_gms_cw_per_time / params_dict['mw_cw']

        #Max rate of conversion of glucose to cell wall material is the maxrate of
        # grams of cell wall used per step/timestep / (gms cw produced per gms glucose used)/ (MW glucose)
        params_dict['yield_in_moles'] = params_dict['yield_c']*params_dict['mw_glucose']/params_dict['mw_cw']
        #params_dict['kc1_gluc'] =  max_moles_cw_per_time/yield_in_moles
        #params_dict['ku1_gluc'] = params_dict['kc1_gluc']


    return params_dict, config

# ----------------------------------------------------------------------------

def get_filepath(params):
    """
    Parameters
    ----------
    params : dict
        Dictionary of parameters to be used in simulation.

    Returns
    -------
    folder_string : str
        Path for the folder where results will be stored.
    param_string : str
        Filename to be used to label stored results, contains parameter info.

    """
    # folder_string = "ihc={}_ext={}_dy={}_sl={}_dt={:.3}_ft={}".format(
    #                                                         params['init_segs_count'],
    #                                                         params['init_sub_e_dist'],
    #                                                         params['dy'],
    #                                                         params['sl'],
    #                                                         params['dt'],
    #                                                         params['final_time'])
    # folder_string = "oldD2Tip_Fus_tipRe_brRate1e9_resBr4_noBkDiffLowGluc2_bkPatchy_Trsloc_4init"
    # folder_string = 'recalibration_02242022'
    #folder_string = "noFusion_tipRel_homogenousEnv_convert"
    #folder_string = "NoFusion_NoTipRel_homogenousEnv_initGluc20mm_branch0_3_brCost1x_seg=400"
    try:
        folder_string = params['output_path']
    except NameError:
        #folder_string = "NoFusion_NoTipRel_homogenousEnv_initGluc20mm_branch0_3_brCost1x_t1"
        folder_string = "test_Bill_branch3"
    # file_string = "{}_b={:.3e}_ieg={}_deg={}_iig={:.3e}_dig={}_vw={}_kyu={},{:.3e},{}_kyc={:.3e},{:.3e},{}_kyg={},{:.3e},{}".format(
    #     folder_string,
    #     params['branch_rate'],
    #     params['init_sub_e_gluc'], params['diffusion_e_gluc'],
    #     params['init_sub_i_gluc'], params['diffusion_i_gluc'],
    #     # params['vel_gluc'],
    #     params['vel_wall'],
    #     params['ku1_gluc'], params['Ku2_gluc'], params['yield_u'],
    #     params['kc1_gluc'], params['Kc2_gluc'], params['yield_c'],
    #     params['kg1_wall'], params['Kg2_wall'])
    #file_string = "NoFusion_tipRel_patch3Env_initGluc2um_branch0_3_brCost1x_t1"
    file_string = "NoFusion_AllHyphRelease_homogenousEnv_initGluc20mm_branch0_3_brCost1x_50x50x0.20umGrid"
    #file_string = "Fusion_AllHyphRelease_patchyEnv_initGluc20mm_branch0_3_brCost1x_200x200x0.20umRandomGrid"
    try:
        file_string = params['output_path']
    except NameError:
        #folder_string = "NoFusion_NoTipRel_homogenousEnv_initGluc20mm_branch0_3_brCost1x_t1"
        file_string = "test_Bill_branch3"
    
    return folder_string, file_string


# ----------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ----------------------------------------------------------------------------

def output_hyphal_coordinates(segments, hyphal_coord_file):
    thisfile = open(hyphal_coord_file, 'w')
    for i in range(np.shape(segments)[0]): 
        print(*segments[i],sep=', ',file=thisfile)
    thisfile.close()

def output_extern_concs(sub_e, extern_conc_file):
    thisfile = open(extern_conc_file, 'w')
    max_i, max_j = np.shape(sub_e) 
    #for i in range(max_i): 
    #    for j in range(max_j-1): 
    #        print(sub_e[i,j],sep = ', ',file=thisfile)
    #    print(sub_e[i,max_j-1],'\n',file=thisfile)
    for i in range(max_i): 
        print(sub_e[i],sep = ', ',file=thisfile)

    thisfile.close()

def plot_fungus(mycelia, num_total_segs, curr_time, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    hy : list
        List of class instances containing information about each hyphae segment.
    curr_time : double
        The current time of simulation, in days.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot fungal mycelia network.
    Color of a segment corresponds to internal substrate concentration.

    """
    # cur_len = len(hy)
    idx_to_display = np.where(mycelia['seg_vol'] > 0)[0]
    si_conc = mycelia['cw_i'][idx_to_display]/mycelia['seg_vol'][idx_to_display] *1.0e12
    

    #si = si_conc[idx_to_display].flatten()
    si = si_conc.flatten()
    
    x1 = mycelia['xy1'][idx_to_display, 0].tolist()
    x2 = mycelia['xy2'][idx_to_display, 0].tolist()
    y1 = mycelia['xy1'][idx_to_display, 1].tolist()
    y2 = mycelia['xy2'][idx_to_display, 1].tolist()


    # x1 = mycelia['xy1'][:num_total_segs, 0].tolist()
    # x2 = mycelia['xy2'][:num_total_segs, 0].tolist()
    # y1 = mycelia['xy1'][:num_total_segs, 1].tolist()
    # y2 = mycelia['xy2'][:num_total_segs, 1].tolist()
    # si = mycelia['cw_i'][:num_total_segs].flatten()

    if any(si < 1.0e-9):
        #min_value = min(si[(si > 1.0e-9)])
        si[np.where(si < 1.0e-9)] = 1.0e-09
    si = np.log10(si)


    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])
    
    segments_xyz_concs = []

    for xi1, yi1, xi2, yi2, concsi in zip(x1, y1, x2, y2, si_conc):
        segments_xyz_concs.append((xi1, yi1, xi2, yi2, concsi[0]))

    # Generated plot
    fig, ax = plt.subplots(dpi=600)

    top = cm.get_cmap('Oranges_r', 128) # r means reversed version
    bottom = cm.get_cmap('Blues', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Plot linesegments with coloring according to internal substrate conc.
    lc = mc.LineCollection(segments, array=si, cmap=cm.jet)#orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # plt.scatter(x1,y1,s=0.1)

    # Colorbar
    ax.add_collection(lc)
    fc = fig.colorbar(lc)
    fc.set_label('Cell Wall Components\n Log Conc. (Molar)')
    fc.outline.set_visible(False)

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time

    hyphal_coord_file = "Results/{}/Run{}/{}_t={:0.2f}_hyphal_coordinates_run{}.txt".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        run)
    output_hyphal_coordinates(segments_xyz_concs, hyphal_coord_file)
#    print(hyphal_coord_file)
#    thisfile = open(hyphal_coord_file, 'w')
#    print(type(thisfile))
#    coord_file = open(thisfile, 'w')
#    for i in range(np.shape(segments)[0]): 
#        print(segments[i],file=coord_file)
#    coord_file.close()
    # Set labels, title, margins, etc.
    # ax.set_ylabel('dm')
    # ax.set_xlabel('dm')
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    #ax.set_title('Mycelia Network \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    ax.set_title('Time = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.axis('equal')
    ax.margins(0.1)
    # breakpoint()

    # Show the plot
    sns.despine()
    #plt.show()

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_mycelia_cellwall_{}.png".format(param_string,
                                                                     run,
                                                                     param_string,
                                                                     curr_time,
                                                                     run)
    fig.savefig(fig_name)
    plt.close()

# ----------------------------------------------------------------------------

def plot_fungus_gluc(mycelia, num_total_segs, curr_time, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    hy : list
        List of class instances containing information about each hyphae segment.
    curr_time : double
        The current time of simulation, in days.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot fungal mycelia network.
    Color of a segment corresponds to internal substrate concentration.

    """
    # cur_len = len(hy)
    idx_to_display = np.where(mycelia['seg_vol'] > 0)[0]
    si_conc = mycelia['gluc_i'][idx_to_display]/mycelia['seg_vol'][idx_to_display] *1.0e12

    si = si_conc.flatten()
    x1 = mycelia['xy1'][idx_to_display, 0].tolist()
    x2 = mycelia['xy2'][idx_to_display, 0].tolist()
    y1 = mycelia['xy1'][idx_to_display, 1].tolist()
    y2 = mycelia['xy2'][idx_to_display, 1].tolist()
 
 
    
    #si = mycelia['gluc_i'][idx_to_display].flatten()
    # x1 = mycelia['xy1'][:num_total_segs, 0].tolist()
    # x2 = mycelia['xy2'][:num_total_segs, 0].tolist()
    # y1 = mycelia['xy1'][:num_total_segs, 1].tolist()
    # y2 = mycelia['xy2'][:num_total_segs, 1].tolist()
    # si = mycelia['gluc_i'][:num_total_segs].flatten()

    if any(si < 1.0e-9):
        #min_value = min(si[(si > 1.0e-9)])
        si[np.where(si < 1.0e-9)] = 1.0e-09
    si = np.log10(si)

    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])

    # Generated plot
    fig, ax = plt.subplots(dpi=600)

    top = cm.get_cmap('Oranges_r', 128) # r means reversed version
    bottom = cm.get_cmap('Blues', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Plot linesegments with coloring according to internal substrate conc.
    lc = mc.LineCollection(segments, array=si, cmap=cm.jet)#orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # plt.scatter(x1,y1,s=0.1)

    # Colorbar
    ax.add_collection(lc)
    fc = fig.colorbar(lc)
    fc.set_label('Glucose\n Log Conc. (Molar)')
    fc.outline.set_visible(False)

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time


    # Set labels, title, margins, etc.
    # ax.set_ylabel('dm')
    # ax.set_xlabel('dm')
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    #ax.set_title('Mycelia Network \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    ax.set_title('Time = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.axis('equal')
    ax.margins(0.1)
    # breakpoint()

    # Show the plot
    sns.despine()
    #plt.show()

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_mycelia_gluc_{}.png".format(param_string,
                                                                     run,
                                                                     param_string,
                                                                     curr_time,
                                                                     run)
    fig.savefig(fig_name)
    plt.close()

def plot_fungus_generic(mycelia, num_total_segs, curr_time, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    hy : list
        List of class instances containing information about each hyphae segment.
    curr_time : double
        The current time of simulation, in days.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot fungal mycelia network.
    Color of a segment corresponds to internal substrate concentration.

    """
    # cur_len = len(hy)

    x1 = mycelia['xy1'][:num_total_segs, 0].tolist()
    x2 = mycelia['xy2'][:num_total_segs, 0].tolist()
    y1 = mycelia['xy1'][:num_total_segs, 1].tolist()
    y2 = mycelia['xy2'][:num_total_segs, 1].tolist()
    ssi = mycelia['can_branch'][:num_total_segs].flatten()

    #si[np.where(si == 0.0)] = 1.0e-14
    si = np.zeros(ssi.shape)
    si[np.where(ssi == True)] = 500.0
    si[np.where(ssi == False)] = -500.0
    #breakpoint()

    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])

    # Generated plot
    fig, ax = plt.subplots(dpi=600)

    top = cm.get_cmap('Oranges_r', 128) # r means reversed version
    bottom = cm.get_cmap('Blues', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Plot linesegments with coloring according to internal substrate conc.
    lc = mc.LineCollection(segments, array=si, cmap=cm.jet)#orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # plt.scatter(x1,y1,s=0.1)

    # Colorbar
    ax.add_collection(lc)
    fc = fig.colorbar(lc)
    fc.set_label('Generic')
    fc.outline.set_visible(False)

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time


    # Set labels, title, margins, etc.
    # ax.set_ylabel('dm')
    # ax.set_xlabel('dm')
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    #ax.set_title('Mycelia Network \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    ax.set_title('Time = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.axis('equal')
    ax.margins(0.1)
    # breakpoint()

    # Show the plot
    sns.despine()
    #plt.show()

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_mycelia_gluc_{}.png".format(param_string,
                                                                     run,
                                                                     param_string,
                                                                     curr_time,
                                                                     run)
    fig.savefig(fig_name)
    plt.close()
    
def plot_fungus_treha(mycelia, num_total_segs, curr_time, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    hy : list
        List of class instances containing information about each hyphae segment.
    curr_time : double
        The current time of simulation, in days.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot fungal mycelia network.
    Color of a segment corresponds to internal substrate concentration.

    """
    idx_to_display = np.where(mycelia['seg_vol'] > 0)[0]
    si_conc = mycelia['treha_i'][idx_to_display]/mycelia['seg_vol'][idx_to_display] *1.0e12

    si = si_conc.flatten()
    x1 = mycelia['xy1'][idx_to_display, 0].tolist()
    x2 = mycelia['xy2'][idx_to_display, 0].tolist()
    y1 = mycelia['xy1'][idx_to_display, 1].tolist()
    y2 = mycelia['xy2'][idx_to_display, 1].tolist()

    # x1 = mycelia['xy1'][:num_total_segs, 0].tolist()
    # x2 = mycelia['xy2'][:num_total_segs, 0].tolist()
    # y1 = mycelia['xy1'][:num_total_segs, 1].tolist()
    # y2 = mycelia['xy2'][:num_total_segs, 1].tolist()
    # si = mycelia['gluc_i'][:num_total_segs].flatten()

    if any(si < 1.0e-9):
        #min_value = min(si[(si > 1.0e-9)])
        si[np.where(si < 1.0e-9)] = 1.0e-09
    si = np.log10(si)

    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])

    # Generated plot
    fig, ax = plt.subplots(dpi=600)

    top = cm.get_cmap('Oranges_r', 128) # r means reversed version
    bottom = cm.get_cmap('Blues', 128)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Plot linesegments with coloring according to internal substrate conc.
    #offs = (0.0, 0.0)
    #lc = mc.LineCollection(segments, offsets=offs, array=si, cmap=cm.jet)#orange_blue)
    lc = mc.LineCollection(segments, array=si, cmap=cm.jet)#orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # plt.scatter(x1,y1,s=0.1)

    # Colorbar
    ax.add_collection(lc)
    fc = fig.colorbar(lc)
    fc.set_label('Trehalose\n Log Conc. (Molar)')
    fc.outline.set_visible(False)

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time


    # Set labels, title, margins, etc.
    # ax.set_ylabel('dm')
    # ax.set_xlabel('dm')
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    #ax.set_title('Mycelia Network \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    ax.set_title('Time = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    ax.axis('equal')
    ax.margins(0.1)
    # breakpoint()

    # Show the plot
    sns.despine()
    #plt.show()

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_mycelia_treha_{}.png".format(param_string,
                                                                     run,
                                                                     param_string,
                                                                     curr_time,
                                                                     run)
    fig.savefig(fig_name)
    plt.close()

# ----------------------------------------------------------------------------

# def plot_externalsub(sub_e, yticks, yticklabels, curr_time, sub_e_max, plot_type, folder_string, param_string, params, run):
#     """
#     Parameters
#     ----------
#     sub_e : 2D numpy array
#         Matrix containing external nutrient concentration values at discritized grid points.
#     yticks : list
#         Helps determine how many labels appear of x- and y-axes.
#     yticklabels : list
#         Values to appear on the x- and y-axes.
#     curr_time : double
#         The current time of simulation, in days.
#     sub_e_max : double
#         Largest possible value for external substrate concentration.
#     param_string : str
#         Used to create filename of saved plot.

#     Returns
#     -------
#     None.

#     Purpose
#     -------
#     Plot external nutrient concentration.

#     """
#     fig, ax = pl.subplots(dpi=600)
#     # For the orange-blue color map
#     N = 1024#512#256
#     top = cm.get_cmap('Oranges_r', N)#128) # r means reversed version
#     bottom = cm.get_cmap('Blues', N)#128)# combine it all
#     newcolors = np.vstack((top(np.linspace(0, 1, N)),
#                            bottom(np.linspace(0, 1, N))))# create a new colormaps with a name of OrangeBlue
#     orange_blue = ListedColormap(newcolors, name='OrangeBlue')

#     # Convert units
#     if params['plot_units_time'] == 'days':
#         plot_time = curr_time / (60*60*24)
#     elif params['plot_units_time'] == 'hours':
#         plot_time = curr_time / (60*60)
#     elif params['plot_units_time'] == 'minutes':
#         plot_time = curr_time / 60
#     elif params['plot_units_time'] == 'seconds':
#         plot_time = curr_time

#     # Plot
#     if plot_type == 'Se':
#         ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, vmax=sub_e_max, xticklabels=yticklabels, yticklabels=yticklabels)
#     elif plot_type == 'Ce':
#         ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, xticklabels=yticklabels, yticklabels=yticklabels)
#     ax.set_yticks(yticks)
#     ax.set_xticks(yticks)
#     if plot_type == 'Se':
#         ax.collections[0].colorbar.set_label("External Nutrient Concentration")
#         ax.set_title('External Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
#     elif plot_type == 'Ce':
#         ax.collections[0].colorbar.set_label("Chemical Inhibitor Concentration")
#         ax.set_title('Chemical Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

#     ax.set_ylabel('{}'.format(params['plot_units_space']))
#     ax.set_xlabel('{}'.format(params['plot_units_space']))
#     ax.invert_yaxis()
#     ax.invert_xaxis()
#     ax.axis('equal')
#     ax.margins(0.1)
#     plt.show()
#     fig_name = "Results/{}/Run{}/{}_t={:0.2f}_external{}_{}_gluc.png".format(param_string,
#                                                                         run,
#                                                                         param_string,
#                                                                         curr_time,
#                                                                         plot_type,
#                                                                         run)
#     fig = ax.get_figure()
#     fig.savefig(fig_name)
    
#-----------------------------------------------------------------------------

def plot_externalsub(sub_e_orig, yticks, y_tick_labels, curr_time, sub_e_max, plot_type, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    sub_e : 2D numpy array
        Matrix containing external nutrient mole values at discritized grid points.
    yticks : list
        Helps determine how many labels appear of x- and y-axes.
    yticklabels : list
        Values to appear on the x- and y-axes.
    curr_time : double
        The current time of simulation, in days.
    sub_e_max : double
        Largest possible value for external substrate concentration.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot external nutrient concentration.

    """
    sub_e = sub_e_orig.copy()
    # Convert to molar quantities for display
    idx = np.where(sub_e > 0)
    nidx = np.where(sub_e <= 0)
    # Set concentrations to log concentrations
    sub_e[idx] = np.log10(sub_e[idx]/params['vol_grid']*1e12)
    # set the other values to the min value
    #if np.any(idx):
    #    sub_e[nidx] = np.min(sub_e[idx])-1
    # else:
    sub_e[nidx] = -20.0 # set log concentration to -30.0
    sub_e_max = np.max(sub_e)
    
    
    # For the orange-blue color map
    top = cm.get_cmap('Oranges_r', 256) # r means reversed version
    bottom = cm.get_cmap('Blues', 256)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 256)),
                           bottom(np.linspace(0, 1, 256))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time
    # breakpoint()
    # Plot
    if plot_type == 'Se':
        ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmax=sub_e_max, xticklabels=y_tick_labels, yticklabels=y_tick_labels)
    # elif plot_type == 'Ce':
    #     ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, xticklabels=yticklabels, yticklabels=yticklabels)
    # breakpoint()
    ax.set_yticks(yticks)
    ax.set_xticks(yticks)
    if plot_type == 'Se':
        ax.collections[0].colorbar.set_label("External Glucose\n Log Conc. (Molar)")
        ax.set_title('External Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    elif plot_type == 'Ce':
        ax.collections[0].colorbar.set_label("Chemical Inhibitor Concentration")
        ax.set_title('Chemical Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.set_xticklabels(y_tick_labels)
    ax.set(xticklabels=y_tick_labels)
    ax.set_yticklabels(y_tick_labels)
    ax.set(yticklabels=y_tick_labels)
    
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    ax.invert_yaxis()
    #ax.axis('equal')
    ax.margins(0.1)
    # ax.set(yticklabels=[])
    # ax.set(xticklabels=[])
    # ax.invert_xaxis()
    #plt.show()
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_external_gluc_{}_{}.png".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        plot_type,
                                                                        run)
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(fig_name, bbox_inches="tight")
    plt.close()

def plot_externalsub_hyphae(sub_e_orig, mycelia, num_total_segs, yticks, y_tick_labels, curr_time, sub_e_max, plot_type, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    sub_e : 2D numpy array
        Matrix containing external nutrient concentration values at discritized grid points.
    yticks : list
        Helps determine how many labels appear of x- and y-axes.
    yticklabels : list
        Values to appear on the x- and y-axes.
    curr_time : double
        The current time of simulation, in days.
    sub_e_max : double
        Largest possible value for external substrate concentration.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot external nutrient concentration.

    """
    # Convert to molar quantities for display
    sub_e = sub_e_orig.copy()
    idx = np.where(sub_e > 0)
    nidx = np.where(sub_e <= 0)
    sub_e[idx] = np.log10(sub_e[idx]/params['vol_grid']*1e12)
    # set the other values to the min value
    #if np.any(idx):
    #    sub_e[nidx] = np.min(sub_e[idx])-1
    # else:
    sub_e[nidx] = -20.0 # set log concentration to -30.0
    sub_e_max = np.max(sub_e)

    # For the orange-blue color map
    top = cm.get_cmap('Oranges_r', 256) # r means reversed version
    bottom = cm.get_cmap('Blues', 256)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 256)),
                           bottom(np.linspace(0, 1, 256))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time
    # breakpoint()
    # Plot
    if plot_type == 'Se':
        ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmax=sub_e_max, xticklabels=y_tick_labels, yticklabels=y_tick_labels)
    # elif plot_type == 'Ce':
    #     ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, xticklabels=yticklabels, yticklabels=yticklabels)
    # breakpoint()
    ax.set_yticks(yticks)
    ax.set_xticks(yticks)
    if plot_type == 'Se':
        ax.collections[0].colorbar.set_label("External Glucose\n Log Conc. (Molar)")
        ax.set_title('External Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    elif plot_type == 'Ce':
        ax.collections[0].colorbar.set_label("Chemical Inhibitor Concentration")
        ax.set_title('Chemical Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.set_xticklabels(y_tick_labels)
    ax.set(xticklabels=y_tick_labels)
    ax.set_yticklabels(y_tick_labels)
    ax.set(yticklabels=y_tick_labels)
    
    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    ax.invert_yaxis()
    #ax.axis('equal')
    ax.margins(0.1)
    # ax.set(yticklabels=[])
    # ax.set(xticklabels=[])
    # ax.invert_xaxis()
    #plt.show()

    # Now plot hyphae:
    ngrids = len(sub_e)
    max_xy = np.max(y_tick_labels)

    idx_to_display = np.where(mycelia['seg_vol'] > 0)[0]
    si_conc = mycelia['gluc_i'][idx_to_display]/mycelia['seg_vol'][idx_to_display] *1.0e12
    
    si = si_conc.flatten()
    x1 = (mycelia['xy1'][idx_to_display, 0]*ngrids/2/max_xy + ngrids/2).tolist()
    x2 = (mycelia['xy2'][idx_to_display, 0]*ngrids/2/max_xy + ngrids/2).tolist()
    y1 = (mycelia['xy1'][idx_to_display, 1]*ngrids/2/max_xy + ngrids/2).tolist()
    y2 = (mycelia['xy2'][idx_to_display, 1]*ngrids/2/max_xy + ngrids/2).tolist()

    if any(si < 1.0e-9):
        #min_value = min(si[(si > 1.0e-9)])
        si[np.where(si < 1.0e-9)] = 1.0e-09
    si = np.log10(si)
   
    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])

    # Plot linesegments with coloring according to internal substrate conc.
    #offset = (1.0, 1.0)
    #lc = mc.LineCollection(segments, offsets = offset, array=si, cmap=orange_blue)
    lc = mc.LineCollection(segments, array=si, cmap=orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # End plot of hyphae

    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_external_gluc_hyphae_{}_{}.png".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        plot_type,
                                                                        run)
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(fig_name, bbox_inches="tight")
    plt.close()


def plot_externalsub_treha(sub_e_orig, yticks, yticklabels, curr_time, sub_e_max, plot_type, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    sub_e : 2D numpy array
        Matrix containing external nutrient concentration values at discritized grid points.
    yticks : list
        Helps determine how many labels appear of x- and y-axes.
    yticklabels : list
        Values to appear on the x- and y-axes.
    curr_time : double
        The current time of simulation, in days.
    sub_e_max : double
        Largest possible value for external substrate concentration.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot external nutrient concentration.

    """
    # Convert to molar quantities for display
    #sub_e = np.log10(sub_e/params['vol_grid']*1e12) 
    #sub_e_max = np.max(sub_e[np.where(np.isfinite(sub_e))])

    #sub_e[np.where(np.isinf(sub_e))] = np.min(sub_e[np.where(np.isfinite(sub_e))])-1
    #sub_e[np.where(np.isinf(sub_e))] = 10*sub_e_max
    
    sub_e = sub_e_orig.copy()
    # Convert to molar quantities for display
    idx = np.where(sub_e > 0)
    nidx = np.where(sub_e <= 0)
    sub_e[idx] = np.log10(sub_e[idx]/params['vol_grid']*1e12)
    # set the other values to the min value
    # sub_e[nidx] = np.min(sub_e[idx])-1
    sub_e[nidx] = -20.0
    sub_e_max = np.max(sub_e)


    # For the orange-blue color map
    top = cm.get_cmap('Oranges_r', 256) # r means reversed version
    bottom = cm.get_cmap('Blues', 256)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 256)),
                           bottom(np.linspace(0, 1, 256))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time

    extern_conc_file = "Results/{}/Run{}/{}_t={:0.2f}_external_logConcentrations_run{}.txt".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        run)
    output_extern_concs(sub_e, extern_conc_file)
    # breakpoint()
    # Plot
    if plot_type == 'Se':
        #ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)
        #ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)
        ax = sns.heatmap(sub_e, cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)

    # elif plot_type == 'Ce':
    #     ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, xticklabels=yticklabels, yticklabels=yticklabels)
    # breakpoint()
    ax.set_yticks(yticks)
    ax.set_xticks(yticks)
    
    ax.set_yticklabels(yticklabels)
    ax.set_xticklabels(yticklabels)
    if plot_type == 'Se':
        ax.collections[0].colorbar.set_label("External Trehalose\n Log Conc. (Molar)")
        ax.set_title('External Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    elif plot_type == 'Ce':
        ax.collections[0].colorbar.set_label("Chemical Inhibitor Concentration")
        ax.set_title('Chemical Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    ##ax.invert_yaxis()
    #ax.axis('equal')
    ax.margins(1.9)
    # ax.set(yticklabels=[])
    # ax.set(xticklabels=[])
    # ax.invert_xaxis()
    #plt.show()
    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_external{}_{}.png".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        plot_type,
                                                                        run)
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(fig_name, bbox_inches="tight")
    plt.close(fig)
# ----------------------------------------------------------------------------

def plot_externalsub_treha_hyphae(sub_e_orig, mycelia, num_total_segs, yticks, yticklabels, curr_time, sub_e_max, plot_type, folder_string, param_string, params, run):
    """
    Parameters
    ----------
    sub_e : 2D numpy array
        Matrix containing external nutrient concentration values at discritized grid points.
    yticks : list
        Helps determine how many labels appear of x- and y-axes.
    yticklabels : list
        Values to appear on the x- and y-axes.
    curr_time : double
        The current time of simulation, in days.
    sub_e_max : double
        Largest possible value for external substrate concentration.
    param_string : str
        Used to create filename of saved plot.

    Returns
    -------
    None.

    Purpose
    -------
    Plot external nutrient concentration.

    """
    # Convert to molar quantities for display
    sub_e = sub_e_orig.copy()
    idx = np.where(sub_e > 0)
    nidx = np.where(sub_e <= 0)
    sub_e[idx] = np.log10(sub_e[idx]/params['vol_grid']*1e12)
    # set the other values to the min value
    #if np.any(idx):
    #    sub_e[nidx] = np.min(sub_e[idx])-1
    #else:
    sub_e[nidx] = -20.0 # set log concentration to -30.0
    sub_e_max = np.max(sub_e)
    



    # For the orange-blue color map
    top = cm.get_cmap('Oranges_r', 256) # r means reversed version
    bottom = cm.get_cmap('Blues', 256)# combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 256)),
                           bottom(np.linspace(0, 1, 256))))# create a new colormaps with a name of OrangeBlue
    orange_blue = ListedColormap(newcolors, name='OrangeBlue')

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time

    # breakpoint()
    # Plot
    if plot_type == 'Se':
        #ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)
        ax = sns.heatmap(np.rot90(sub_e, k=1, axes=(0,1)), cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)
        #ax = sns.heatmap(sub_e, cmap=orange_blue, vmax=sub_e_max)#, xticklabels=yticklabels, yticklabels=yticklabels)
    # elif plot_type == 'Ce':
    #     ax = sns.heatmap(np.transpose(sub_e), cmap=orange_blue, vmin=0, xticklabels=yticklabels, yticklabels=yticklabels)
    # breakpoint()
    ax.set_yticks(yticks)
    ax.set_xticks(yticks)
    
    ax.set_yticklabels(yticklabels)
    ax.set_xticklabels(yticklabels)
    if plot_type == 'Se':
        ax.collections[0].colorbar.set_label("External Trehalose\n Log Conc. (Molar)")
        ax.set_title('External Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")
    elif plot_type == 'Ce':
        ax.collections[0].colorbar.set_label("Chemical Inhibitor Concentration")
        ax.set_title('Chemical Domain \nTime = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    ax.set_ylabel('{}'.format(params['plot_units_space']))
    ax.set_xlabel('{}'.format(params['plot_units_space']))
    ##ax.invert_yaxis()
    #ax.axis('equal')
    ax.margins(1.9)
    # ax.set(yticklabels=[])
    # ax.set(xticklabels=[])
    # ax.invert_xaxis()
    #plt.show()

    # Now plot hyphae:
    ngrids = len(sub_e)
    max_xy = np.max(yticklabels)

    idx_to_display = np.where(mycelia['seg_vol'] > 0)[0]
    si_conc = mycelia['treha_i'][idx_to_display]/mycelia['seg_vol'][idx_to_display] *1.0e12

    si = si_conc.flatten()
    x1 = (mycelia['xy1'][idx_to_display, 0]*ngrids/2/max_xy + ngrids/2).tolist()
    x2 = (mycelia['xy2'][idx_to_display, 0]*ngrids/2/max_xy + ngrids/2).tolist()
    y1 = (mycelia['xy1'][idx_to_display, 1]*ngrids/2/max_xy + ngrids/2).tolist()
    y2 = (mycelia['xy2'][idx_to_display, 1]*ngrids/2/max_xy + ngrids/2).tolist()

    if any(si < 1.0e-9):
        #min_value = min(si[(si > 1.0e-9)])
        si[np.where(si < 1.0e-9)] = 1.0e-09
    si = np.log10(si)
    
    segments = []
    for xi1, yi1, xi2, yi2 in zip(x1, y1, x2, y2):
        segments.append([(xi1, yi1), (xi2, yi2)])

    # Plot linesegments with coloring according to internal substrate conc.
    #offset = (1.0, 1.0)
    #lc = mc.LineCollection(segments, offsets = offset, array=si, cmap=orange_blue)
    lc = mc.LineCollection(segments, array=si, cmap=orange_blue)
    lc.set_linewidth(1)
    ax.add_collection(lc)

    # End plot of hyphae

    fig_name = "Results/{}/Run{}/{}_t={:0.2f}_external{}_hyphae_{}.png".format(param_string,
                                                                        run,
                                                                        param_string,
                                                                        curr_time,
                                                                        plot_type,
                                                                        run)
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(fig_name, bbox_inches="tight")
    plt.close(fig)
# ----------------------------------------------------------------------------

def plot_stat(count_times, count_stat, stat_type, folder_string, param_string, params, run):

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_times = count_times
    elif params['plot_units_time'] == 'hours':
        plot_times = 24*count_times
    elif params['plot_units_time'] == 'minutes':
        plot_times = 60*24*count_times
    elif params['plot_units_time'] == 'seconds':
        plot_times = 60*60*24*count_times

    fig, ax = plt.subplots()
    ax.plot(plot_times, count_stat)
    # ax.set_title(stat_type)
    ax.set_xlabel('Time ({})'.format(params['plot_units_time']))
    ax.set_ylabel(stat_type)
    sns.despine()
    #plt.show()

    if stat_type == 'Num. of Branches':
        key_word = 'stat_b'
    elif stat_type == 'Num. of Tips':
        key_word = 'stat_t'
    elif stat_type == 'Branching Density':
        key_word = 'stat_d'
    elif stat_type == 'Radii of Mycelia ({})'.format(params['plot_units_space']):
        key_word = 'stat_r'
    elif stat_type == 'RMS Radii of Mycelia ({})'.format(params['plot_units_space']):
        key_word = 'stat_rms_r'
    elif stat_type == 'RMS Radii of Mycelia Tips ({})'.format(params['plot_units_space']):
        key_word = 'stat_rms_r'
    # Save the plot
    fig_name = "Results/{}/Run{}/{}_{}_{}.png".format(param_string,
                                                      run,
                                                      param_string,
                                                      key_word,
                                                      run)
    fig.savefig(fig_name)
    plt.close()

# ----------------------------------------------------------------------------

def plot_avg_treha_annulus(count_stat,max_count_stat,min_count_stat, stat_type, folder_string, param_string, curr_time, params, run):

    fig, ax = plt.subplots()
    xlabel = range(30,30*len(count_stat),30)
    # breakpoint()
    ax.plot(xlabel,count_stat[1:]/max_count_stat)
    # ax.set_title(stat_type)
    # ax.set_xlabel('Time ({})'.format(params['plot_units_time']))
    ax.set_ylabel(stat_type)
    sns.despine()
    #plt.show()

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time


    max_conc = (max_count_stat/params['vol_grid']*1e12)
    min_conc = (min_count_stat/params['vol_grid']*1e12) 
    ax.set_title('Max conc = {:0.2e} Min conc = {:0.2e}'.format(max_conc, min_conc),fontweight="bold")
    plt.suptitle('Time = {:0.2f} {}'.format(plot_time, params['plot_units_time']),fontweight="bold")

    if stat_type == 'Num. of Branches':
        key_word = 'stat_b'
    elif stat_type == 'Num. of Tips':
        key_word = 'stat_t'
    elif stat_type == 'Branching Density':
        key_word = 'stat_d'
    elif stat_type == 'Radii of Mycelia ({})'.format(params['plot_units_space']):
        key_word = 'stat_r'

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_{}_{}_avgTrehaAnnulus.png".format(param_string,
                                                      run,
                                                      param_string,
                                                      curr_time,
                                                      run)
    fig.savefig(fig_name)
    plt.close()
    
def plot_max_treha_annulus(count_stat,max_count_stat, stat_type, folder_string, param_string, current_time, params, run):

    fig, ax = plt.subplots()
    xlabel = range(30,30*len(count_stat),30)
    # breakpoint()
    ax.plot(xlabel,count_stat[1:]/max_count_stat)
    # ax.set_title(stat_type)
    # ax.set_xlabel('Time ({})'.format(params['plot_units_time']))
    ax.set_ylabel(stat_type)
    sns.despine()
    #plt.show()

    if stat_type == 'Num. of Branches':
        key_word = 'stat_b'
    elif stat_type == 'Num. of Tips':
        key_word = 'stat_t'
    elif stat_type == 'Branching Density':
        key_word = 'stat_d'
    elif stat_type == 'Radii of Mycelia ({})'.format(params['plot_units_space']):
        key_word = 'stat_r'

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_{}_{}_maxTrehaAnnulus.png".format(param_string,
                                                      run,
                                                      param_string,
                                                      current_time,
                                                      run)
    fig.savefig(fig_name)
    plt.close()

def plot_min_treha_annulus(count_stat,max_count_stat, stat_type, folder_string, param_string, current_time, params, run):

    fig, ax = plt.subplots()
    xlabel = range(30,30*len(count_stat),30)
    # breakpoint()
    ax.plot(xlabel,count_stat[1:]/max_count_stat)
    # ax.set_title(stat_type)
    # ax.set_xlabel('Time ({})'.format(params['plot_units_time']))
    ax.set_ylabel(stat_type)
    sns.despine()
    #plt.show()

    if stat_type == 'Num. of Branches':
        key_word = 'stat_b'
    elif stat_type == 'Num. of Tips':
        key_word = 'stat_t'
    elif stat_type == 'Branching Density':
        key_word = 'stat_d'
    elif stat_type == 'Radii of Mycelia ({})'.format(params['plot_units_space']):
        key_word = 'stat_r'

    # Save the plot
    fig_name = "Results/{}/Run{}/{}_{}_{}_minTrehaAnnulus.png".format(param_string,
                                                      run,
                                                      param_string,
                                                      current_time,
                                                      run)
    fig.savefig(fig_name)
    
# ----------------------------------------------------------------------------

def plot_errorbar_stat(count_times, avg_stat, std_stat, stat_type, folder_string, param_string, params, num_runs):

    # Convert units
    if params['plot_units_time'] == 'days':
        plot_times = count_times
    elif params['plot_units_time'] == 'hours':
        plot_times = 24*count_times
    elif params['plot_units_time'] == 'minutes':
        plot_times = 60*24*count_times
    elif params['plot_units_time'] == 'seconds':
        plot_times = 60*60*24*count_times

    fig, ax = plt.subplots()
    ax.errorbar(plot_times, avg_stat, std_stat)
    # ax.set_title(stat_type)
    ax.set_xlabel('Time ({})'.format(params['plot_units_time']))
    ax.set_ylabel(stat_type)
    sns.despine()
    #plt.show()

    if stat_type == 'Avg. Num. of Branches ({} Iterations)'.format(num_runs):
        key_word = 'avg_b'
    elif stat_type == 'Avg. Num. of Tips ({} Iterations)'.format(num_runs):
        key_word = 'avg_t'
    elif stat_type == 'Avg. Branching Density ({} Iterations)'.format(num_runs):
        key_word = 'avg_d'
    elif stat_type == 'Avg. Radii in {} ({} Iterations)'.format(params['plot_units_space'], num_runs):
        key_word = 'avg_r'

    # Save the plot
    fig_name = "Results/{}/Avg{}/{}_{}_avg{}.png".format(param_string,
                                                         num_runs,
                                                         param_string,
                                                         key_word,
                                                         num_runs)
    fig.savefig(fig_name)
    plt.close()

# ----------------------------------------------------------------------------

def plot_biomassdensity(radius_i, biomass_density, curr_time):
    """
    Parameters
    ----------
    radius_i : list
        List of smaller annuli radii value in which density is commputed.
    biomass_density : list
        Density of hyphae segment in an annulus with inner radii corresponding to radius_i.
    curr_time : double
        The current time of simulation, in days.

    Returns
    -------
    None.

    Purpose
    -------
    Plot biomass density at different distances from center of colony.
    For an annulus with inner radius r1 and outer radius r2,
        biomass density = (num. of segments in annulus) / (pi*(r2^2-r1^2))

    """
    fig1, ax1 = plt.subplots()
    ax1.plot(radius_i[1:len(radius_i)], biomass_density[1:len(radius_i)], marker='o')
    ax1.set_title('Biomass Density \nTime = {:.1f} Days'.format(curr_time), fontweight='bold')
    ax1.set_xlabel('Distance To Center (mm)')
    ax1.set_ylabel('Biomass Density')
    sns.despine()
    #plt.show()
    plt.close()


def plot_tipdensity(radius_i, tip_density, curr_time):
    """
    Parameters
    ----------
    radius_i : list
        List of smaller annuli radii value in which density is commputed.
    tip_density : list
        Density of hyphae tips in an annulus with inner radii corresponding to radius_i.
    curr_time : double
        The current time of simulation, in days.
    scale_val : double
        Parameter descrbing units used, scale_val=1 for mm or scale_val=1000 for µm.

    Returns
    -------
    None.

    Purpose
    -------
    Plot hyphal tip density at different distances from center of colony.
    For an annulus with inner radius r1 and outer radius r2,
        tip density = (num. of tips in annulus) / (pi*(r2^2-r1^2))

    """
    fig2, ax2 = plt.subplots()
    ax2.plot(radius_i[1:len(radius_i)], tip_density[1:len(radius_i)], marker='o')
    ax2.set_title('Hyphal Tip Density \nTime = {:.1f} Days'.format(curr_time), fontweight='bold')
    ax2.set_xlabel('Distance To Center (mm)')
    ax2.set_ylabel('Hyphal Tip Density')
    sns.despine()
    #plt.show()
    plt.close()

##############################################################################

def plot_hist(mycelia, curr_time,num_total_segs, param_string, params, run):
    
    if params['plot_units_time'] == 'days':
        plot_time = curr_time / (60*60*24)
    elif params['plot_units_time'] == 'hours':
        plot_time = curr_time / (60*60)
    elif params['plot_units_time'] == 'minutes':
        plot_time = curr_time / 60
    elif params['plot_units_time'] == 'seconds':
        plot_time = curr_time
    fig, ax = plt.subplots()
    # breakpoint()
    ax.hist(mycelia['dist_from_center'][:num_total_segs], range=[0, 1000], bins=100)
    fig_name = "Results/{}/Run{}/{}_{}_{}.png".format(param_string,
                                                      run,
                                                      curr_time,
                                                      param_string,                                                    
                                                      run)
    fig.savefig(fig_name)
    
##############################################################################

def plot_density_annulus(density_per_unit_annulus, num_total_segs, param_string, params, run):
    
    fig, ax = plt.subplots()
    ax.plot(range(2000), density_per_unit_annulus)
    fig_name = "Results/{}/Run{}/{}_{}_density.png".format(param_string,
                                                      run,
                                                      param_string,                                                    
                                                      run)
    fig.savefig(fig_name)
    
##############################################################################

def plot_treha_conc_annulus(avg_treha_annulus, num_total_segs, param_string, current_time, params, run):
    
    fig, ax = plt.subplots()
    ax.plot(range(2000), avg_treha_annulus)
    fig_name = "Results/{}/Run{}/{}_{}_{}_treha_conc.png".format(param_string,
                                                      run,
                                                      param_string,
                                                      current_time,
                                                      run)
    fig.savefig(fig_name)


# ----------------------------------------------------------------------------
# BIOMASS TRACKING FUNCTIONS
# ----------------------------------------------------------------------------

def calculate_total_biomass(mycelia, num_total_segs, params):
    """
    Calculate total biomass of the fungal colony.
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
    num_total_segs : int
        Current total number of segments.
    params : dict
        Parameter dictionary containing density and cross-sectional area.
    
    Returns
    -------
    total_biomass : float
        Total biomass in grams.
    """
    # Get all active segments (non-null segments)
    non_null_segs = np.where(mycelia['branch_id'][:num_total_segs] > -1)[0]
    
    # Sum of all segment lengths
    total_length = np.sum(mycelia['seg_length'][non_null_segs])
    
    # Volume = length * cross-sectional area
    total_volume = total_length * params['cross_area']
    
    # Biomass = volume * density
    # hy_density is in g/cm³, cross_area is in cm², seg_length is in cm
    total_biomass = total_volume * params['hy_density']
    
    return total_biomass


def plot_biomass_vs_time(time_array, biomass_array, folder_string, param_string, params, run):
    """
    Plot total biomass over time.
    
    Parameters
    ----------
    time_array : array
        Time values in seconds.
    biomass_array : array
        Corresponding biomass values in grams.
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for filename.
    params : dict
        Parameter dictionary.
    run : int
        Run number.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Convert time units
    if params['plot_units_time'] == 'days':
        plot_times = time_array / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = time_array / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = time_array / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = time_array
        time_label = 'Time (seconds)'
    
    ax.plot(plot_times, biomass_array, linewidth=2, marker='o', markersize=4, color='#2E7D32')
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel('Total Biomass (g)', fontsize=12, fontweight='bold')
    ax.set_title('Total Fungal Biomass vs Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add text with final biomass
    final_biomass = biomass_array[-1]
    ax.text(0.95, 0.05, f'Final: {final_biomass:.6f} g', 
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_biomass_vs_time_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Biomass plot saved to: {fig_name}")
    print(f"Final biomass: {final_biomass:.6f} g")


def plot_biomass_comparison(all_times, all_biomass, run_labels, folder_string, param_string, params):
    """
    Compare biomass growth across multiple runs.
    
    Parameters
    ----------
    all_times : list of arrays
        Time arrays for each run.
    all_biomass : list of arrays
        Biomass arrays for each run.
    run_labels : list of str
        Labels for each run/condition.
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for filename.
    params : dict
        Parameter dictionary.
    """
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_times)))
    
    for i, (times, biomass, label) in enumerate(zip(all_times, all_biomass, run_labels)):
        if params['plot_units_time'] == 'days':
            plot_times = times / (60*60*24)
            time_label = 'Time (days)'
        elif params['plot_units_time'] == 'hours':
            plot_times = times / (60*60)
            time_label = 'Time (hours)'
        elif params['plot_units_time'] == 'minutes':
            plot_times = times / 60
            time_label = 'Time (minutes)'
        else:
            plot_times = times
            time_label = 'Time (seconds)'
        
        ax.plot(plot_times, biomass, linewidth=2, marker='o', 
                markersize=3, label=label, color=colors[i])
    
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel('Total Biomass (g)', fontsize=12, fontweight='bold')
    ax.set_title('Biomass Comparison Across Conditions', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Print which grew most
    final_biomasses = [biomass[-1] for biomass in all_biomass]
    max_idx = np.argmax(final_biomasses)
    min_idx = np.argmin(final_biomasses)
    
    # Add text box with statistics
    stats_text = f'Highest: {run_labels[max_idx]}\n({final_biomasses[max_idx]:.6f} g)\n'
    stats_text += f'Lowest: {run_labels[min_idx]}\n({final_biomasses[min_idx]:.6f} g)'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            ha='left', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    num_runs = len(all_times)
    fig_name = f"Results/{param_string}/Avg{num_runs}/{param_string}_biomass_comparison_avg{num_runs}.png"
    fig.savefig(fig_name)
    plt.close()
    
    print(f"\nBiomass comparison plot saved to: {fig_name}")
    print(f"Highest final biomass: {run_labels[max_idx]} with {final_biomasses[max_idx]:.6f} g")
    print(f"Lowest final biomass: {run_labels[min_idx]} with {final_biomasses[min_idx]:.6f} g")


def plot_errorbar_biomass(count_times, avg_biomass, std_biomass, folder_string, param_string, params, num_runs):
    """
    Plot average biomass with error bars across multiple runs.
    
    Parameters
    ----------
    count_times : array
        Time values.
    avg_biomass : array
        Average biomass values.
    std_biomass : array
        Standard deviation of biomass values.
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for filename.
    params : dict
        Parameter dictionary.
    num_runs : int
        Number of runs averaged.
    """
    # Convert units
    if params['plot_units_time'] == 'days':
        plot_times = count_times / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = count_times / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = count_times / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = count_times
        time_label = 'Time (seconds)'

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.errorbar(plot_times, avg_biomass, std_biomass, linewidth=2, 
                marker='o', markersize=4, capsize=5, color='#2E7D32')
    
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Avg. Total Biomass ({num_runs} Iterations) (g)', fontsize=12, fontweight='bold')
    ax.set_title(f'Average Biomass vs Time ({num_runs} Runs)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add text with final values
    final_avg = avg_biomass[-1]
    final_std = std_biomass[-1]
    ax.text(0.95, 0.05, f'Final: {final_avg:.6f} ± {final_std:.6f} g', 
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()

    # Save the plot
    fig_name = "Results/{}/Avg{}/{}_avg_biomass_avg{}.png".format(
        param_string, num_runs, param_string, num_runs)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Average biomass plot saved to: {fig_name}")
    print(f"Final average biomass: {final_avg:.6f} ± {final_std:.6f} g")


# ----------------------------------------------------------------------------
# LONGEST BRANCH TRACKING FUNCTIONS
# ----------------------------------------------------------------------------

def calculate_longest_branch(mycelia, num_total_segs):
    """
    Calculate the length of the longest continuous branch.
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
    num_total_segs : int
        Current total number of segments.
    
    Returns
    -------
    longest_branch_length : float
        Length of the longest branch.
    """
    # Get all active segments
    non_null_segs = np.where(mycelia['branch_id'][:num_total_segs] > -1)[0]
    
    if len(non_null_segs) == 0:
        return 0.0
    
    # Get unique branch IDs
    branch_ids = np.unique(mycelia['branch_id'][non_null_segs])
    
    longest_branch_length = 0.0
    
    # For each branch, sum the lengths of all its segments
    for branch_id in branch_ids:
        branch_segs = np.where(mycelia['branch_id'][:num_total_segs] == branch_id)[0]
        branch_length = np.sum(mycelia['seg_length'][branch_segs])
        if branch_length > longest_branch_length:
            longest_branch_length = branch_length
    
    return longest_branch_length


def calculate_Leuc_to_farthest_tip(mycelia, num_total_segs, origin=None):
    """
    Calculate L_euc: the Euclidean (straight-line) distance from the origin 
    to the farthest hyphal tip.
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
        Required keys:
        - 'branch_id': segment branch IDs (-1 indicates inactive/null segment)
        - 'is_tip': boolean array indicating which segments are tips
        - 'xy2': endpoint coordinates of each segment
    num_total_segs : int
        Current total number of segments in the mycelium.
    origin : tuple or ndarray, optional
        (x, y) coordinates of the origin point. If None, uses (0, 0).
        Default is None.
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'Leuc_max': float, the Euclidean distance to the farthest tip
        - 'farthest_tip_idx': int, segment index of the farthest tip
        - 'farthest_tip_xy': ndarray, (x, y) coordinates of the farthest tip
        - 'all_tip_Leuc': dict, mapping tip indices to their Euclidean distances
        - 'num_tips': int, total number of tips found
    
    Notes
    -----
    This computes the straight-line distance, which is always less than or equal 
    to L_net (the network path distance). The ratio L_net/L_euc gives the tortuosity.
    """
    # Set default origin
    if origin is None:
        origin = np.array([0.0, 0.0])
    else:
        origin = np.array(origin)
    
    # Get all active segments (branch_id > -1)
    active_segs = np.where(mycelia['branch_id'][:num_total_segs].flatten() > -1)[0]
    
    if len(active_segs) == 0:
        return {
            'Leuc_max': 0.0,
            'farthest_tip_idx': -1,
            'farthest_tip_xy': origin.copy(),
            'all_tip_Leuc': {},
            'num_tips': 0
        }
    
    # Find all tip segments
    tip_mask = mycelia['is_tip'][:num_total_segs].flatten().astype(bool)
    tip_indices = np.where(tip_mask & (mycelia['branch_id'][:num_total_segs].flatten() > -1))[0]
    
    if len(tip_indices) == 0:
        return {
            'Leuc_max': 0.0,
            'farthest_tip_idx': -1,
            'farthest_tip_xy': origin.copy(),
            'all_tip_Leuc': {},
            'num_tips': 0
        }
    
    # Calculate Euclidean distance from origin to each tip
    all_tip_Leuc = {}
    max_Leuc = 0.0
    farthest_tip_idx = -1
    
    for tip_idx in tip_indices:
        tip_xy = mycelia['xy2'][tip_idx]
        # Euclidean distance: sqrt((x - x0)^2 + (y - y0)^2)
        Leuc = np.sqrt((tip_xy[0] - origin[0])**2 + (tip_xy[1] - origin[1])**2)
        all_tip_Leuc[tip_idx] = Leuc
        
        if Leuc > max_Leuc:
            max_Leuc = Leuc
            farthest_tip_idx = tip_idx
    
    # Get coordinates of farthest tip
    if farthest_tip_idx >= 0:
        farthest_tip_xy = mycelia['xy2'][farthest_tip_idx].copy()
    else:
        farthest_tip_xy = origin.copy()
    
    return {
        'Leuc_max': max_Leuc,
        'farthest_tip_idx': farthest_tip_idx,
        'farthest_tip_xy': farthest_tip_xy,
        'all_tip_Leuc': all_tip_Leuc,
        'num_tips': len(tip_indices)
    }


def calculate_Lnet_to_farthest_tip(mycelia, num_total_segs, origin_segments=None):
    """
    Calculate L_net: the shortest path distance through the weighted hyphal network
    from the origin (inoculum) to the farthest tip.
    
    For non-fusing mycelia (trees): There is only one unique path from the inoculum 
    to any tip, so L_net is simply the sum of segment lengths in the lineage.
    
    For fusing mycelia (networks): Because loops exist, there are multiple ways to 
    reach a tip. This function uses Dijkstra's algorithm to find the shortest path.
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
        Required keys:
        - 'branch_id': segment branch IDs (-1 indicates inactive/null segment)
        - 'nbr_idxs': list of neighbor segment indices for each segment
        - 'seg_length': length of each segment (edge weights)
        - 'is_tip': boolean array indicating which segments are tips
        - 'xy2': endpoint coordinates of each segment
    num_total_segs : int
        Current total number of segments in the mycelium.
    origin_segments : list or None, optional
        Indices of origin segments (inoculum location). If None, uses segments 
        closest to (0,0) as origin. Default is None.
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'Lnet_max': float, the L_net distance to the farthest reachable tip
        - 'farthest_tip_idx': int, segment index of the farthest tip
        - 'farthest_tip_xy': ndarray, (x,y) coordinates of the farthest tip
        - 'euclidean_distance': float, straight-line distance from origin to farthest tip
        - 'tortuosity': float, ratio of Lnet to Euclidean distance (>= 1.0)
        - 'all_tip_Lnet': dict, mapping tip indices to their L_net values
        - 'num_tips': int, total number of tips found
        - 'num_reachable_tips': int, number of tips reachable from origin
    
    Notes
    -----
    Uses Dijkstra's algorithm for shortest path computation. Edge weights are the
    segment lengths. The graph is constructed treating each segment as a node,
    with edges to its neighbors (from 'nbr_idxs').
    """
    import heapq
    
    # Get all active segments (branch_id > -1)
    active_segs = np.where(mycelia['branch_id'][:num_total_segs].flatten() > -1)[0]
    
    if len(active_segs) == 0:
        return {
            'Lnet_max': 0.0,
            'farthest_tip_idx': -1,
            'farthest_tip_xy': np.array([0.0, 0.0]),
            'euclidean_distance': 0.0,
            'tortuosity': 1.0,
            'all_tip_Lnet': {},
            'num_tips': 0,
            'num_reachable_tips': 0
        }
    
    # Determine origin segments (inoculum location)
    if origin_segments is None:
        # Find segments closest to origin (0, 0)
        # Use xy1 coordinates to find the starting point
        distances_to_origin = np.sqrt(mycelia['xy1'][active_segs, 0]**2 + 
                                       mycelia['xy1'][active_segs, 1]**2)
        min_dist = np.min(distances_to_origin)
        # Origin segments are those very close to the minimum distance point
        tolerance = 1e-6
        origin_mask = distances_to_origin <= min_dist + tolerance
        origin_segments = active_segs[origin_mask].tolist()
        
        # If no clear origin, use segment 0 if it's active
        if len(origin_segments) == 0:
            if 0 in active_segs:
                origin_segments = [0]
            else:
                origin_segments = [active_segs[0]]
    
    # Find all tip segments
    tip_mask = mycelia['is_tip'][:num_total_segs].flatten().astype(bool)
    tip_indices = np.where(tip_mask & (mycelia['branch_id'][:num_total_segs].flatten() > -1))[0]
    
    if len(tip_indices) == 0:
        return {
            'Lnet_max': 0.0,
            'farthest_tip_idx': -1,
            'farthest_tip_xy': np.array([0.0, 0.0]),
            'euclidean_distance': 0.0,
            'tortuosity': 1.0,
            'all_tip_Lnet': {},
            'num_tips': 0,
            'num_reachable_tips': 0
        }
    
    # -------------------------------------------------------------------------
    # Dijkstra's Algorithm for shortest paths from origin to all segments
    # -------------------------------------------------------------------------
    # Distance dictionary: dist[seg_idx] = shortest path distance from origin
    dist = {seg: float('inf') for seg in active_segs}
    
    # Priority queue: (distance, segment_index)
    pq = []
    
    # Initialize origin segments with distance = their own segment length
    # (the path to reach the end of the origin segment is its own length)
    for origin_seg in origin_segments:
        seg_len = float(mycelia['seg_length'][origin_seg])
        dist[origin_seg] = seg_len
        heapq.heappush(pq, (seg_len, origin_seg))
    
    # Process the priority queue
    while pq:
        current_dist, current_seg = heapq.heappop(pq)
        
        # Skip if we've already found a shorter path
        if current_dist > dist[current_seg]:
            continue
        
        # Get neighbors of current segment
        neighbors = mycelia['nbr_idxs'][current_seg]
        if neighbors is None:
            continue
        
        for neighbor in neighbors:
            # Skip inactive segments
            if neighbor >= num_total_segs:
                continue
            if mycelia['branch_id'][neighbor] <= -1:
                continue
            
            # Calculate new distance: current path + neighbor's segment length
            neighbor_len = float(mycelia['seg_length'][neighbor])
            new_dist = current_dist + neighbor_len
            
            # Update if we found a shorter path
            if new_dist < dist.get(neighbor, float('inf')):
                dist[neighbor] = new_dist
                heapq.heappush(pq, (new_dist, neighbor))
    
    # -------------------------------------------------------------------------
    # Find the farthest tip and collect all tip L_net values
    # -------------------------------------------------------------------------
    all_tip_Lnet = {}
    max_Lnet = 0.0
    farthest_tip_idx = -1
    num_reachable = 0
    
    for tip_idx in tip_indices:
        tip_dist = dist.get(tip_idx, float('inf'))
        if tip_dist < float('inf'):
            all_tip_Lnet[tip_idx] = tip_dist
            num_reachable += 1
            if tip_dist > max_Lnet:
                max_Lnet = tip_dist
                farthest_tip_idx = tip_idx
    
    # Get coordinates and calculate Euclidean distance and tortuosity
    if farthest_tip_idx >= 0:
        farthest_tip_xy = mycelia['xy2'][farthest_tip_idx].copy()
        euclidean_dist = np.sqrt(farthest_tip_xy[0]**2 + farthest_tip_xy[1]**2)
        tortuosity = max_Lnet / euclidean_dist if euclidean_dist > 0 else 1.0
    else:
        farthest_tip_xy = np.array([0.0, 0.0])
        euclidean_dist = 0.0
        tortuosity = 1.0
    
    return {
        'Lnet_max': max_Lnet,
        'farthest_tip_idx': farthest_tip_idx,
        'farthest_tip_xy': farthest_tip_xy,
        'euclidean_distance': euclidean_dist,
        'tortuosity': tortuosity,
        'all_tip_Lnet': all_tip_Lnet,
        'num_tips': len(tip_indices),
        'num_reachable_tips': num_reachable
    }


def get_Lnet_to_all_tips(mycelia, num_total_segs, origin_segments=None):
    """
    Calculate L_net (shortest path distance) from origin to ALL tips.
    
    This is a convenience wrapper that returns a sorted list of (tip_idx, Lnet) tuples.
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
    num_total_segs : int
        Current total number of segments in the mycelium.
    origin_segments : list or None, optional
        Indices of origin segments. If None, auto-detects origin.
    
    Returns
    -------
    sorted_tips : list of tuples
        List of (tip_index, Lnet_distance, xy_coordinates) sorted by Lnet (descending).
    """
    result = calculate_Lnet_to_farthest_tip(mycelia, num_total_segs, origin_segments)
    
    # Create list of (tip_idx, Lnet, xy) tuples
    tip_data = []
    for tip_idx, lnet in result['all_tip_Lnet'].items():
        xy = mycelia['xy2'][tip_idx].copy()
        tip_data.append((tip_idx, lnet, xy))
    
    # Sort by Lnet distance (descending - farthest first)
    tip_data.sort(key=lambda x: x[1], reverse=True)
    
    return tip_data


def plot_Lnet_vs_time(time_array, Lnet_array, folder_string, param_string, params, run):
    """
    Plot L_net (shortest path to farthest tip) vs time for a single run.
    
    Parameters
    ----------
    time_array : ndarray
        Array of time values (in seconds).
    Lnet_array : ndarray
        Array of L_net values corresponding to each time point.
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for file naming.
    params : dict
        Simulation parameters.
    run : int
        Run number.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Convert time units
    if params['plot_units_time'] == 'days':
        plot_times = time_array / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = time_array / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = time_array / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = time_array
        time_label = 'Time (seconds)'
    
    ax.plot(plot_times, Lnet_array, linewidth=2, marker='o', markersize=4, color='#2E7D32')
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel(f'L_net to Farthest Tip ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title('Shortest Path Distance to Farthest Tip (L_net) vs Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    final_Lnet = Lnet_array[-1]
    ax.text(0.95, 0.05, f'Final L_net: {final_Lnet:.4f}', 
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_Lnet_vs_time_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"L_net plot saved to: {fig_name}")


def plot_Lnet_vs_euclidean(time_array, Lnet_array, euclidean_array, folder_string, param_string, params, run):
    """
    Compare L_net (network distance) vs Euclidean distance over time.
    
    Parameters
    ----------
    time_array : ndarray
        Array of time values (in seconds).
    Lnet_array : ndarray
        Array of L_net values.
    euclidean_array : ndarray
        Array of Euclidean distances to farthest tip.
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for file naming.
    params : dict
        Simulation parameters.
    run : int
        Run number.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # Convert time units
    if params['plot_units_time'] == 'days':
        plot_times = time_array / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = time_array / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = time_array / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = time_array
        time_label = 'Time (seconds)'
    
    # Left plot: Both distances vs time
    ax1.plot(plot_times, Lnet_array, linewidth=2, marker='o', markersize=3, 
             color='#2E7D32', label='L_net (Network)')
    ax1.plot(plot_times, euclidean_array, linewidth=2, marker='s', markersize=3, 
             color='#1565C0', label='Euclidean')
    ax1.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax1.set_ylabel(f'Distance ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax1.set_title('L_net vs Euclidean Distance', fontsize=14, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Right plot: Tortuosity over time
    tortuosity = np.array(Lnet_array) / np.array([e if e > 0 else 1 for e in euclidean_array])
    ax2.plot(plot_times, tortuosity, linewidth=2, marker='^', markersize=3, color='#7B1FA2')
    ax2.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, label='Perfect straight line')
    ax2.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax2.set_ylabel('Tortuosity (L_net / Euclidean)', fontsize=12, fontweight='bold')
    ax2.set_title('Path Tortuosity Over Time', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')
    
    final_tort = tortuosity[-1]
    ax2.text(0.95, 0.95, f'Final Tortuosity: {final_tort:.3f}', 
            transform=ax2.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='plum', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_Lnet_comparison_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"L_net comparison plot saved to: {fig_name}")


def plot_route_factor_vs_time(time_array, route_factor_array, folder_string, param_string, params, run):
    """
    Plot Route Factor (q = L_net / L_euc) vs time for a single run.
    
    The route factor quantifies path efficiency:
    - q = 1.0: perfectly straight path (most efficient)
    - q > 1.0: winding/tortuous path (less efficient)
    
    Parameters
    ----------
    time_array : ndarray
        Array of time values (in seconds).
    route_factor_array : ndarray
        Array of route factor values (q = L_net / L_euc).
    folder_string : str
        Folder path for saving.
    param_string : str
        Parameter string for file naming.
    params : dict
        Simulation parameters.
    run : int
        Run number.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Convert time units
    if params['plot_units_time'] == 'days':
        plot_times = time_array / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = time_array / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = time_array / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = time_array
        time_label = 'Time (seconds)'
    
    ax.plot(plot_times, route_factor_array, linewidth=2, marker='o', markersize=4, color='#7B1FA2')
    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='Perfect efficiency (q=1)')
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel('Route Factor (q = L_net / L_euc)', fontsize=12, fontweight='bold')
    ax.set_title('Route Factor (Path Efficiency) vs Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    
    # Set y-axis to start at 1.0 (minimum possible value)
    y_min = min(0.95, min(route_factor_array) - 0.05)
    y_max = max(route_factor_array) * 1.1
    ax.set_ylim(y_min, y_max)
    
    final_q = route_factor_array[-1]
    avg_q = np.mean(route_factor_array)
    stats_text = f'Final q: {final_q:.3f}\nAvg q: {avg_q:.3f}'
    ax.text(0.95, 0.95, stats_text, 
            transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='plum', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_route_factor_vs_time_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Route factor plot saved to: {fig_name}")


def plot_longest_branch_vs_time(time_array, longest_branch_array, folder_string, param_string, params, run):
    """
    Plot longest branch length vs time for a single run.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Convert time units
    if params['plot_units_time'] == 'days':
        plot_times = time_array / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = time_array / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = time_array / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = time_array
        time_label = 'Time (seconds)'
    
    ax.plot(plot_times, longest_branch_array, linewidth=2, marker='o', markersize=4, color='#1976D2')
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Longest Branch Length ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title('Longest Branch Length vs Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    final_length = longest_branch_array[-1]
    ax.text(0.95, 0.05, f'Final: {final_length:.4f}', 
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_longest_branch_vs_time_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Longest branch plot saved to: {fig_name}")


def plot_longest_branch_comparison(all_times, all_longest_branch, run_labels, folder_string, param_string, params):
    """
    Compare longest branch across multiple runs.
    """
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_times)))
    
    for i, (times, longest, label) in enumerate(zip(all_times, all_longest_branch, run_labels)):
        if params['plot_units_time'] == 'days':
            plot_times = times / (60*60*24)
        elif params['plot_units_time'] == 'hours':
            plot_times = times / (60*60)
        elif params['plot_units_time'] == 'minutes':
            plot_times = times / 60
        else:
            plot_times = times
        
        ax.plot(plot_times, longest, linewidth=2, marker='o', 
                markersize=3, label=label, color=colors[i])
    
    ax.set_ylabel(f'Longest Branch Length ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title('Longest Branch Comparison Across Conditions', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    final_lengths = [longest[-1] for longest in all_longest_branch]
    max_idx = np.argmax(final_lengths)
    min_idx = np.argmin(final_lengths)
    
    stats_text = f'Longest: {run_labels[max_idx]}\nShortest: {run_labels[min_idx]}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            ha='left', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    num_runs = len(all_times)
    fig_name = f"Results/{param_string}/Avg{num_runs}/{param_string}_longest_branch_comparison_avg{num_runs}.png"
    fig.savefig(fig_name)
    plt.close()
    
    print(f"\nLongest branch comparison plot saved to: {fig_name}")


def plot_errorbar_longest_branch(count_times, avg_longest, std_longest, folder_string, param_string, params, num_runs):
    """
    Plot average longest branch with error bars across multiple runs.
    """
    # Convert units
    if params['plot_units_time'] == 'days':
        plot_times = count_times / (60*60*24)
        time_label = 'Time (days)'
    elif params['plot_units_time'] == 'hours':
        plot_times = count_times / (60*60)
        time_label = 'Time (hours)'
    elif params['plot_units_time'] == 'minutes':
        plot_times = count_times / 60
        time_label = 'Time (minutes)'
    else:
        plot_times = count_times
        time_label = 'Time (seconds)'

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.errorbar(plot_times, avg_longest, std_longest, linewidth=2, 
                marker='o', markersize=4, capsize=5, color='#1976D2')
    
    ax.set_xlabel(time_label, fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Avg. Longest Branch ({num_runs} Iterations)', fontsize=12, fontweight='bold')
    ax.set_title(f'Average Longest Branch vs Time ({num_runs} Runs)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    final_avg = avg_longest[-1]
    final_std = std_longest[-1]
    ax.text(0.95, 0.05, f'Final: {final_avg:.4f} ± {final_std:.4f}', 
            transform=ax.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()

    fig_name = "Results/{}/Avg{}/{}_avg_longest_branch_avg{}.png".format(
        param_string, num_runs, param_string, num_runs)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Average longest branch plot saved to: {fig_name}")


# ----------------------------------------------------------------------------
# NUTRIENTS AND COMBINED VISUALIZATION FUNCTIONS
# ----------------------------------------------------------------------------

def plot_biomass_vs_nutrients_consumed(biomass_array, nutrients_consumed_array, folder_string, param_string, params, run):
    """
    Plot biomass vs nutrients consumed for a single run.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    nutrients_mmol = nutrients_consumed_array * 1e12
    
    ax.plot(nutrients_mmol, biomass_array, linewidth=2, marker='o', markersize=4, color='#D32F2F')
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Total Biomass (g)', fontsize=12, fontweight='bold')
    ax.set_title('Biomass vs Nutrients Consumed', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    if nutrients_mmol[-1] > 0:
        efficiency = biomass_array[-1] / nutrients_mmol[-1]
    else:
        efficiency = 0.0
    
    stats_text = f'Final Biomass: {biomass_array[-1]:.6f} g\nEfficiency: {efficiency:.6f} g/mmol'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            ha='left', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_biomass_vs_nutrients_consumed_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Biomass vs nutrients plot saved to: {fig_name}")


def plot_biomass_vs_nutrients_comparison(all_biomass, all_nutrients_consumed, run_labels, folder_string, param_string, params):
    """
    Compare biomass vs nutrients across multiple runs.
    """
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_biomass)))
    
    for i, (biomass, nutrients, label) in enumerate(zip(all_biomass, all_nutrients_consumed, run_labels)):
        nutrients_mmol = nutrients * 1e12
        ax.plot(nutrients_mmol, biomass, linewidth=2, marker='o', 
                markersize=3, label=label, color=colors[i])
    
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Total Biomass (g)', fontsize=12, fontweight='bold')
    ax.set_title('Biomass vs Nutrients Comparison', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    sns.despine()
    plt.tight_layout()
    
    num_runs = len(all_biomass)
    fig_name = f"Results/{param_string}/Avg{num_runs}/{param_string}_biomass_vs_nutrients_comparison_avg{num_runs}.png"
    fig.savefig(fig_name)
    plt.close()
    
    print(f"\nBiomass vs nutrients comparison plot saved to: {fig_name}")


def plot_errorbar_biomass_vs_nutrients(avg_nutrients_consumed, avg_biomass, std_biomass, folder_string, param_string, params, num_runs):
    """
    Plot average biomass vs nutrients with error bars.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    nutrients_mmol = avg_nutrients_consumed * 1e12
    
    ax.errorbar(nutrients_mmol, avg_biomass, std_biomass, linewidth=2, 
                marker='o', markersize=4, capsize=5, color='#D32F2F')
    
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Avg. Total Biomass ({num_runs} Iterations) (g)', fontsize=12, fontweight='bold')
    ax.set_title(f'Average Biomass vs Nutrients ({num_runs} Runs)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    sns.despine()
    plt.tight_layout()

    fig_name = "Results/{}/Avg{}/{}_avg_biomass_vs_nutrients_avg{}.png".format(
        param_string, num_runs, param_string, num_runs)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Average biomass vs nutrients plot saved to: {fig_name}")


def plot_longest_branch_vs_nutrients_consumed(longest_branch_array, nutrients_consumed_array, folder_string, param_string, params, run):
    """
    Plot longest branch vs nutrients consumed for a single run.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    nutrients_mmol = nutrients_consumed_array * 1e12
    
    ax.plot(nutrients_mmol, longest_branch_array, linewidth=2, marker='o', markersize=4, color='#0097A7')
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Longest Branch Length ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title('Longest Branch vs Nutrients Consumed', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    sns.despine()
    plt.tight_layout()
    
    fig_name = "Results/{}/Run{}/{}_longest_branch_vs_nutrients_run{}.png".format(
        param_string, run, param_string, run)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Longest branch vs nutrients plot saved to: {fig_name}")


def plot_longest_branch_vs_nutrients_comparison(all_longest_branch, all_nutrients_consumed, run_labels, folder_string, param_string, params):
    """
    Compare longest branch vs nutrients across multiple runs.
    """
    fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_longest_branch)))
    
    for i, (longest, nutrients, label) in enumerate(zip(all_longest_branch, all_nutrients_consumed, run_labels)):
        nutrients_mmol = nutrients * 1e12
        ax.plot(nutrients_mmol, longest, linewidth=2, marker='o', 
                markersize=3, label=label, color=colors[i])
    
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Longest Branch Length ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title('Longest Branch vs Nutrients Comparison', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    sns.despine()
    plt.tight_layout()
    
    num_runs = len(all_longest_branch)
    fig_name = f"Results/{param_string}/Avg{num_runs}/{param_string}_longest_branch_vs_nutrients_comparison_avg{num_runs}.png"
    fig.savefig(fig_name)
    plt.close()
    
    print(f"\nLongest branch vs nutrients comparison plot saved to: {fig_name}")


def plot_errorbar_longest_branch_vs_nutrients(avg_nutrients_consumed, avg_longest_branch, std_longest_branch, folder_string, param_string, params, num_runs):
    """
    Plot average longest branch vs nutrients with error bars.
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    nutrients_mmol = avg_nutrients_consumed * 1e12
    
    ax.errorbar(nutrients_mmol, avg_longest_branch, std_longest_branch, linewidth=2, 
                marker='o', markersize=4, capsize=5, color='#0097A7')
    
    ax.set_xlabel('Nutrients Consumed (mmol)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Avg. Longest Branch ({num_runs} Iterations) ({params["plot_units_space"]})', fontsize=12, fontweight='bold')
    ax.set_title(f'Average Longest Branch vs Nutrients ({num_runs} Runs)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    sns.despine()
    plt.tight_layout()

    fig_name = "Results/{}/Avg{}/{}_avg_longest_branch_vs_nutrients_avg{}.png".format(
        param_string, num_runs, param_string, num_runs)
    fig.savefig(fig_name)
    plt.close()
    
    print(f"Average longest branch vs nutrients plot saved to: {fig_name}")


# ----------------------------------------------------------------------------
# MASS CONSERVATION TRACKING FUNCTIONS
# ----------------------------------------------------------------------------

def calculate_total_mass_with_metabolites(mycelia, num_total_segs, sub_e_gluc, sub_e_treha, params):
    """
    Calculate current mass state across all components.
    
    Parameters
    ----------
    mycelia : dict
        Dictionary containing mycelia structural information
    num_total_segs : int
        Current total number of segments in the mycelium
    sub_e_gluc : 2D numpy array
        External glucose grid
    sub_e_treha : 2D numpy array
        External trehalose grid
    params : dict
        Parameter dictionary containing:
        - 'yield_c': Fraction of glucose converted to cell wall (default 0.18)
        - 'convert_metabolite': Fraction of glucose converted to trehalose (default 0.10)
        - 'cross_area': Cross-sectional area of hyphae (um^2)
    
    Returns
    -------
    mass_dict : dict
        Dictionary containing:
        - external_glucose_mmol: Glucose in external grid (mmol)
        - internal_glucose_mmol: Glucose inside hyphae (mmol)
        - cell_wall_mmol: Cell wall material in hyphae (mmol)
        - internal_trehalose_mmol: Trehalose inside hyphae (mmol)
        - external_trehalose_mmol: Trehalose in external grid (mmol)
        - total_glucose_mmol: Sum of external + internal glucose (mmol)
        - total_trehalose_mmol: Sum of internal + external trehalose (mmol)
        - total_tracked_mmol: Sum of all tracked components (mmol)
        - yield_c, convert_metabolite, unmodeled_fraction: Parameters used
    """
    
    # Get yield parameters
    yield_c = params.get('yield_c', 0.18)
    convert_metabolite = params.get('convert_metabolite', 0.10)
    unmodeled_fraction = 1.0 - yield_c - convert_metabolite
    
    # Calculate external glucose (sum of all grid cells)
    external_glucose_mmol = np.sum(sub_e_gluc)
    
    # Calculate internal glucose (sum in all segments)
    gluc_i = mycelia['gluc_i'][:num_total_segs]
    internal_glucose_mmol = np.sum(gluc_i)
    
    # Calculate cell wall material in all segments
    cw_i = mycelia['cw_i'][:num_total_segs]
    cell_wall_mmol = np.sum(cw_i)
    
    # Calculate internal trehalose
    treha_i = mycelia['treha_i'][:num_total_segs]
    internal_trehalose_mmol = np.sum(treha_i)
    
    # Calculate external trehalose (sum of all grid cells)
    external_trehalose_mmol = np.sum(sub_e_treha)
    
    # Calculate totals
    total_glucose_mmol = external_glucose_mmol + internal_glucose_mmol
    total_trehalose_mmol = internal_trehalose_mmol + external_trehalose_mmol
    
    # Total tracked mass (in glucose equivalents)
    total_tracked_mmol = (total_glucose_mmol + 
                         cell_wall_mmol + 
                         total_trehalose_mmol)
    
    # Return results dictionary
    mass_dict = {
        'external_glucose_mmol': external_glucose_mmol,
        'internal_glucose_mmol': internal_glucose_mmol,
        'cell_wall_mmol': cell_wall_mmol,
        'internal_trehalose_mmol': internal_trehalose_mmol,
        'external_trehalose_mmol': external_trehalose_mmol,
        'total_glucose_mmol': total_glucose_mmol,
        'total_trehalose_mmol': total_trehalose_mmol,
        'total_tracked_mmol': total_tracked_mmol,
        'yield_c': yield_c,
        'convert_metabolite': convert_metabolite,
        'unmodeled_fraction': unmodeled_fraction
    }
    
    return mass_dict


def check_mass_conservation(initial_mass_dict, current_mass_dict):
    """
    Verify mass conservation between time points.
    
    Parameters
    ----------
    initial_mass_dict : dict
        Initial mass state (from calculate_total_mass_with_metabolites)
    current_mass_dict : dict
        Current mass state (from calculate_total_mass_with_metabolites)
    
    Returns
    -------
    conservation_dict : dict
        Dictionary containing conservation analysis:
        - initial_total_glucose: Initial total glucose (mmol)
        - current_total_glucose: Current total glucose (mmol)
        - glucose_consumed: Difference (mmol)
        - glucose_in_cell_wall: Cell wall as glucose equivalent (mmol)
        - glucose_in_trehalose: Trehalose as glucose equivalent (mmol)
        - unmodeled_metabolites: Implicitly calculated (mmol)
        - total_in_glucose_equiv: Sum of all components (mmol)
        - conservation_error_absolute: Initial - Current equiv (mmol)
        - conservation_error_percent: Error as percentage (%)
        - fraction_to_cell_wall: Fraction of consumed glucose to CW
        - fraction_to_trehalose: Fraction of consumed glucose to trehalose
        - fraction_unmodeled: Fraction of consumed glucose unmodeled
    """
    
    # Get yield parameters from current dict (same as initial)
    yield_c = current_mass_dict['yield_c']
    convert_metabolite = current_mass_dict['convert_metabolite']
    unmodeled_fraction = current_mass_dict['unmodeled_fraction']
    
    # Calculate totals in glucose equivalents
    initial_total = initial_mass_dict['total_tracked_mmol']
    
    # Current total glucose
    current_total_glucose = current_mass_dict['total_glucose_mmol']
    
    # Glucose consumed from environment
    glucose_consumed = initial_total - current_mass_dict['total_glucose_mmol']
    
    # Glucose in cell wall (stored in cw_i, need to convert back to glucose equivalent)
    # cw_i is produced such that: yield_c * glucose_consumed = cw_i
    # Therefore: glucose_used_for_cw = cw_i / yield_c
    glucose_in_cell_wall = current_mass_dict['cell_wall_mmol'] / yield_c if yield_c > 0 else 0
    
    # Glucose in trehalose
    # convert_metabolite * glucose_consumed = trehalose_produced
    # Therefore: glucose_used_for_treha = trehalose / convert_metabolite
    glucose_in_trehalose = current_mass_dict['total_trehalose_mmol'] / convert_metabolite if convert_metabolite > 0 else 0
    
    # Unmodeled metabolites (implicitly calculated from mass balance)
    unmodeled_metabolites = glucose_consumed - glucose_in_cell_wall - glucose_in_trehalose
    
    # Total in glucose equivalents
    total_in_glucose_equiv = current_total_glucose + glucose_in_cell_wall + glucose_in_trehalose + unmodeled_metabolites
    
    # Conservation error
    conservation_error_absolute = initial_total - total_in_glucose_equiv
    conservation_error_percent = (conservation_error_absolute / initial_total * 100) if initial_total > 0 else 0
    
    # Calculate fractions (only if glucose was consumed)
    if glucose_consumed > 1e-20:  # Avoid division by zero
        fraction_to_cell_wall = glucose_in_cell_wall / glucose_consumed
        fraction_to_trehalose = glucose_in_trehalose / glucose_consumed
        fraction_unmodeled = unmodeled_metabolites / glucose_consumed
    else:
        fraction_to_cell_wall = 0
        fraction_to_trehalose = 0
        fraction_unmodeled = 0
    
    conservation_dict = {
        'initial_total_glucose': initial_total,
        'current_total_glucose': current_total_glucose,
        'glucose_consumed': glucose_consumed,
        'glucose_in_cell_wall': glucose_in_cell_wall,
        'glucose_in_trehalose': glucose_in_trehalose,
        'unmodeled_metabolites': unmodeled_metabolites,
        'total_in_glucose_equiv': total_in_glucose_equiv,
        'conservation_error_absolute': conservation_error_absolute,
        'conservation_error_percent': conservation_error_percent,
        'fraction_to_cell_wall': fraction_to_cell_wall,
        'fraction_to_trehalose': fraction_to_trehalose,
        'fraction_unmodeled': fraction_unmodeled,
        'yield_c': yield_c,
        'convert_metabolite': convert_metabolite
    }
    
    return conservation_dict


def print_mass_conservation_report(conservation_dict, time_point, params):
    """
    Display formatted mass conservation report.
    
    Parameters
    ----------
    conservation_dict : dict
        Conservation analysis (from check_mass_conservation)
    time_point : float
        Current time in simulation (seconds)
    params : dict
        Parameter dictionary containing 'plot_units_time'
    
    Returns
    -------
    None (prints to console)
    """
    
    # Convert time to desired units
    if params.get('plot_units_time') == 'hours':
        display_time = time_point / 3600
        time_unit = 'hours'
    elif params.get('plot_units_time') == 'days':
        display_time = time_point / (3600 * 24)
        time_unit = 'days'
    elif params.get('plot_units_time') == 'minutes':
        display_time = time_point / 60
        time_unit = 'minutes'
    else:
        display_time = time_point
        time_unit = 'seconds'
    
    # Print report
    print("\n" + "="*70)
    print(f"MASS CONSERVATION REPORT - Time: {display_time:.2f} {time_unit}")
    print("="*70)
    
    print(f"Initial total glucose:           {conservation_dict['initial_total_glucose']:.6e} mmol")
    print()
    
    print("Current state (glucose equivalents):")
    print(f"  Glucose (free):                {conservation_dict['current_total_glucose']:.6e} mmol")
    print(f"  Cell wall (as glucose):        {conservation_dict['glucose_in_cell_wall']:.6e} mmol")
    print(f"  Trehalose (as glucose):        {conservation_dict['glucose_in_trehalose']:.6e} mmol")
    print(f"  Unmodeled metabolites:         {conservation_dict['unmodeled_metabolites']:.6e} mmol")
    print(f"  TOTAL:                         {conservation_dict['total_in_glucose_equiv']:.6e} mmol")
    print()
    
    print(f"Conservation error:              {conservation_dict['conservation_error_absolute']:.6e} mmol ({conservation_dict['conservation_error_percent']:.4f}%)")
    print()
    
    # Only show breakdown if glucose was consumed
    if conservation_dict['glucose_consumed'] > 1e-20:
        print("Glucose consumption breakdown:")
        print(f"  Total consumed:                {conservation_dict['glucose_consumed']:.6e} mmol")
        print(f"  To cell wall:                  {conservation_dict['fraction_to_cell_wall']:.4f} ({conservation_dict['fraction_to_cell_wall']*100:.2f}%)")
        print(f"  To trehalose:                  {conservation_dict['fraction_to_trehalose']:.4f} ({conservation_dict['fraction_to_trehalose']*100:.2f}%)")
        print(f"  To unmodeled:                  {conservation_dict['fraction_unmodeled']:.4f} ({conservation_dict['fraction_unmodeled']*100:.2f}%)")
    else:
        print("No glucose consumed yet (or very small amount)")
    
    print("="*70 + "\n")


def print_initial_mass_state(initial_mass_dict, params):
    """
    Display initial mass state before simulation starts.
    
    Parameters
    ----------
    initial_mass_dict : dict
        Initial mass state (from calculate_total_mass_with_metabolites)
    params : dict
        Parameter dictionary containing 'plot_units_time'
    
    Returns
    -------
    None (prints to console)
    """
    
    print("\n" + "="*70)
    print("INITIAL MASS CALCULATION")
    print("="*70)
    print(f"External glucose:                {initial_mass_dict['external_glucose_mmol']:.6e} mmol")
    print(f"Internal glucose:                {initial_mass_dict['internal_glucose_mmol']:.6e} mmol")
    print(f"Cell wall material:              {initial_mass_dict['cell_wall_mmol']:.6e} mmol")
    print(f"Internal trehalose:              {initial_mass_dict['internal_trehalose_mmol']:.6e} mmol")
    print(f"External trehalose:              {initial_mass_dict['external_trehalose_mmol']:.6e} mmol")
    print(f"Total tracked mass:              {initial_mass_dict['total_tracked_mmol']:.6e} mmol")
    print()
    print(f"Yield parameters:")
    print(f"  yield_c (to cell wall):        {initial_mass_dict['yield_c']:.4f}")
    print(f"  convert_metabolite (to treha): {initial_mass_dict['convert_metabolite']:.4f}")
    print(f"  unmodeled fraction:            {initial_mass_dict['unmodeled_fraction']:.4f}")
    print("="*70 + "\n")


# ----------------------------------------------------------------------------
# DIFFUSION CONSERVATION VERIFICATION FUNCTIONS
# ----------------------------------------------------------------------------

def check_diffusion_conservation(sub_e_before, sub_e_after, mask=None, tolerance=1e-10):
    """
    Verify that diffusion step conserved mass.
    
    Diffusion should only redistribute mass spatially, not create or destroy it.
    
    Parameters
    ----------
    sub_e_before : 2D numpy array
        External substrate concentration before diffusion (mmol per grid cell)
    sub_e_after : 2D numpy array
        External substrate concentration after diffusion (mmol per grid cell)
    mask : 2D boolean array, optional
        Mask indicating which grid cells are active (for patchy environments)
        If None, all grid cells are considered
    tolerance : float, optional
        Acceptable relative error (default 1e-10 = 0.00000001%)
        
    Returns
    -------
    is_conserved : bool
        True if mass is conserved within tolerance
    report_dict : dict
        Detailed information about the conservation check:
        - 'mass_before': Total mass before diffusion
        - 'mass_after': Total mass after diffusion
        - 'absolute_error': Absolute difference
        - 'relative_error': Relative error (fraction)
        - 'relative_error_percent': Relative error as percentage
        - 'tolerance': Tolerance used
        - 'is_conserved': Whether conservation check passed
        - 'num_active_cells': Number of grid cells checked
        
    Notes
    -----
    - This checks NUMERICAL conservation (algorithm correctness)
    - Separate from biological mass balance (uptake, metabolism, growth)
    - Small errors (~1e-12) are expected due to floating-point arithmetic
    - Larger errors suggest a bug in the diffusion implementation
    
    Examples
    --------
    >>> # Before diffusion
    >>> mass_before = np.sum(sub_e_gluc)
    >>> 
    >>> # Apply diffusion
    >>> sub_e_gluc_new = nf.finite_volume_diffusion_2d(sub_e_gluc, dt, params, mask)
    >>> 
    >>> # Check conservation
    >>> is_ok, report = hf.check_diffusion_conservation(sub_e_gluc, sub_e_gluc_new, mask)
    >>> if not is_ok:
    >>>     print(f"WARNING: Diffusion not conserved! Error: {report['relative_error_percent']:.4f}%")
    """
    
    # Apply mask if provided
    if mask is not None:
        mass_before = np.sum(sub_e_before[mask])
        mass_after = np.sum(sub_e_after[mask])
        num_active_cells = np.sum(mask)
    else:
        mass_before = np.sum(sub_e_before)
        mass_after = np.sum(sub_e_after)
        num_active_cells = sub_e_before.size
    
    # Calculate errors
    absolute_error = abs(mass_after - mass_before)
    
    if mass_before > 0:
        relative_error = absolute_error / mass_before
    else:
        # If there was no mass to begin with, check if any was created
        relative_error = 0.0 if absolute_error == 0 else float('inf')
    
    relative_error_percent = relative_error * 100.0
    
    # Check if conserved
    is_conserved = relative_error <= tolerance
    
    # Compile report
    report_dict = {
        'mass_before': mass_before,
        'mass_after': mass_after,
        'absolute_error': absolute_error,
        'relative_error': relative_error,
        'relative_error_percent': relative_error_percent,
        'tolerance': tolerance,
        'is_conserved': is_conserved,
        'num_active_cells': num_active_cells
    }
    
    return is_conserved, report_dict


def print_diffusion_conservation_warning(report_dict, substrate_name='glucose', current_time=None):
    """
    Print a formatted warning if diffusion failed to conserve mass.
    
    Parameters
    ----------
    report_dict : dict
        Report dictionary from check_diffusion_conservation()
    substrate_name : str, optional
        Name of substrate being diffused (default 'glucose')
    current_time : float, optional
        Current simulation time in seconds (for context)
    """
    
    if report_dict['is_conserved']:
        # No warning needed
        return
    
    print("Warning: diffusion mass not conserved")

    if current_time is not None:
        time_hours = current_time / 3600.0
        print(f"Time: {time_hours:.2f} hours ({current_time:.1f} seconds)")
    
    print(f"Substrate: {substrate_name}")
    print(f"Mass before diffusion: {report_dict['mass_before']:.12e} mmol")
    print(f"Mass after diffusion:  {report_dict['mass_after']:.12e} mmol")
    print(f"Absolute error:        {report_dict['absolute_error']:.12e} mmol")
    print(f"Relative error:        {report_dict['relative_error_percent']:.6f}%")
    print(f"Tolerance:             {report_dict['tolerance']*100:.6f}%")
    print(f"Active grid cells:     {report_dict['num_active_cells']}")
    
    if report_dict['relative_error_percent'] > 1.0:
        print("Error is large (>1%)")
    elif report_dict['relative_error_percent'] > 0.01:
        print("Error is moderate (>0.01%)")
    else:
        print("Error is small but above tolerance")


# ----------------------------------------------------------------------------
# BRANCH SPATIAL RATIO FUNCTIONS
# ----------------------------------------------------------------------------

def calculate_branch_spatial_ratios(mycelia, num_total_segs, origin=None):
    """
    Calculate spatial ratios (ρ_b = r_b / R_max) for all branches in the network.
    
    For any given time step in the simulation:
    1. Identify the origin (O): The starting (x, y) coordinates of the initial spore/inoculum
    2. Calculate R_max: Euclidean distance from origin to the farthest hyphal tip
    3. For each branch b: Calculate r_b (distance from origin to branch base) and ρ_b = r_b / R_max
    
    Parameters
    ----------
    mycelia : dict
        Stores structural information of mycelia colony.
        Required keys:
        - 'branch_id': segment branch IDs (-1 indicates inactive/null segment)
        - 'is_tip': boolean array indicating which segments are tips
        - 'xy1': start coordinates of each segment
        - 'xy2': endpoint coordinates of each segment
        - 'seg_id': segment ID within each branch
        - 'nbr_idxs': neighbor indices for each segment
    num_total_segs : int
        Current total number of segments in the mycelium.
    origin : tuple or ndarray, optional
        (x, y) coordinates of the origin point. If None, uses (0, 0).
        Default is None.
    
    Returns
    -------
    result : dict
        Dictionary containing:
        - 'origin': ndarray, (x, y) coordinates of the origin
        - 'R_max': float, Euclidean distance to the farthest tip
        - 'farthest_tip_idx': int, segment index of the farthest tip
        - 'farthest_tip_xy': ndarray, (x, y) coordinates of the farthest tip
        - 'branch_data': list of dicts, each containing:
            - 'branch_id': int, the branch identifier
            - 'base_xy': ndarray, (x, y) coordinates of the branch base
            - 'r_b': float, Euclidean distance from origin to branch base
            - 'rho_b': float, spatial ratio (r_b / R_max)
            - 'base_segment_idx': int, index of the first segment of the branch
            - 'parent_segment_idx': int, index of the parent segment (or -1 for initial branches)
        - 'num_branches': int, total number of active branches
    """
    # Set default origin
    if origin is None:
        origin = np.array([0.0, 0.0])
    else:
        origin = np.array(origin)
    
    # Get all active segments (branch_id > -1)
    active_segs = np.where(mycelia['branch_id'][:num_total_segs].flatten() > -1)[0]
    
    if len(active_segs) == 0:
        return {
            'origin': origin,
            'R_max': 0.0,
            'farthest_tip_idx': -1,
            'farthest_tip_xy': origin.copy(),
            'branch_data': [],
            'num_branches': 0
        }
    
    # Step 1: Calculate R_max (distance to farthest tip)
    tip_mask = mycelia['is_tip'][:num_total_segs].flatten().astype(bool)
    tip_indices = np.where(tip_mask & (mycelia['branch_id'][:num_total_segs].flatten() > -1))[0]
    
    R_max = 0.0
    farthest_tip_idx = -1
    farthest_tip_xy = origin.copy()
    
    if len(tip_indices) > 0:
        for tip_idx in tip_indices:
            tip_xy = mycelia['xy2'][tip_idx]
            dist = np.sqrt((tip_xy[0] - origin[0])**2 + (tip_xy[1] - origin[1])**2)
            if dist > R_max:
                R_max = dist
                farthest_tip_idx = tip_idx
                farthest_tip_xy = tip_xy.copy()
    
    # Step 2: Find all unique branch IDs
    branch_ids = np.unique(mycelia['branch_id'][active_segs].flatten())
    branch_ids = branch_ids[branch_ids > -1]  # Filter out any -1 values
    
    # Step 3: For each branch, find the base and calculate r_b and ρ_b
    branch_data = []
    
    for branch_id in branch_ids:
        branch_id = int(branch_id)
        
        # Find all segments belonging to this branch
        branch_segs = np.where(mycelia['branch_id'][:num_total_segs].flatten() == branch_id)[0]
        
        if len(branch_segs) == 0:
            continue
        
        # Find the first segment of this branch (seg_id == 0)
        seg_ids = mycelia['seg_id'][branch_segs].flatten()
        first_seg_idx = branch_segs[np.argmin(seg_ids)]
        
        # The base of the branch is at xy1 of the first segment
        base_xy = mycelia['xy1'][first_seg_idx].copy()
        
        # Find the parent segment (the neighbor that's not on this branch)
        parent_segment_idx = -1
        nbr_idxs = mycelia['nbr_idxs'][first_seg_idx]
        if nbr_idxs is not None:
            for nbr in nbr_idxs:
                if mycelia['branch_id'][nbr][0] != branch_id:
                    parent_segment_idx = nbr
                    break
        
        # Calculate r_b (distance from origin to branch base)
        r_b = np.sqrt((base_xy[0] - origin[0])**2 + (base_xy[1] - origin[1])**2)
        
        # Calculate ρ_b = r_b / R_max
        if R_max > 0:
            rho_b = r_b / R_max
        else:
            rho_b = 0.0
        
        branch_data.append({
            'branch_id': branch_id,
            'base_xy': base_xy,
            'r_b': r_b,
            'rho_b': rho_b,
            'base_segment_idx': first_seg_idx,
            'parent_segment_idx': parent_segment_idx
        })
    
    return {
        'origin': origin,
        'R_max': R_max,
        'farthest_tip_idx': farthest_tip_idx,
        'farthest_tip_xy': farthest_tip_xy,
        'branch_data': branch_data,
        'num_branches': len(branch_data)
    }


def save_branch_spatial_ratios_csv(spatial_ratio_result, folder_string, param_string, 
                                    current_time, run, params):
    """
    Save branch spatial ratio data to a CSV file.
    
    Parameters
    ----------
    spatial_ratio_result : dict
        Output from calculate_branch_spatial_ratios()
    folder_string : str
        Folder path for output
    param_string : str
        Parameter string for filename
    current_time : float
        Current simulation time in seconds
    run : int
        Run number
    params : dict
        Simulation parameters
    
    Returns
    -------
    output_filepath : str
        Path to the saved CSV file
    """
    import csv
    import os
    
    # Create output filename
    output_filepath = "Results/{}/Run{}/{}_branch_spatial_ratios_run{}_t={:0.2f}.csv".format(
        param_string, run, param_string, run, current_time)
    
    # Ensure directory exists
    output_dir = os.path.dirname(output_filepath)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Convert time to appropriate units for display
    if params.get('plot_units_time') == 'days':
        plot_time = current_time / (60*60*24)
        time_unit = 'days'
    elif params.get('plot_units_time') == 'hours':
        plot_time = current_time / (60*60)
        time_unit = 'hours'
    elif params.get('plot_units_time') == 'minutes':
        plot_time = current_time / 60
        time_unit = 'minutes'
    else:
        plot_time = current_time
        time_unit = 'seconds'
    
    # Write CSV file
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header metadata
        writer.writerow(['# Branch Spatial Ratio Data'])
        writer.writerow(['# Time', f'{plot_time:.4f}', time_unit])
        writer.writerow(['# Origin_X', f'{spatial_ratio_result["origin"][0]:.6f}'])
        writer.writerow(['# Origin_Y', f'{spatial_ratio_result["origin"][1]:.6f}'])
        writer.writerow(['# R_max (Maximum Radius)', f'{spatial_ratio_result["R_max"]:.6f}'])
        writer.writerow(['# Farthest_Tip_Index', f'{spatial_ratio_result["farthest_tip_idx"]}'])
        writer.writerow(['# Farthest_Tip_X', f'{spatial_ratio_result["farthest_tip_xy"][0]:.6f}'])
        writer.writerow(['# Farthest_Tip_Y', f'{spatial_ratio_result["farthest_tip_xy"][1]:.6f}'])
        writer.writerow(['# Number_of_Branches', f'{spatial_ratio_result["num_branches"]}'])
        writer.writerow([])  # Empty row for separation
        
        # Write column headers
        writer.writerow(['branch_id', 'base_x', 'base_y', 'r_b', 'rho_b', 
                        'base_segment_idx', 'parent_segment_idx'])
        
        # Write data for each branch
        for branch in spatial_ratio_result['branch_data']:
            writer.writerow([
                branch['branch_id'],
                f'{branch["base_xy"][0]:.6f}',
                f'{branch["base_xy"][1]:.6f}',
                f'{branch["r_b"]:.6f}',
                f'{branch["rho_b"]:.6f}',
                branch['base_segment_idx'],
                branch['parent_segment_idx']
            ])
    
    print(f"Branch spatial ratios saved to: {output_filepath}")

    return output_filepath


def save_Lnet_metrics_csv(Lnet_result, Leuc_result, folder_string, param_string,
                          current_time, run, params):
    """
    Save L_net, L_euc, and route factor metrics to a CSV file at each output timestep.

    Parameters
    ----------
    Lnet_result : dict
        Output from calculate_Lnet_to_farthest_tip()
    Leuc_result : dict
        Output from calculate_Leuc_to_farthest_tip()
    folder_string : str
        Folder path for output (not used, kept for compatibility)
    param_string : str
        Parameter string for filename
    current_time : float
        Current simulation time in seconds
    run : int
        Run number
    params : dict
        Simulation parameters

    Returns
    -------
    output_filepath : str
        Path to the saved CSV file
    """
    import csv
    import os

    # Create output filename
    output_filepath = "Results/{}/Run{}/{}_Lnet_metrics_run{}_t={:0.2f}.csv".format(
        param_string, run, param_string, run, current_time)

    # Ensure directory exists
    output_dir = os.path.dirname(output_filepath)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Convert time to appropriate units for display
    if params.get('plot_units_time') == 'days':
        plot_time = current_time / (60*60*24)
        time_unit = 'days'
    elif params.get('plot_units_time') == 'hours':
        plot_time = current_time / (60*60)
        time_unit = 'hours'
    elif params.get('plot_units_time') == 'minutes':
        plot_time = current_time / 60
        time_unit = 'minutes'
    else:
        plot_time = current_time
        time_unit = 'seconds'

    # Calculate route factor and tortuosity
    route_factor = Lnet_result['Lnet_max'] / Leuc_result['Leuc_max'] if Leuc_result['Leuc_max'] > 0 else 1.0

    # Write CSV file
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header metadata
        writer.writerow(['# L_net, L_euc, and Route Factor Metrics'])
        writer.writerow(['# Time', f'{plot_time:.6f}', time_unit])
        writer.writerow(['# L_euc (Euclidean distance)', f'{Leuc_result["Leuc_max"]:.6f}'])
        writer.writerow(['# L_net (Network distance)', f'{Lnet_result["Lnet_max"]:.6f}'])
        writer.writerow(['# Route Factor (q = L_net/L_euc)', f'{route_factor:.6f}'])
        writer.writerow(['# Tortuosity (L_net/L_euc)', f'{Lnet_result["tortuosity"]:.6f}'])
        writer.writerow(['# Farthest_Tip_Index', f'{Lnet_result["farthest_tip_idx"]}'])
        writer.writerow(['# Farthest_Tip_X', f'{Lnet_result["farthest_tip_xy"][0]:.6f}'])
        writer.writerow(['# Farthest_Tip_Y', f'{Lnet_result["farthest_tip_xy"][1]:.6f}'])
        writer.writerow(['# Number_of_Tips', f'{Lnet_result["num_tips"]}'])
        writer.writerow(['# Number_of_Reachable_Tips', f'{Lnet_result["num_reachable_tips"]}'])
        writer.writerow([])  # Empty row for separation

        # Write column headers
        writer.writerow(['tip_index', 'Leuc', 'Lnet', 'is_reachable'])

        # Write data for each tip
        for tip_idx in sorted(Leuc_result['all_tip_Leuc'].keys()):
            Leuc_val = Leuc_result['all_tip_Leuc'][tip_idx]
            Lnet_val = Lnet_result['all_tip_Lnet'].get(tip_idx, float('inf'))
            is_reachable = 1 if Lnet_val < float('inf') else 0
            writer.writerow([
                tip_idx,
                f'{Leuc_val:.6f}',
                f'{Lnet_val:.6f}' if Lnet_val < float('inf') else 'inf',
                is_reachable
            ])

    print(f"L_net metrics saved to: {output_filepath}")

    return output_filepath


def get_spatial_ratio_summary(spatial_ratio_result):
    """
    Get a summary of the spatial ratio distribution across branches.
    
    Parameters
    ----------
    spatial_ratio_result : dict
        Output from calculate_branch_spatial_ratios()
    
    Returns
    -------
    summary : dict
        Dictionary containing:
        - 'R_max': Maximum radius
        - 'num_branches': Number of branches
        - 'rho_min': Minimum spatial ratio
        - 'rho_max': Maximum spatial ratio
        - 'rho_mean': Mean spatial ratio
        - 'rho_std': Standard deviation of spatial ratios
        - 'rho_median': Median spatial ratio
    """
    if len(spatial_ratio_result['branch_data']) == 0:
        return {
            'R_max': 0.0,
            'num_branches': 0,
            'rho_min': 0.0,
            'rho_max': 0.0,
            'rho_mean': 0.0,
            'rho_std': 0.0,
            'rho_median': 0.0
        }
    
    rho_values = np.array([b['rho_b'] for b in spatial_ratio_result['branch_data']])
    
    return {
        'R_max': spatial_ratio_result['R_max'],
        'num_branches': spatial_ratio_result['num_branches'],
        'rho_min': np.min(rho_values),
        'rho_max': np.max(rho_values),
        'rho_mean': np.mean(rho_values),
        'rho_std': np.std(rho_values),
        'rho_median': np.median(rho_values)
    }
