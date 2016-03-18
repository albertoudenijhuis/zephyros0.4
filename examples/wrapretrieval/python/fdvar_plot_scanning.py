#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapretrieval")); import retrieval
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/examples/wrapradarfilter/python")); import fun_plot_scanning


import os; oldcwd = os.getcwd(); os.chdir(os.path.dirname(__file__))
import sys, os; sys.path.append(os.path.expanduser("../additional_output")); import additional_output
os.chdir(oldcwd)

import numpy as np


for myway in [
    'wave',
    'vector',
    'grid',																											
    #'rankine_vortex',
    #'lamb_oseen_vortex',
    #'turbulence_mann1998_small_scales',
    ]:

    measurements_filename = "../../wrapradarfilter/data_scanning/scanning_"+myway+"_additional_output.zout"
    additional_output_filename = "../fdvar_data/scanning_"+myway+"_additional_output.zout"


    #if not os.path.exists(additional_output_filename):
    if True:
        observations = {}
        
        #general configuration files
        observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
        
        observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/retrieval/algorithm/windvectors_fdvar_horizontal_hdir_solution.cfg;"
        observations['additional_output_filename']  =  additional_output_filename
        observations['measurements_filename']       =  measurements_filename

        ao = retrieval.retrieval(observations)
    else:
        ao = additional_output.ao_dct(additional_output_filename)

    opts = {
        'Doppler_velocity_ms_min': -10.,
        'Doppler_velocity_ms_max': 10.,
    }

    fun_plot_scanning.plot_scanning(additional_output_filename, 'fdvar_plots_scanning/scanning_'+myway+'_', opts)




