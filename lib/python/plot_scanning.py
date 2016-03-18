#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapradarfilter")); import radarfilter
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/examples/wrapradarfilter/python")); import fun_plot_scanning


import numpy as np
from copy import deepcopy
import pickle
import gzip

for myway in [
    'vector',
    'grid',
    'wave',
    #'rankine_vortex',
    #'lamb_oseen_vortex',
    #'turbulence_mann1998_small_scales',
    ]:


    observations = {}
    observations['additional_output_filename']  =  "../data_scanning/scanning_"+myway+"_additional_output.zout"

    if not os.path.exists(observations['additional_output_filename']):
        #general configuration files
        observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"

        observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
        
        observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
        observations['cfg_filename']                += "../../../input_files/simulation/scatterers/hogan_statocumulus_cloud.cfg;"
        
        observations['cfg_filename']                +=  "../../../input_files/simulation/wind/"+myway+".cfg;"

        dr = 30.
        di = 10
        elevation_angle = 0.
        
        observations['dr']                      = dr
        observations['r1_slice']                = slice(0., 15.e3 + 1., dr * di)
        observations['phi_slice']               = slice(0., 2. * np.pi, 2. * np.pi / 72.)
        phi, r1                                 = np.mgrid[observations['phi_slice'], observations['r1_slice']]
        observations['azel_r1_m']               = deepcopy(np.ndarray.flatten(r1))
        observations['azel_alpha_rad']          = deepcopy(np.ndarray.flatten(phi))
        del r1, phi
        observations['n']                       = len(observations['azel_r1_m'])

        observations['enu_radar_location_x']    = np.zeros(observations['n'])
        observations['enu_radar_location_y']    = np.zeros(observations['n'])
        observations['enu_radar_location_z']    = np.zeros(observations['n']) + 500. 
        observations['enu_radar_time_t']        = np.zeros(observations['n'])

        observations['azel_r2_m']         = observations['azel_r1_m'] + observations['dr']
        observations['azel_gamma_rad']    = np.zeros(observations['n']) + np.deg2rad(elevation_angle)
        observations['beam_FWHM0_rad']    = np.zeros(observations['n']) + np.deg2rad(2.1)
        observations['beam_FWHM1_rad']    = np.zeros(observations['n']) + np.deg2rad(2.1)
        observations['dt']                = np.ones(observations['n']) 

        ao = radarfilter.calc_radar_meas(observations)
    else:
        ao = additional_output.ao_dct(observations['additional_output_filename'])


    fun_plot_scanning.plot_scanning(observations['additional_output_filename'], 'scanning_plots/scanning_'+myway+'_')


