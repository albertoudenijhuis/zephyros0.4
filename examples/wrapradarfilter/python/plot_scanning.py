#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapradarfilter")); import radarfilter
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
from copy import deepcopy
import pickle
import gzip

for myway in [
    'vector',
    'grid',
    'wave',
    'rankine_vortex',
    'lamb_oseen_vortex',
    'turbulence_mann1998_small_scales',
    ]:


    observations = {}
    observations['additional_output_filename']  =  "../data/scanning_"+myway+"_additional_output.zout"

    if not os.path.exists(observations['additional_output_filename']):
        #general configuration files
        observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"

        observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
        
        observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
        observations['cfg_filename']                += "../../../input_files/simulation/radarfilter/standard.cfg;"
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


    #raw_input("che??")
    import scipy.interpolate
    from matplotlib.mlab import griddata
    import matplotlib.pyplot as plt
    import matplotlib

    #from mpl_toolkits.basemap import Basemap, cm

    #from matplotlib import rc
    #rc('text',usetex=True)



    fontsize0 = 14
    matplotlib.rc('xtick', labelsize=fontsize0) 
    matplotlib.rc('ytick', labelsize=fontsize0) 


    for plot in [
        'dBZ_hh',
        'Doppler_velocity_hh_ms',
        'Doppler_spectral_width_hh_ms',
        ]:
        

        #~ else:
        fig = plt.figure(figsize=(5,5))
    
        if plot == 'dBZ_hh':
            plt_cmap = plt.cm.RdBu_r
            z = np.array(ao['dBZ_hh'])
            plt.title(r"dBZ", fontsize=fontsize0)
            vmin, vmax = z.min(), z.max()
            #vmin, vmax = -10., 10.
        if plot == 'Doppler_velocity_hh_ms':
            plt_cmap = plt.cm.RdBu_r
            z = np.array(ao['Doppler_velocity_hh_ms'])
            plt.title(r"Doppler mean velocity [m/s]", fontsize=fontsize0)
            tmp = np.max(np.abs(((1, z.min(), z.max()))))
            vmin, vmax = -tmp, tmp
            #vmin, vmax = -10., 10.

        if plot == 'Doppler_spectral_width_hh_ms':
            plt_cmap = plt.cm.rainbow
            z = np.array(np.sqrt(ao['Doppler_spectral_width_hh_ms']**2.))
            plt.title(r"spectral width [m/s]", fontsize=fontsize0)
            vmin, vmax = 0., np.max((1., z.max()))
            #vmin, vmax = 0., 5.
        
        x = np.array(ao['center_x']) * 1.e-3
        y = np.array(ao['center_y']) * 1.e-3
        xmin = -15.; xmax = 15.; ymin = -15.; ymax = 15.
        
        # Set up a regular grid of interpolation points
        xi, yi = np.linspace(xmin, xmax, 50), np.linspace(ymin, ymax, 200)
        xi, yi = np.meshgrid(xi, yi)

        plt.xlabel('x [km]', fontsize=fontsize0)
        plt.ylabel('y [km]', fontsize=fontsize0)

        zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear', fill_value=np.nan)

        ri      = np.sqrt((xi ** 2.) + (yi ** 2.))
        #angi    = np.rad2deg(np.arctan2(xi, yi))
        zi = np.where((ri < 15.), zi, np.nan)

        plt.imshow(zi, origin='lower',
                   extent=[xmin, xmax, ymin, ymax], cmap=plt_cmap, interpolation='none', vmin=vmin, vmax=vmax)
        plt.colorbar(shrink=0.6)

        plt.tight_layout()

        plt.savefig('scanning_plots/scanning_'+myway + "_"+plot+"_.png")
        plt.close(fig)
