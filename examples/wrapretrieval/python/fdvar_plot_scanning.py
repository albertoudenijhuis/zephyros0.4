#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapretrieval")); import retrieval

import os; oldcwd = os.getcwd(); os.chdir(os.path.dirname(__file__))
import sys, os; sys.path.append(os.path.expanduser("../additional_output")); import additional_output
os.chdir(oldcwd)

import numpy as np


import scipy.interpolate
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib


#for myway in [
#    'vector',
    #'grid',
    #'wave',
    #'rankine_vortex',
    #'lamb_oseen_vortex',
    #'turbulence_mann1998_small_scales',
#    ]:

    #measurements_filename = "../../wrapradarfilter/data/scanning_"+myway+"_additional_output.zout"
    #additional_output_filename = "../fdvar_data/scanning_"+myway+"_additional_output.zout"


#for development
if True:
    #myway = 'vector'
    #measurements_filename = "../../wrapradarfilter/data/scanning_vector_additional_output.zout"
    
    
    myway = 'staring'
    measurements_filename = "../../wrapradarfilter/data/staring_parametric_additional_output.zout"
    additional_output_filename = "../fdvar_data/staring_dev.zout"

    

    #if not os.path.exists(additional_output_filename):
    if True:
        observations = {}
        
        #general configuration files
        observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
        
        observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
        observations['cfg_filename']                += "../../../input_files/retrieval/radarfilter/standard.cfg;"
        observations['cfg_filename']                +=  "../../../input_files/retrieval/algorithm/windvectors_fdvar_dev.cfg;"
        observations['additional_output_filename']  =  additional_output_filename
        observations['measurements_filename']       =  measurements_filename

        ao = retrieval.retrieval(observations)
    else:
        ao = additional_output.ao_dct(additional_output_filename)





    #from mpl_toolkits.basemap import Basemap, cm

    #from matplotlib import rc
    #rc('text',usetex=True)



    fontsize0 = 16
    matplotlib.rc('xtick', labelsize=fontsize0) 
    matplotlib.rc('ytick', labelsize=fontsize0) 


    for plot in [
        'dBZ_hh',
        'Doppler_velocity_hh_ms',
        'Doppler_spectral_width_hh_ms',
        ]:
        

        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(1,1,1)
            
        if plot == 'dBZ_hh':
            plt_cmap = plt.cm.RdBu_r
            z = np.array(ao['dBZ_hh'])
            plt.title(r"dBZ")
            #vmin, vmax = z.min(), z.max()
            #vmin, vmax = -10., 10.
        if plot == 'Doppler_velocity_hh_ms':
            plt_cmap = plt.cm.RdBu_r
            z = np.array(ao['Doppler_velocity_hh_ms'])
            plt.title(r"Doppler mean velocity [m/s]")
            #vmin, vmax = z.min(), z.max()
            vmin, vmax = -10., 10.

        if plot == 'Doppler_spectral_width_hh_ms':
            plt_cmap = plt.cm.rainbow
            z = np.array(np.sqrt(ao['Doppler_spectral_width_hh_ms']**2.))
            plt.title(r"spectral width [m/s]")
            #vmin, vmax = z.min(), z.max()
            vmin, vmax = 0., 5.
        
        x = np.array(ao['center_x']) * 1.e-3
        y = np.array(ao['center_y']) * 1.e-3
        xmin = -15.; xmax = 15.; ymin = -15.; ymax = 15.
        
        # Set up a regular grid of interpolation points
        xi, yi = np.linspace(xmin, xmax, 50), np.linspace(ymin, ymax, 200)
        xi, yi = np.meshgrid(xi, yi)

        plt.xlabel('x [km]')
        plt.ylabel('y [km]')

        zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear', fill_value=np.nan)
        #zi_u = scipy.interpolate.griddata((x, y), ao['windvector_u'], (xi, yi), method='linear', fill_value=np.nan)
        #zi_v = scipy.interpolate.griddata((x, y), ao['windvector_v'], (xi, yi), method='linear', fill_value=np.nan)

        ri      = np.sqrt((xi ** 2.) + (yi ** 2.))
        #angi    = np.rad2deg(np.arctan2(xi, yi))
        zi = np.where((ri < 15.), zi, np.nan)
        #zi_u = np.where((ri < 15.), zi_u, np.nan)
        #zi_v = np.where((ri < 15.), zi_v, np.nan)



        cb = ax.imshow(zi, origin='lower',
                   extent=[xmin, xmax, ymin, ymax], cmap=plt_cmap, interpolation='none')
        plt.colorbar(cb)

        iskip = 5
        #Q = ax.quiver(
        #    np.ndarray.flatten(xi[::iskip,::iskip]),
        #    np.ndarray.flatten(yi[::iskip,::iskip]),
        #    np.ndarray.flatten(zi_u[::iskip,::iskip]),
        #    np.ndarray.flatten(zi_v[::iskip,::iskip]),)

        plt.savefig('fdvar_plots_scanning/scanning_'+myway + "_"+plot+"_.png")
        plt.close(fig)







    #~ U=griddata(ao('x'),ao('y'),observations('windvector_u'),X,Y,'cubic');
    #~ V=griddata(ao('x'),ao('y'),observations('windvector_v'),X,Y,'cubic');
    #~ 
    #~ maxval = max(abs( ao('lwm_vr')));
    #~ caxis([-maxval, maxval]);
    #~ lh = quiver(X,Y,U,V);
