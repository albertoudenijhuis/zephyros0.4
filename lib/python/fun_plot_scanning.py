#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
import scipy.interpolate
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_scanning(ao_name, fname, input_opts = {}):
    default_opts = {
        'xmin': -15.,
        'xmax': 15.,
        'ymin': -15.,
        'ymax': 15.,
        'tmin': -1.e100,
        'tmax': 1.e100,
        'dBZmin': 0.,
        'dBZmax': 50.,
        'Doppler_velocity_ms_min': -15.,
        'Doppler_velocity_ms_max': 15.,
        'file_ext': 'png',
    }
    
    opts = deepcopy(default_opts)
    opts.update(input_opts)


    fontsize0 = 12
    mpl.rc('xtick', labelsize=fontsize0) 
    mpl.rc('ytick', labelsize=fontsize0) 
    
    
    ao = additional_output.ao_dct(ao_name)
    
    okindices = np.arange(len(ao['enu_radar_time_t']))[
        (opts['tmin'] < ao['enu_radar_time_t'])
        &
        (ao['enu_radar_time_t'] < opts['tmax'])
        ]
        
    for plot in [
        'dBZ_hh',
        'Doppler_velocity_hh_ms',
        'Doppler_spectral_width_hh_ms',
        ]:
        
        if plot in ao.keys():

            fig = plt.figure(figsize=(4,4))
            ax=plt.subplot(111)

               
            #default
            #plt_cmap = plt.cm.get_cmap("viridis", 100)
            plt_cmap = plt.cm.viridis
            #~ plt_cmap.set_under('0.90')
            #~ plt_cmap.set_over('cyan')  
            
            if plot == 'dBZ_hh':
                z = np.array(ao['dBZ_hh'][okindices])
                plt.title(r"dBZ", fontsize=fontsize0)
                #vmin, vmax = z.min(), z.max()
                vmin, vmax = opts['dBZmin'], opts['dBZmax']
            if plot == 'Doppler_velocity_hh_ms':
                plt_cmap = plt.cm.RdBu_r
                z = np.array(ao['Doppler_velocity_hh_ms'][okindices])
                plt.title(r"Doppler mean velocity [m s$^{-1}$]", fontsize=fontsize0)
                tmp = np.max(np.abs(((1, z.min(), z.max()))))
                #vmin, vmax = -tmp, tmp
                vmin, vmax = opts['Doppler_velocity_ms_min'], opts['Doppler_velocity_ms_max']
            if plot == 'Doppler_spectral_width_hh_ms':
                z = np.array(np.sqrt(ao['Doppler_spectral_width_hh_ms'][okindices]**2.))
                plt.title(r"spectral width [m s$^{-1}$]", fontsize=fontsize0)
                vmin, vmax = 0., np.max((1., z.max()))
                #vmin, vmax = 0., 5.


            plt_cmap.set_over('Purple')
            plt_cmap.set_under('cyan') 
            plt_cmap.set_bad('1.')                

            #plt_cmap.set_under('0.90')



            if ('title' in opts.keys()):
                ax.set_title(opts['title'])
                        
            x = np.array(ao['center_x'][okindices]) * 1.e-3
            y = np.array(ao['center_y'][okindices]) * 1.e-3
            ang    = np.rad2deg(np.arctan2(x, y)) % 360.

            # Set up a regular grid of interpolation points
            n = 50
            dx = (opts['xmax'] - opts['xmin']) / n
            xedges = np.linspace(opts['xmin'] - dx/2., opts['xmax']+dx/2., n+1)
            dy = (opts['ymax'] - opts['ymin']) / n
            yedges = np.linspace(opts['ymin'] - dy/2., opts['ymax']+dy/2., n+1)
            xi, yi = np.linspace(opts['xmin'], opts['xmax'], n), np.linspace(opts['ymin'], opts['ymax'], n)
            xi, yi = np.meshgrid(xi, yi)

            plt.xlabel('x [km]', fontsize=fontsize0)
            plt.ylabel('y [km]', fontsize=fontsize0)

            zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear', fill_value=np.nan)
            if 'retrieved_u' in ao.keys():
                zi_u = scipy.interpolate.griddata((x, y), ao['retrieved_u'][okindices], (xi, yi), method='linear', fill_value=np.nan)
                zi_v = scipy.interpolate.griddata((x, y), ao['retrieved_v'][okindices], (xi, yi), method='linear', fill_value=np.nan)


            H, xedges, yedges = np.histogram2d(y, x, bins=(xedges, yedges))
            okdata = H > 0
            #remove badly interpolated points
            #~ ri      = np.sqrt((xi ** 2.) + (yi ** 2.))
            #~ angi    = np.rad2deg(np.arctan2(xi, yi)) % 360.
            
            #~ mydang = 360./ n
            #~ okdata  = np.zeros(ri.shape)
            #~ fct = 1.5
            #~ for myang in np.arange(0, 360., mydang):
                #~ sel = ((myang < ang) & (ang <= (myang + fct * mydang)))
                #~ if np.sum(sel) > 1:
                    #~ myrmin = np.min(1.e-3 * ao['azel_r1_m'][sel])
                    #~ myrmax = np.max(1.e-3 * ao['azel_r2_m'][sel])
                    
                    #~ sel2 = (((myang < angi) & (angi <= (myang + fct * mydang))) &
                           #~ ((myrmin < ri) &  (ri <= myrmax)))
                    #~ okdata = np.where(sel2
                                    #~ , 1., okdata)
            
				
            zi = np.where(okdata, zi, np.nan)
            
            if plot == 'dBZ_hh':				
				#plot contours
				CS = ax.contour(xi, yi, zi, [0., 10., 20., 30., 40., 50.],
					 colors='k',  # negative contours will be dashed by default
					 )
            
            if 'retrieved_u' in ao.keys():
                zi_u = np.where(okdata, zi_u, np.nan)
                zi_v = np.where(okdata, zi_v, np.nan)
                iskip = 5
                Q = ax.quiver(
                    np.ndarray.flatten(xi[::iskip,::iskip]),
                    np.ndarray.flatten(yi[::iskip,::iskip]),
                    np.ndarray.flatten(zi_u[::iskip,::iskip]),
                    np.ndarray.flatten(zi_v[::iskip,::iskip]),)
        
            plt.imshow(zi, origin='lower',
                       extent=[opts['xmin'], opts['xmax'], opts['ymin'], opts['ymax']], cmap=plt_cmap, interpolation='none', vmin=vmin, vmax=vmax)
            plt.colorbar(shrink=0.6, extend='both')



            plt.tight_layout()

            plt.savefig(fname+"_"+plot+"."+opts['file_ext'])
            plt.close(fig)




        #~ cb = ax.imshow(zi, origin='lower',
                   #~ extent=[xmin, xmax, ymin, ymax], cmap=plt_cmap, interpolation='none')
        #~ plt.colorbar(cb)

        #~ iskip = 5
        #~ #Q = ax.quiver(
        #~ #    np.ndarray.flatten(xi[::iskip,::iskip]),
        #~ #    np.ndarray.flatten(yi[::iskip,::iskip]),
        #~ #    np.ndarray.flatten(zi_u[::iskip,::iskip]),
        #~ #    np.ndarray.flatten(zi_v[::iskip,::iskip]),)

        #~ plt.savefig('fdvar_plots_scanning/scanning_'+myway + "_"+plot+"_.png")
        #~ plt.close(fig)




    #~ U=griddata(ao('x'),ao('y'),observations('windvector_u'),X,Y,'cubic');
    #~ V=griddata(ao('x'),ao('y'),observations('windvector_v'),X,Y,'cubic');
    #~ 
    #~ maxval = max(abs( ao('lwm_vr')));
    #~ caxis([-maxval, maxval]);
    #~ lh = quiver(X,Y,U,V);


