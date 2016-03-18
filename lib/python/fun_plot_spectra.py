import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/examples/wrapradarfilter/python")); import calc_hr_spectrum

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl


def plot_spectogram(ao_name, fname, opts = {}):
    fontsize0 = 20
    mpl.rc('xtick', labelsize=fontsize0) 
    mpl.rc('ytick', labelsize=fontsize0) 
    
    
    ao = additional_output.ao_dct(ao_name)

    recalc = {}
    recalc['spectrum_velocity_center'] = np.repeat( [np.linspace(-15.,15., 1000)], len(ao['azel_r1_m']), axis=0)
    recalc = calc_hr_spectrum.ao_to_recalc(ao, recalc)

    for plot in [
        'Doppler_spectrum_dBZ_hh',
        'Doppler_spectrum_dBZ_hv',
        'Doppler_spectrum_dBZ_vh',
        'Doppler_spectrum_dBZ_vv',
        'specific_dBZdr',
        'specific_dBLdr',
        'specific_rho_co',
        'specific_rho_cxh',
        'specific_rho_cxv',
         ]:

        if plot in recalc.keys():
            fig = plt.figure(figsize=(5,5))
            ax=plt.subplot(111)

            #plt_cmap = plt.cm.jet
            plt_cmap = plt.cm.get_cmap("jet", 100)

            plt_cmap.set_over('Black')
            plt_cmap.set_under('LightGray')
            plt_cmap.set_bad('0.')
                         
            if ('title' in opts.keys()):
                ax.set_title(opts['title'])
            
            z = recalc[plot]
            z = np.where(z == calc_hr_spectrum._FillValueminINF, np.nan, z)

            vmin = np.nanmin(z); vmax = np.nanmax(z)
            if plot in opts.keys():
                if 'vmin' in opts[plot].keys():
                    vmin = opts[plot]['vmin']
                    #z = np.where(z < vmin, np.nan, z)
                    
                if 'vmax' in opts[plot].keys():
                    vmax = opts[plot]['vmax']
                    #z = np.where(z > vmax, np.nan, z)
            
            bounds = np.linspace(vmin, vmax,101)
            #norm = mpl.colors.BoundaryNorm(bounds, plt_cmap.N)

            #~ if plot in [
            #~ 'Doppler_spectrum_dBZ_hh',
            #~ 'Doppler_spectrum_dBZ_hv',
            #~ 'Doppler_spectrum_dBZ_vh',
            #~ 'Doppler_spectrum_dBZ_vv',
                #~ ]:
                #~ for i in range(z.shape[0]):
                    #~ z[i,:] = np.where(z[i,:] < (np.nanmax(z[i,:]) - 20.), np.nan, z[i,:])
            
            nspectrum = z.shape[1]
                

            ranges = np.repeat(ao['azel_r1_m'] , nspectrum).reshape( z.shape)
            heights = np.repeat(ao['center_z'] , nspectrum).reshape( z.shape)

            CF = ax.contourf(
                recalc['spectrum_velocity_center'],
                1.e-3 * heights,
                z, 
                bounds,
                cmap=plt_cmap,
                vmin=vmin,
                vmax=vmax,
                extend='both')
                

                
            ax.set_xbound(-5., 5.)
            ax.set_xlabel("Doppler velocity [m/s]") 
            ax.set_ylabel("height [km]") 



            cb = plt.colorbar(  CF,
                                shrink=0.7,
                                #ticks=bounds,
                                spacing='uniform',
            )
            
                                                 

            #cb.set_clim(vmin-10.,vmax+10.)
            
            plt.tight_layout()


            plt.savefig(fname + '_'+plot+'.png')
            plt.close(fig)





def plot_spectrum(ao_name, fname, opts = {}):
    fontsize0 = 20
    mpl.rc('xtick', labelsize=fontsize0) 
    mpl.rc('ytick', labelsize=fontsize0) 
    
    ao = additional_output.ao_dct(ao_name)

    recalc = {}
    recalc['spectrum_velocity_center'] = np.repeat( [np.linspace(-15.,15., 500)], len(ao['azel_r1_m']), axis=0)
    d = recalc['spectrum_velocity_center'][0,1] - recalc['spectrum_velocity_center'][0,0]
    recalc['spectrum_velocity_ubound'] = recalc['spectrum_velocity_center'] + d/2.
    recalc['spectrum_velocity_lbound'] = recalc['spectrum_velocity_center'] - d/2.
    
    recalc = calc_hr_spectrum.ao_to_recalc(ao, recalc)
    calc_hr_spectrum.smooth_spectra(recalc)


    for plot in [
        'Doppler_spectrum_dBZ_hh',
        'Doppler_spectrum_dBZ_hv',
        'Doppler_spectrum_dBZ_vh',
        'Doppler_spectrum_dBZ_vv',
        'specific_dBZdr',
        'specific_dBLdr',
        'specific_rho_co',
        'specific_rho_cxh',
        'specific_rho_cxv',
         ]:

        if plot in recalc.keys():
            fig = plt.figure(figsize=(5,5))
            ax=plt.subplot(111)
                     
            if ('title' in opts.keys()):
                ax.set_title(opts['title'])
            
            z = recalc[plot][0]
            z = np.where(z == calc_hr_spectrum._FillValueminINF, np.nan, z)
            
            ax.plot(recalc['spectrum_velocity_center'][0,:], z, linewidth=2)

            ax.set_xlabel("Doppler velocity [m/s]") 
            ax.set_ylabel('[dB]')

                
            ax.set_xbound(-5., 5.)

                                                 

            
            plt.tight_layout()


            plt.savefig(fname + '_'+plot+'.png')
            plt.close(fig)






