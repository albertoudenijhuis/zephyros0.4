#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapradarfilter")); import radarfilter
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
from copy import deepcopy
import pickle
import gzip

from scipy import stats

from scipy.interpolate import interp1d
def special_interpolate_via_integral(
        lbound, ubound, values, newx
        ):

    #define interpolation points for the integral function
    int_x = np.hstack((lbound,ubound[-1]))
    tmp = np.hstack( (0., (ubound - lbound) * values))
    int_y = np.cumsum(tmp)

    nonans = False == (np.isnan(int_y) | np.isinf(int_y))
    int_x = np.compress(nonans, int_x)
    int_y = np.compress(nonans, int_y)

    if len(int_x) <= 1:
        return newx * np.nan
    else:           
        int_f = interp1d(int_x, int_y, kind='quadratic', bounds_error=False)
        delta = 1.e-10
        f = lambda x: (int_f(x + delta/2.) - int_f(x - delta/2.)) / delta
        
        newy = f(newx)
        return newy

dBinv = lambda x: 10. ** (x / 10.)
dB = lambda x: 10. * np.log10(x)

allobs = {}

for title in ['stochastic', 'parametric', 'stochastic_inertia_effect']:
    for instrument in ['parsax', 'tara']:
        
        myplotname = instrument + '_' + title
        allobs[myplotname] = {}
        allobs[myplotname]['title'] = title
        allobs[myplotname]['instrument'] = instrument
        
        turbulence = title
        if title == 'stochastic_inertia_effect':
            turbulence = 'stochastic'
            myfilter = 'all_inertial_effect'
        else:
            myfilter = 'all'

        allobs[myplotname]['turbulence'] = turbulence


        myplotdct = {}
        myplotdct['range_min'] = 0.
        myplotdct['range_max'] = 3.e3       #height approx 2 km
        myplotdct['beam_FWHM0_rad'] = np.deg2rad(2.1)
        myplotdct['beam_FWHM1_rad'] = np.deg2rad(2.1)

        i_skip = 1
        if instrument == 'tara':
            myplotdct['resolution'] = 30.
            i_skip = 2 #debug modus
        if instrument == 'parsax':
            myplotdct['resolution'] = 3.
            i_skip = 20 #debug modus
            

        #debug modus
        #if myfilter == 'all_inertial_effect':
        #   myplotdct['range_max'] = myplotdct['resolution'] * 1.1


        observations = {}
        observations['additional_output_filename']  =  "../data_staring/"+myplotname+"_additional_output.zout"

        if not os.path.exists(observations['additional_output_filename']):

            #general configuration files
            observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"

            observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
            
            observations['cfg_filename']                += "../../../input_files/general/instruments/"+instrument+".cfg;"
            observations['cfg_filename']                += "../../../input_files/simulation/radarfilter/"+myfilter+".cfg;"

            observations['cfg_filename']                += "../../../input_files/simulation/scatterers/profile_light_precipitation.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/simulation/wind/profile_turbulence_"+turbulence+".cfg;"

            observations['dr'] = myplotdct['resolution']
            elevation_angle = 45.
            alpha = 0.  #hence therefore in the xz-plane.

            observations['azel_r1_m']               = np.arange(myplotdct['range_min'], myplotdct['range_max'] + 1.e-10, i_skip*observations['dr'])
            observations['n']                       = len(observations['azel_r1_m'])

            observations['enu_radar_location_x']    = np.zeros(observations['n'])
            observations['enu_radar_location_y']    = np.zeros(observations['n'])
            observations['enu_radar_location_z']    = np.zeros(observations['n'])
            observations['enu_radar_time_t']        = np.zeros(observations['n'])
            
            observations['azel_r2_m']         = observations['azel_r1_m'] + observations['dr']
            observations['azel_alpha_rad']    = np.zeros(observations['n']) + np.deg2rad(alpha)
            observations['azel_gamma_rad']    = np.zeros(observations['n']) + np.deg2rad(elevation_angle)
            observations['beam_FWHM0_rad']    = np.zeros(observations['n']) + myplotdct['beam_FWHM0_rad']
            observations['beam_FWHM1_rad']    = np.zeros(observations['n']) + myplotdct['beam_FWHM1_rad']
            observations['dt']                = np.ones(observations['n'])
            
            ao = radarfilter.calc_radar_meas(observations)
        else:
            ao = additional_output.ao_dct(observations['additional_output_filename'])

        allobs[myplotname]['ao'] = ao


        if False:
            #~ #do a recalculation with nicer resolution
            recalc = {}
            recalc['shape'] = np.array(ao['spectrum_velocity_lbound'].shape)
            recalc['shape'][1] = 400
            
            recalc['x'] = np.zeros(recalc['shape'])
            recalc['Doppler_spectrum_dBZ_hh'] = np.zeros(recalc['shape'])
            recalc['Doppler_spectrum_dBZ_hv'] = np.zeros(recalc['shape'])
            recalc['Doppler_spectrum_dBZ_vh'] = np.zeros(recalc['shape'])
            recalc['Doppler_spectrum_dBZ_vv'] = np.zeros(recalc['shape'])
            recalc['specific_rho_co'] = np.zeros(recalc['shape'])
            recalc['specific_rho_cxh'] = np.zeros(recalc['shape'])
            recalc['specific_rho_cxv'] = np.zeros(recalc['shape'])
            
            for i in range(recalc['shape'][0]):
                #recalc['x'][i,:] = np.linspace(ao['spectrum_velocity_lbound'][i,0], ao['spectrum_velocity_ubound'][i,-1], 200)
                recalc['x'][i,:] = np.linspace(-10., 10., recalc['shape'][1])
                for myvar in ['Doppler_spectrum_dBZ_hh', 'Doppler_spectrum_dBZ_hv', 'Doppler_spectrum_dBZ_vh', 'Doppler_spectrum_dBZ_vv']:
                    recalc[myvar][i,:] = dB(special_interpolate_via_integral(
                                            ao['spectrum_velocity_lbound'][i,:],
                                            ao['spectrum_velocity_ubound'][i,:],
                                            dBinv(ao[myvar][i,:]),
                                            recalc['x'][i,:]))
                    #ignore low reflection
                    recalc[myvar][i,:] = np.where(recalc[myvar][i,:] > (np.nanmax(recalc[myvar][i,:]) - 10.), recalc[myvar][i,:], np.nan)
              
                recalc['specific_rho_co'][i,:] = special_interpolate_via_integral(
                                ao['spectrum_velocity_lbound'][i,:],
                                ao['spectrum_velocity_ubound'][i,:],
                                ao['specific_rho_co'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_vv'][i,:])),
                                recalc['x'][i,:]) / np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_vv'][i,:]))
                recalc['specific_rho_cxh'][i,:] = special_interpolate_via_integral(
                                ao['spectrum_velocity_lbound'][i,:],
                                ao['spectrum_velocity_ubound'][i,:],
                                ao['specific_rho_cxh'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_hv'][i,:])),
                                recalc['x'][i,:]) /  np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_hv'][i,:]))
                recalc['specific_rho_cxv'][i,:] = special_interpolate_via_integral(
                                ao['spectrum_velocity_lbound'][i,:],
                                ao['spectrum_velocity_ubound'][i,:],
                                    ao['specific_rho_cxv'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_vv'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_vh'][i,:])),
                                    recalc['x'][i,:]) /  np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_vv'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_vh'][i,:]))

            recalc['specific_dBZdr'] = recalc['Doppler_spectrum_dBZ_hh'] - recalc['Doppler_spectrum_dBZ_vv']
            recalc['specific_dBLdr'] = recalc['Doppler_spectrum_dBZ_hv'] - recalc['Doppler_spectrum_dBZ_vv']



            import scipy.interpolate
            from matplotlib.mlab import griddata
            import matplotlib.pyplot as plt
            import matplotlib


            fontsize0 = 20
            matplotlib.rc('xtick', labelsize=fontsize0) 
            matplotlib.rc('ytick', labelsize=fontsize0) 


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
            
                fig = plt.figure(figsize=(5,5))
                ax=plt.subplot(111)

                plt_cmap = plt.cm.rainbow
                z = recalc[plot]

                nspectrum = z.shape[1]

                vmin, vmax = z.min(), z.max()

                
                ranges = np.repeat(ao['azel_r1_m'] , nspectrum).reshape( z.shape)

                CF = ax.contourf(recalc['x'], 1.e-3 * ranges, z, cmap=plt_cmap)
                ax.set_xbound(-10., 10.)

                ax.set_xlabel("radial velocity [m/s]") 
                ax.set_ylabel("radar range [km]") 

                plt.colorbar(CF, shrink=0.7)

                plt.tight_layout()


                plt.savefig('staring_plots/staring_'+myplotname+'_'+plot+'.png')
                plt.close(fig)



if True:
    for instrument in ['parsax', 'tara']:

		import scipy.interpolate
		from matplotlib.mlab import griddata
		import matplotlib.pyplot as plt
		import matplotlib


		fontsize0 = 20
		fontsize1 = 10
		matplotlib.rc('xtick', labelsize=fontsize0) 
		matplotlib.rc('ytick', labelsize=fontsize0) 

		fig = plt.figure(figsize=(5,5))
		ax=plt.subplot(111)

		ax.set_xlabel("spectral width [m/s]") 
		ax.set_ylabel("height [km]") 

		for myplotname in sorted(allobs.keys()):
			if allobs[myplotname]['instrument'] == instrument:
				myao = allobs[myplotname]['ao']
				lbl = allobs[myplotname]['instrument'] + ', ' + allobs[myplotname]['title']	
			
				ax.plot(myao['Doppler_spectral_width_hh_ms'], 1.e-3 * myao['center_z'], label=lbl, linewidth=2)

		ax.legend(fontsize=fontsize1, frameon=False)

		plt.tight_layout()
		#ax.set_xbound(0., 1.)


		plt.savefig('staring_plots/staring_'+instrument+'.png')
		plt.close(fig)

