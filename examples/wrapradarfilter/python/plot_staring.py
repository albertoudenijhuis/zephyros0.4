#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapradarfilter")); import radarfilter
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
from copy import deepcopy
import pickle
import gzip

from scipy import stats

from calc_hr_spectrum import *

        
myplotdct = {}
myplotdct['range_min'] = 0.
myplotdct['range_max'] = 3.e3       #height approx 2 km
myplotdct['beam_FWHM0_rad'] = np.deg2rad(2.1)
myplotdct['beam_FWHM1_rad'] = np.deg2rad(2.1)


observations = {}
observations['additional_output_filename']  =  "../data_staring/staring_additional_output.zout"

if not os.path.exists(observations['additional_output_filename']):

	#general configuration files
	observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
	observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
	observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"

	observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
	
	observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
	observations['cfg_filename']                += "../../../input_files/simulation/radarfilter/all.cfg;"

	observations['cfg_filename']                += "../../../input_files/simulation/scatterers/profile_light_precipitation.cfg;"
	observations['cfg_filename']                +=  "../../../input_files/simulation/wind/profile_turbulence_parametric.cfg;"

	observations['dr'] = 30.
	i_skip = 1
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




if True:
	#~ #do a recalculation with nicer resolution

	import scipy.interpolate
	from matplotlib.mlab import griddata
	import matplotlib.pyplot as plt
	import matplotlib


	fontsize0 = 20
	matplotlib.rc('xtick', labelsize=fontsize0) 
	matplotlib.rc('ytick', labelsize=fontsize0) 

	recalc = ao_to_recalc(ao)


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


		CF = ax.contourf(recalc['spectrum_velocity_center'], 1.e-3 * ranges, z, cmap=plt_cmap)
		ax.set_xbound(-10., 10.)

		ax.set_xlabel("radial velocity [m/s]") 
		ax.set_ylabel("radar range [km]") 

		plt.colorbar(CF, shrink=0.7)

		plt.tight_layout()


		plt.savefig('staring_plots/staring_'+plot+'.png')
		plt.close(fig)


