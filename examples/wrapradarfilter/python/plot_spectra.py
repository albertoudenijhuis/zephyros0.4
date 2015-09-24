#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapradarfilter")); import radarfilter
import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import calc_hr_spectrum

import numpy as np
from copy import deepcopy
import pickle
import gzip

from scipy import stats
import fun_plot_spectra


myplotdcts = {}

myplotdcts['spectra_tara_45deg_elevation_rain_parametric_turbulence'] = \
    {'windfield': 'turbulence_parametric',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'rain',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': 'rain + parametric turbulence',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_parametric_turbulence'] = \
    {'windfield': 'turbulence_parametric',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_1mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '1 mm droplet + parametric turbulence',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_parametric_turbulence'] = \
    {'windfield': 'turbulence_parametric',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_8mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '8 mm droplet + parametric turbulence',
        }

myplotdcts['spectra_tara_45deg_elevation_rain_stochastic_turbulence'] = \
    {'windfield': 'turbulence_mann1998_small_scales',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'rain',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': 'rain + stochastic turbulence',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_stochastic_turbulence'] = \
    {'windfield': 'turbulence_mann1998_small_scales',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_1mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '1 mm droplet + stochastic turbulence',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_stochastic_turbulence'] = \
    {'windfield': 'turbulence_mann1998_small_scales',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_8mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '8 mm droplet + stochastic turbulence',
        }

myplotdcts['spectra_tara_45deg_elevation_rain_vector'] = \
    {'windfield': 'vector',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'rain',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': 'rain + vector',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_vector'] = \
    {'windfield': 'vector',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_1mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '1 mm droplet + vector',
        }
myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_vector'] = \
    {'windfield': 'vector',
     'range': 1.e3,
     'elevation_angle_deg': 45.,
     'scatterers': 'droplet_8mm',
     'radar_beam_FWHM0_rad': np.deg2rad(2.1), 'radar_beam_FWHM1_rad': np.deg2rad(2.1),
     'plotlabel': '8 mm droplet + vector',
        }


#delete some, not ready yet...
del myplotdcts['spectra_tara_45deg_elevation_rain_parametric_turbulence']
#del myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_parametric_turbulence']
del myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_parametric_turbulence']
del myplotdcts['spectra_tara_45deg_elevation_rain_stochastic_turbulence']
del myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_stochastic_turbulence']
del myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_stochastic_turbulence']
del myplotdcts['spectra_tara_45deg_elevation_rain_vector']
del myplotdcts['spectra_tara_45deg_elevation_droplet_1mm_vector']
del myplotdcts['spectra_tara_45deg_elevation_droplet_8mm_vector']

if True:
    for myplotname in myplotdcts.keys():
        myplotdct = myplotdcts[myplotname]
        print myplotname

        observations = {}
        observations['additional_output_filename']  =  "../data/"+myplotname+"_additional_output.zout"

        #if (os.path.exists(observations['additional_output_filename']) == False):
        if True:
            #general configuration files
            observations['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"
            observations['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
            
            observations['cfg_filename']                += "../../../input_files/general/instruments/tara.cfg;"
            observations['cfg_filename']                += "../../../input_files/simulation/radarfilter/all.cfg;"
            observations['cfg_filename']                += "../../../input_files/simulation/scatterers/"+myplotdct['scatterers']+".cfg;"
            observations['cfg_filename']                +=  "../../../input_files/simulation/wind/"+myplotdct['windfield']+".cfg;"

            dr = 30.
            #alpha = 0.  #hence therefore in the xz-plane.
            alpha = 90.  #hence therefore in the yz-plane.
            observations['dr']                      = dr

            observations['azel_r1_m']               = np.array([myplotdct['range']])
            observations['n']                       = len(observations['azel_r1_m'])

            observations['enu_radar_location_x']    = np.zeros(observations['n'])
            observations['enu_radar_location_y']    = np.zeros(observations['n'])
            observations['enu_radar_location_z']    = np.zeros(observations['n'])
            observations['enu_radar_time_t']        = np.zeros(observations['n'])

            observations['azel_r2_m']               = observations['azel_r1_m'] + observations['dr']
            observations['azel_alpha_rad']          = np.zeros(observations['n']) + np.deg2rad(alpha)
            observations['azel_gamma_rad']          = np.zeros(observations['n']) + np.deg2rad(myplotdct['elevation_angle_deg'])
            observations['beam_FWHM0_rad']          = np.zeros(observations['n']) + myplotdct['radar_beam_FWHM0_rad']
            observations['beam_FWHM1_rad']          = np.zeros(observations['n']) + myplotdct['radar_beam_FWHM1_rad']
            observations['dt']                      = np.ones(observations['n'])

            ao = radarfilter.calc_radar_meas(observations)
        else:
            ao = additional_output.ao_dct(observations['additional_output_filename'])

        fun_plot_spectra.plot_spectrum(observations['additional_output_filename'], 'spectra_plots/spectra_')


#~ import scipy.interpolate
#~ from matplotlib.mlab import griddata
#~ import matplotlib.pyplot as plt
#~ import matplotlib
#~ fontsize0 = 20
#~ matplotlib.rc('xtick', labelsize=fontsize0) 
#~ matplotlib.rc('ytick', labelsize=fontsize0) 
#~ 
#~ if True:
    #~ for plot in [
        #~ 'Doppler_spectrum_dBZ_hh',
        #~ 'Doppler_spectrum_dBZ_hv',
        #~ 'Doppler_spectrum_dBZ_vh',
        #~ 'Doppler_spectrum_dBZ_vv',
        #~ 'specific_dBZdr',
        #~ 'specific_dBLdr',
        #~ 'specific_rho_co',
        #~ 'specific_rho_cxh',
        #~ 'specific_rho_cxv',
         #~ ]:
#~ 
        #~ fig = plt.figure(figsize=(5,5))
        #~ ax=plt.subplot(111)
#~ 
        #~ for myplotname in sorted(myplotdcts.keys()):
            #~ ao = myplotdcts[myplotname]['ao']
    #~ 
            #~ recalc = calc_hr_spectrum.ao_to_recalc(ao)
            #~ 
            #~ if plot in recalc.keys():
                #~ x2 = recalc['x']
                #~ z2 = recalc[plot]
            #~ else:
                #~ z =  ao[plot][0,:]
                #~ z = np.where(z == 0., np.nan, z)
#~ 
                #~ x2 = np.zeros(len(z) * 2)
                #~ z2 = np.zeros(len(z) * 2)
                #~ z2[::2]     = z[::]
                #~ z2[1::2]    = z[::]
                #~ x2[::2]     = ao['spectrum_velocity_lbound'][0,:]
                #~ x2[1::2]    = ao['spectrum_velocity_ubound'][0,:]
           #~ 
            #~ 
            #~ ax.plot(x2, z2, label=myplotdcts[myplotname]['plotlabel'], linewidth=2)
#~ 
#~ 
#~ 
        #~ #if plot in ['specific_Zdr', 'specific_Ldr',]:
        #~ ax.set_ylabel('[dB]')
#~ 
#~ 
        #~ ax.set_xlabel("radial velocity [m/s]") 
        #~ #ax.set_ylabel("radar range [km]") 
        #~ 
        #~ 
        #~ leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize0)
#~ 
        #~ plt.tight_layout()
        #~ plt.savefig('spectra_plots/spectra_'+plot+'.png', bbox_extra_artists=(leg,),  bbox_inches='tight')
        #~ plt.close(fig)


if False:
    fig = plt.figure()
    ax=plt.subplot(111)

    myplotname = 'spectra_tara_90deg_elevation_rain_turbulence'
    
    ao = myplotdcts[myplotname]['ao']
    x = ao['analysis_psd0_meas0_discrete_D_equiv_mm']
    y = ao['analysis_psd0_meas0_unweighted_canting_angle_wrt_vertical_variance']
    y = (180. / np.pi) * np.sqrt(y)
    
    ax.plot(x, y) #, label='canting angle spread')

    ax.set_xlabel("equivolumetric diameter [mm]") 
    ax.set_ylabel('canting angle spread [deg]')

    ax.legend(frameon=False)

    plt.savefig('spectra_plots/test.png')
    plt.close(fig)
