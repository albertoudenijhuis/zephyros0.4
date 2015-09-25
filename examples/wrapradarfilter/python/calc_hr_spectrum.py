#!/usr/bin/env python2.7

_FillValueminINF = -9.999e5
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm

#pr: percentile rank

def special_interpolate_via_integral(
        pr_lbound, pr_ubound, integrated_values,
        new_pr_lbound, new_pr_ubound,
        interpolation_kind = 'linear'
        ):

    #define interpolation points for the integral function
    pr = np.hstack((0., pr_lbound[0], pr_ubound))
    tmp = np.hstack( (0., 0., integrated_values) )
    tmp = np.nan_to_num(tmp)
    cdf = np.cumsum(tmp)

    if len(cdf) < 2:
        return new_pr_lbound + np.nan
    else:
        int_cdf = interp1d(pr, cdf, kind=interpolation_kind, bounds_error=False)
        
        new_integrated_values = np.array(int_cdf(new_pr_ubound) - int_cdf(new_pr_lbound))
        new_integrated_values = np.where(np.isnan(new_integrated_values), 0., new_integrated_values)
        return new_integrated_values

def dBinv(x):
    y = 10. ** (x / 10.)
    if np.isscalar(x):
        return 0. if (x == _FillValueminINF) else y
    else:
        return np.where(x == _FillValueminINF, 0., y)
def dB(x):
    y = 10. * np.log10(x)
    if np.isscalar(x):
        return _FillValueminINF if np.isnan(y) else y
    else:
        return np.where(np.isnan(y), _FillValueminINF, y)

def pnans(x):
    if np.isscalar(x):
        return np.nan if (x == _FillValueminINF) else x
    else:
        return np.where(x == _FillValueminINF, np.nan, x)
    

def update_Zdr_Ldr(ao):
    ao['specific_dBZdr'] = ao['Doppler_spectrum_dBZ_hh'] - ao['Doppler_spectrum_dBZ_vv']
    ao['specific_dBLdr'] = ao['Doppler_spectrum_dBZ_hv'] - ao['Doppler_spectrum_dBZ_vv']
    
    crit = (ao['Doppler_spectrum_dBZ_hh'] == _FillValueminINF) | (ao['Doppler_spectrum_dBZ_vv'] == _FillValueminINF)
    ao['specific_dBZdr'] = np.where(crit, np.nan, ao['specific_dBZdr'])
    
    crit = (ao['Doppler_spectrum_dBZ_hv'] == _FillValueminINF) | (ao['Doppler_spectrum_dBZ_vv'] == _FillValueminINF)    
    ao['specific_dBLdr'] = np.where(crit, np.nan, ao['specific_dBLdr'])


def ao_to_recalc(ao, recalc = None,  interpolation_kind = 'linear'):

    if recalc == None:
        n = 100
        #do a recalculation with a different resolution
        recalc = {}
        recalc['spectrum_velocity_center'] = np.zeros([ao['Doppler_spectrum_dBZ_hh'].shape[0], n]) + np.nan
        recalc['spectrum_velocity_lbound'] = np.zeros([ao['Doppler_spectrum_dBZ_hh'].shape[0], n]) + np.nan
        recalc['spectrum_velocity_ubound'] = np.zeros([ao['Doppler_spectrum_dBZ_hh'].shape[0], n]) + np.nan
        for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):
            recalc['spectrum_velocity_center'][i,:] = np.linspace(np.min(ao['spectrum_velocity_lbound'][i,0]), np.max(ao['spectrum_velocity_ubound'][i,-1]), n)

            d = recalc['spectrum_velocity_center'][i,1] - recalc['spectrum_velocity_center'][i,0]
            recalc['spectrum_velocity_lbound'][i,:] = recalc['spectrum_velocity_center'][i,:] - d/2.
            recalc['spectrum_velocity_ubound'][i,:] = recalc['spectrum_velocity_center'][i,:] + d/2.
        
    for myvar in [
        'Doppler_spectrum_dBZ_hh',
        'Doppler_spectrum_dBZ_hv',
        'Doppler_spectrum_dBZ_vh',
        'Doppler_spectrum_dBZ_vv',
        'specific_rho_co',
        'specific_rho_cxh',
        'specific_rho_cxv',
        ]:
        if myvar in ao.keys():
            recalc[myvar] = np.zeros(recalc['spectrum_velocity_center'].shape)


    for myvar in [
        'Doppler_spectrum_dBZ_hh',
        'Doppler_spectrum_dBZ_hv',
        'Doppler_spectrum_dBZ_vh',
        'Doppler_spectrum_dBZ_vv'
        ]:
        for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):
            if myvar in ao.keys():
                tmp = {}
                tmp['integrated_dBinv_'+myvar] = (
                (ao['spectrum_velocity_ubound'][i,:] - ao['spectrum_velocity_lbound'][i,:]) *
                dBinv(ao[myvar][i,:])
                )
                                
                #function to translate x, to ppt
                myf = lambda x: norm.cdf(x, ao['Doppler_velocity_hh_ms'][i], ao['Doppler_spectral_width_hh_ms'][i])             
                recalc[myvar][i,:] = dB(special_interpolate_via_integral(
                                myf(ao['spectrum_velocity_lbound'][i,:]),
                                myf(ao['spectrum_velocity_ubound'][i,:]),
                                tmp['integrated_dBinv_'+myvar],
                                myf(recalc['spectrum_velocity_lbound'][i,:]),
                                myf(recalc['spectrum_velocity_ubound'][i,:]),
                                interpolation_kind))
                                
                recalc[myvar][i,:] = recalc[myvar][i,:] * (recalc['spectrum_velocity_ubound'][i,:] - recalc['spectrum_velocity_lbound'][i,:])
                #ignore low reflection
                recalc[myvar][i,:] = np.where(recalc[myvar][i,:] > (np.nanmax(recalc[myvar][i,:]) - 25.), recalc[myvar][i,:], _FillValueminINF)

                del tmp['integrated_dBinv_'+myvar]
           
    update_Zdr_Ldr(recalc)
    
    #~ 
    for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):          
        #function to translate x, to ppt
        myf = lambda x: norm.cdf(x, ao['Doppler_velocity_hh_ms'][i], ao['Doppler_spectral_width_hh_ms'][i])             

        if 'specific_rho_co' in ao.keys():
            x = np.hstack([0., myf(ao['spectrum_velocity_center'][i,:]), 1.])
            z = np.hstack([ao['specific_rho_co'][i,0], ao['specific_rho_co'][i,:], ao['specific_rho_co'][i,-1]])
            int_ = interp1d(x, z, kind=interpolation_kind, bounds_error=False)
            recalc['specific_rho_co'][i,:] = int_(myf(recalc['spectrum_velocity_center'][i,:]))
            del int_
        if 'specific_rho_cxh' in ao.keys():
            x = np.hstack([0., myf(ao['spectrum_velocity_center'][i,:]), 1.])
            z = np.hstack([ao['specific_rho_cxh'][i,0], ao['specific_rho_cxh'][i,:], ao['specific_rho_cxh'][i,-1]])
            int_ = interp1d(x, z, kind=interpolation_kind, bounds_error=False)
            recalc['specific_rho_cxh'][i,:] = int_(myf(recalc['spectrum_velocity_center'][i,:]))
            del int_
        if 'specific_rho_cxv' in ao.keys():
            x = np.hstack([0., myf(ao['spectrum_velocity_center'][i,:]), 1.])
            z = np.hstack([ao['specific_rho_cxv'][i,0], ao['specific_rho_cxv'][i,:], ao['specific_rho_cxv'][i,-1]])
            int_ = interp1d(x, z, kind=interpolation_kind, bounds_error=False)
            recalc['specific_rho_cxv'][i,:] = int_(myf(recalc['spectrum_velocity_center'][i,:]))
            del int_
            


        
    return recalc



def smooth_spectra(ao, window_len=5):
    for myvar in [
        'Doppler_spectrum_dBZ_hh',
        'Doppler_spectrum_dBZ_hv',
        'Doppler_spectrum_dBZ_vh',
        'Doppler_spectrum_dBZ_vv'
        ]:
        for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):
            if myvar in ao.keys():
                z = ao[myvar][i,:]            
                z = pnans(z)
                w = np.ones(window_len)
                ao[myvar][i,:] = np.convolve(w/w.sum(), z,mode='same')
                
    update_Zdr_Ldr(ao)                
