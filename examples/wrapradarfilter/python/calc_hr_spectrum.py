#!/usr/bin/env python2.7

_FillValueminINF = -9.999e5
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm

def special_interpolate_via_integral(
        lbound, ubound, values, newx, interpolation_kind = 'cubic'
        ):

    #define interpolation points for the integral function
    int_x = np.hstack((lbound,ubound[-1]))
    tmp = np.hstack( (0., (ubound - lbound) * values))
    tmp = np.nan_to_num(tmp)
    int_y = np.cumsum(tmp)

    nonans = False == (np.isnan(int_y) | np.isinf(int_y))
    int_x = np.compress(nonans, int_x)
    int_y = np.compress(nonans, int_y)

    if len(int_y) < 2:
        return newx + np.nan
    else:
        int_f = interp1d(int_x, int_y, kind=interpolation_kind, bounds_error=False)
        delta = 1.e-10
        f = lambda x: (int_f(x + delta/2.) - int_f(x - delta/2.)) / delta
        
        newy = np.array(f(newx))
        newy = np.where(np.isnan(newy), 0., newy)   
        return newy

def dBinv(x):
    y = 10. ** (x / 10.)
    if np.isscalar(x):
        return 0. if (y == _FillValueminINF) else y
    else:
        return np.where(y == _FillValueminINF, 0., y)
def dB(x):
    y = 10. * np.log10(x)
    if np.isscalar(x):
        return _FillValueminINF if np.isnan(y) else y
    else:
        return np.where(np.isnan(y), _FillValueminINF, y)

def ao_to_recalc(ao, recalc = None,  interpolation_kind = 'cubic'):


    if recalc == None:
        n = 100
        #do a recalculation with a different resolution
        recalc = {}
        recalc['spectrum_velocity_center'] = np.zeros([ao['Doppler_spectrum_dBZ_hh'].shape[0], n]) + np.nan
        for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):
            recalc['spectrum_velocity_center'][i,:] = np.linspace(np.min(ao['spectrum_velocity_lbound'][i,0]), np.max(ao['spectrum_velocity_ubound'][i,-1]), n)

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
                #function to translate x, to ppt
                myf = lambda x: norm.cdf(x, ao['Doppler_velocity_hh_ms'][i], ao['Doppler_spectral_width_hh_ms'][i])             
                recalc[myvar][i,:] = dB(special_interpolate_via_integral(
                                myf(ao['spectrum_velocity_lbound'][i,:]),
                                myf(ao['spectrum_velocity_ubound'][i,:]),
                                dBinv(ao[myvar][i,:]),
                                myf(recalc['spectrum_velocity_center'][i,:]), interpolation_kind))
                #ignore low reflection
                recalc[myvar][i,:] = np.where(recalc[myvar][i,:] > (np.nanmax(recalc[myvar][i,:]) - 25.), recalc[myvar][i,:], _FillValueminINF)

    recalc['specific_dBZdr'] = recalc['Doppler_spectrum_dBZ_hh'] - recalc['Doppler_spectrum_dBZ_vv']
    recalc['specific_dBLdr'] = recalc['Doppler_spectrum_dBZ_hv'] - recalc['Doppler_spectrum_dBZ_vv']
    
    crit = (recalc['Doppler_spectrum_dBZ_hh'] == _FillValueminINF) | (recalc['Doppler_spectrum_dBZ_vv'] == _FillValueminINF)
    recalc['specific_dBZdr'] = np.where(crit, np.nan, recalc['specific_dBZdr'])
    
    crit = (recalc['Doppler_spectrum_dBZ_hv'] == _FillValueminINF) | (recalc['Doppler_spectrum_dBZ_vv'] == _FillValueminINF)    
    recalc['specific_dBLdr'] = np.where(crit, np.nan, recalc['specific_dBLdr'])
    
    #~ 
    #~ for i in range(ao['Doppler_spectrum_dBZ_hh'].shape[0]):          
        #~ recalc['specific_rho_co'][i,:] = special_interpolate_via_integral(
                        #~ ao['spectrum_velocity_lbound'][i,:],
                        #~ ao['spectrum_velocity_ubound'][i,:],
                        #~ ao['specific_rho_co'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_vv'][i,:])),
                        #~ recalc['spectrum_velocity_center'][i,:], interpolation_kind) / np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_vv'][i,:]))
        #~ recalc['specific_rho_cxh'][i,:] = special_interpolate_via_integral(
                        #~ ao['spectrum_velocity_lbound'][i,:],
                        #~ ao['spectrum_velocity_ubound'][i,:],
                        #~ ao['specific_rho_cxh'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_hv'][i,:])),
                        #~ recalc['spectrum_velocity_center'][i,:], interpolation_kind) /  np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_hh'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_hv'][i,:]))
        #~ recalc['specific_rho_cxv'][i,:] = special_interpolate_via_integral(
                        #~ ao['spectrum_velocity_lbound'][i,:],
                        #~ ao['spectrum_velocity_ubound'][i,:],
                        #~ ao['specific_rho_cxv'][i,:] * np.sqrt(dBinv(ao['Doppler_spectrum_dBZ_vv'][i,:]) * dBinv(ao['Doppler_spectrum_dBZ_vh'][i,:])),
                        #~ recalc['spectrum_velocity_center'][i,:], interpolation_kind) /  np.sqrt(dBinv(recalc['Doppler_spectrum_dBZ_vv'][i,:]) * dBinv(recalc['Doppler_spectrum_dBZ_vh'][i,:]))
        #~ 
    #~ 
    return recalc
