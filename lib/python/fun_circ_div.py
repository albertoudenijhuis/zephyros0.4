#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
from copy import deepcopy

def calc_circulation_divergence(ao_name, input_opts = {}):
    default_opts = {
        'rangeintervals': 20,
        'rmin_km': 0.,
        'rmax_km': 15.,
    }
    opts = deepcopy(default_opts)
    opts.update(input_opts)
    
    ao = additional_output.ao_dct(ao_name)
    
    res = {}
    
    
    if ('retrieved_u' in ao.keys()) and ('retrieved_v' in ao.keys()):
        lst_ranges = np.linspace(opts['rmin_km'] * 1.e3, opts['rmax_km'] * 1.e3, opts['rangeintervals'] + 1)
        ao_range = (ao['azel_r1_m'] + ao['azel_r2_m']) / 2.
        
        res['range_min_km']             = lst_ranges[:-1] * 1.e-3
        res['range_max_km']             = lst_ranges[1:] * 1.e-3
        res['range_km']                 = (res['range_min_km'] + res['range_max_km']) / 2.
        res['circulation_inproduct']    = np.zeros(res['range_km'].shape)
        res['divergence_inproduct']     = np.zeros(res['range_km'].shape)
        
        #walk through range intervals
        for i_int in range(opts['rangeintervals']):
            range0 = lst_ranges[i_int]
            range1 = lst_ranges[i_int+1]
            sel = (range0 < ao_range) & (ao_range < range1)
            
            #obtain results
            thisdata = {}
            thisdata['range']           = np.compress(sel, ao_range)
            thisdata['x']               = np.compress(sel, ao['center_x'])
            thisdata['y']               = np.compress(sel, ao['center_y'])
            thisdata['azel_alpha_rad']  = np.compress(sel, ao['azel_alpha_rad'])
            thisdata['u']               = np.compress(sel, ao['retrieved_u'])
            thisdata['v']               = np.compress(sel, ao['retrieved_v'])
            
            thisdata['norm_x']  = thisdata['x'] / np.sqrt((thisdata['x']**2.) + (thisdata['y']**2.))
            thisdata['norm_y']  = thisdata['y'] / np.sqrt((thisdata['x']**2.) + (thisdata['y']**2.))
            
            #calculate inproducts           
            thisdata['circulation_inproduct']   = (thisdata['u'] * thisdata['norm_x']) + (thisdata['v'] * thisdata['norm_y'])
            thisdata['divergence_inproduct']    = (thisdata['u'] * thisdata['norm_y']) + (thisdata['v'] * thisdata['norm_x'])

            #delete nans
            thisdata['circulation_inproduct'] = np.compress(np.isnan(thisdata['circulation_inproduct']) == False, thisdata['circulation_inproduct'])
            thisdata['divergence_inproduct'] = np.compress(np.isnan(thisdata['divergence_inproduct']) == False, thisdata['divergence_inproduct'])

            #calculate weighting function
            n_ang = 100
            lst_ang = np.linspace(0, 2 * np.pi, n_ang+1)
            thisdata['weights']         = np.zeros(thisdata['range'].shape)
            for i_ang in range(n_ang):
                myindicies = np.argmin(
                    (lst_ang[i_ang] < thisdata['azel_alpha_rad']) &
                    (thisdata['azel_alpha_rad'] < lst_ang[i_ang + 1])
                    )
                thisdata['weights'][myindicies] = np.sum(myindicies)
            thisdata['weights'] = thisdata['weights'] / np.sum(thisdata['weights'])

            res['circulation_inproduct'][i_int] = np.sum(thisdata['weights'] * thisdata['circulation_inproduct'])
            res['divergence_inproduct'][i_int] = np.sum(thisdata['weights'] * thisdata['divergence_inproduct'])

            del thisdata

    return res
    
    
