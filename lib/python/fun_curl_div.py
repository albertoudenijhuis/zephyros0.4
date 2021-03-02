#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/additional_output")); import additional_output

import numpy as np
from copy import deepcopy

def calc_curl_and_divergence(ao_name, input_opts = {}):
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
        lst_ranges_m = np.linspace(opts['rmin_km'] * 1.e3, opts['rmax_km'] * 1.e3, opts['rangeintervals'] + 1)
        ao_range_m = (ao['azel_r1_m'] + ao['azel_r2_m']) / 2.
        
        res['range_min_km']     = lst_ranges_m[:-1] * 1.e-3
        res['range_max_km']     = lst_ranges_m[1:] * 1.e-3
        res['range_km']         = (res['range_min_km'] + res['range_max_km']) / 2.
        res['curl']    			= np.zeros(res['range_km'].shape)
        res['divergence']     	= np.zeros(res['range_km'].shape)
        
        #walk through range intervals
        for i_int in range(opts['rangeintervals']):
            range0_m = lst_ranges_m[i_int]
            range1_m = lst_ranges_m[i_int+1]
            sel_range = (range0_m < ao_range_m) & (ao_range_m < range1_m)
            sel_nonans = (np.isnan(ao['retrieved_u']) == False) & (np.isnan(ao['retrieved_v']) == False)
            
            n_ang = 100
            lst_ang = np.linspace(0, 2 * np.pi, n_ang+1, True)
            
            for i_ang in range(n_ang):
                ang0 = lst_ang[i_ang]
                ang1 = lst_ang[i_ang + 1]
                
                myindices = (
                    sel_range &
                    sel_nonans &
                    (ang0 < ao['azel_alpha_rad']) &
                    (ao['azel_alpha_rad'] < ang1)
                    )
                avg_u = np.average(ao['retrieved_u'][myindices])
                avg_v = np.average(ao['retrieved_v'][myindices])
                
                if np.isnan(avg_u) or np.isnan(avg_v):
                    myindices = (
                        sel_range &
                        sel_nonans
                        )
                    avg_u = np.average(ao['retrieved_u'][myindices])
                    avg_v = np.average(ao['retrieved_v'][myindices])
                                        
                
                if np.isnan(avg_u) or np.isnan(avg_v):
                    myindices = (
                        sel_nonans
                        )
                    avg_u = np.average(ao['retrieved_u'][myindices])
                    avg_v = np.average(ao['retrieved_v'][myindices])
                                        
                
                curl_tmp_val = (avg_u * (np.cos(ang1) - np.cos(ang0))) - (avg_v * (np.sin(ang1) - np.sin(ang0)))
                div_tmp_val = (-1. * avg_u * (np.sin(ang1) - np.sin(ang0))) + (avg_v * (np.cos(ang1) - np.cos(ang0)))
            
                res['curl'][i_int]			+= curl_tmp_val
                res['divergence'][i_int]  	+= div_tmp_val     
            

    return res
    
    
