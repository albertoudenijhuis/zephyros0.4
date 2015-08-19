#!/usr/bin/env python

#written by A. Oude Nijhuis, albertoudenijhuis@gmail.com

import wrapradarfilter

import os; oldcwd = os.getcwd(); os.chdir(os.path.dirname(__file__))
import sys, os; sys.path.append(os.path.expanduser("../additional_output")); import additional_output
os.chdir(oldcwd)

def calc_radar_meas(observations):

    nrad = len(observations['azel_r1_m'])
    txtfile = ""
    
    for variable in [
        'enu_radar_location_x',
        'enu_radar_location_y',
        'enu_radar_location_z',
        'enu_radar_time_t',
        'azel_r1_m',
        'azel_r2_m',
        'azel_alpha_rad',
        'azel_gamma_rad',
        'beam_FWHM0_rad',
        'beam_FWHM1_rad',
        'dt',
        ]:
            
        if variable not in observations.keys():
            print "variable '{}' was not found in dictionary!".format(variable)
            exit(0)
            
        txtfile += "!!     {}     ".format(variable)
        txtfile += "{}     {}     ".format(1, nrad)
        for i in range(nrad): txtfile += "{:.5e}     ".format(observations[variable][i])
        txtfile += "\n"

    f = open('./tmpmeasurements.z', 'w')
    f.write(txtfile)
    f.close()
    
    wrapradarfilter.wrapradarfilter( 
        observations['cfg_filename'],
        observations['additional_output_filename'],
        './tmpmeasurements.z',
    )

    #alternative
    #mydir = os.path.dirname(__file__)
    #cmd = mydir+"/radarfilter.o '{}' '{}' '{}'".format(observations['cfg_filename'], observations['additional_output_filename'], './tmpmeasurements.z',)
    #os.system(cmd)

    #alternative with debugging
    #mydir = os.path.dirname(__file__)
    #cmd = "valgrind --leak-check=yes "+mydir+"/radarfilter.o '{}' '{}' '{}'".format(observations['cfg_filename'], observations['additional_output_filename'], './tmpmeasurements.z',)
    #print cmd
    #os.system(cmd)
    
    
    
    os.remove('./tmpmeasurements.z')
    ao = additional_output.ao_dct(observations['additional_output_filename'])
     
    return ao
