#!/usr/bin/env python

#written by A. Oude Nijhuis, albertoudenijhuis@gmail.com

#import wrapretrieval

import os; oldcwd = os.getcwd(); os.chdir(os.path.dirname(__file__))
import sys, os; sys.path.append(os.path.expanduser("../lib/python")); import additional_output as ao
os.chdir(oldcwd)

def retrieval(observations):
    #~ wrapretrieval.wrapretrieval( 
        #~ observations['cfg_filename'],
        #~ observations['additional_output_filename'],
        #~ observations['measurements_filename'],
    #~ )
     
    #alternative
    mydir = os.path.dirname(__file__)
    cmd = mydir+"/wrapretrieval.o '{}' '{}' '{}'".format(observations['cfg_filename'], observations['additional_output_filename'], observations['measurements_filename'])
    print cmd
    os.system(cmd)

    #alternative with debugging
    #~ mydir = os.path.dirname(__file__)
    #~ #cmd = "valgrind --leak-check=yes "+mydir+"/wrapretrieval.o '{}' '{}' '{}'".format(observations['cfg_filename'], observations['additional_output_filename'], observations['measurements_filename'])
    #~ cmd = "valgrind "+mydir+"/wrapretrieval.o '{}' '{}' '{}'".format(observations['cfg_filename'], observations['additional_output_filename'], observations['measurements_filename'])
    #~ print cmd
    #~ os.system(cmd)
         
    ao_dct  = ao.ao_dct(observations['additional_output_filename'])
     
    return ao_dct
