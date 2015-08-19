#!/usr/bin/env python

#written by A. Oude Nijhuis, albertoudenijhuis@gmail.com

import wrapretrieval

import os; oldcwd = os.getcwd(); os.chdir(os.path.dirname(__file__))
import sys, os; sys.path.append(os.path.expanduser("../additional_output")); import additional_output
os.chdir(oldcwd)

def retrieval(observations):
    wrapretrieval.wrapretrieval( 
        observations['cfg_filename'],
        observations['additional_output_filename'],
        observations['measurements_filename'],
    )
     
    ao = additional_output.ao_dct(observations['additional_output_filename'])
     
    return ao
