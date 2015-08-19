#!/usr/bin/env python

#written by A. Oude Nijhuis, albertoudenijhuis@gmail.com

import wrapwindfield

def windfield(o):
    no  	= o['n']
            
    (
    o['u'],
    o['v'],
    o['w'],
    ) \
     = wrapwindfield.windfield(
			o['cfg_filename'],
			o['additional_output_filename'],
            o['x'],
            o['y'],
            o['z'],
            o['t'],
            no,
            no,
            no,
            )[1:]
    
    return True
