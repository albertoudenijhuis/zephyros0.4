#!/usr/bin/env python
import re
import numpy as np

def ao_dct(file1):
    dct = {}
    for linestr in open(file1,'r').readlines():
        if linestr[:2] == '!!':
            rawdata = linestr[2:]
            rawdata = rawdata.split()
            
            #varname = linestr[2:32].strip()
            #rawdata = re.findall('(\-?[0-9]{1,10}\.?[0-9]*E?e?\+?\-?[0-9]*|nan|NaN|inf)' , linestr[2:] )
            varname = rawdata[0].strip()
            ndim    = int(rawdata[1])
            dims    = map(int,rawdata[2:2+ndim])
            rawdata2 = map(float, rawdata[2+ndim:])
            #check shape
            if len(rawdata2) == np.product(dims):
                dct[varname] = np.array(rawdata2).reshape(dims, order='C')
            else:
                print "Problem with reading out additional output. Variable ", varname, " has length ", len(rawdata2), " and dimensions ", dims
    return dct
