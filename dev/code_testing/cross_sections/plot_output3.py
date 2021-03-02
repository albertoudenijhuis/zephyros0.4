#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

mywavel = {}
mywavel['a'] = '3.298 GHz, 9.09 cm (TARA)'
mywavel['b'] = '9.487 GHz (IDRA), 3.16 cm (IDRA)'
mywavel['c'] = '35 GHz, 8.56 mm (cloud radar)'

myname = {}
myname['a'] = 'tara'
myname['b'] = 'idra'
myname['c'] = 'cloudradar'

for sel_abc in ['a', 'b', 'c']:
    
    data = []
    for rawdata in open("output3_"+sel_abc+".txt",'r').readlines():
        rawdata = map(float, rawdata.split())
        data.append(rawdata)   

    data = np.array(data)



    #nice, we obtained the data. Now we are going to make contour plots.
    diameter_mm = data[:,0]



    names = {}
    names[2]    = 'spherical Rayleigh backscattering'
    names[6]    = 'De Wolf (1990), spherical, $\sigma_{hh}$'
    names[12]   = 'Mishchenko (2000), spherical, $\sigma_{hh}$'
    names[18]   = 'De Wolf (1990), spheroid, $\sigma_{hh}$'
    names[19]   = 'De Wolf (1990), spheroid, $\sigma_{hv}$'
    names[20]   = 'De Wolf (1990), spheroid, $\sigma_{vv}$'
    names[24]   = 'Mishchenko (2000), spheroid, $\sigma_{hh}$'
    names[25]   = 'Mishchenko (2000), spheroid, $\sigma_{hv}$'
    names[26]   = 'Mishchenko (2000), spheroid, $\sigma_{vv}$'
    names[29]   = 'spherical Mie backscattering (Bohren-Huffman)'



    for plot in [
        'spherical', 
        'spheroid_hh', 
        'spheroid_hv',
        'spheroid_vv', 
        'spheroid_Zdr', 
        #'spheroid_Ldr', 
        ]:
        
        fontsize0 = 20
        matplotlib.rc('xtick', labelsize=fontsize0) 
        matplotlib.rc('ytick', labelsize=fontsize0) 
        fig = plt.figure(figsize=(8,8))
        ax=plt.subplot(111)

        mydct = {}
        for mynr in [2, 29, 6, 12, 18, 24, 19, 25,29,20,26]:
            mytmpname = names[mynr]
            mydct[mytmpname]   = data[:,mynr]

            #~ mydct[myr]['spheroid_Zdr']   = mydct[myr]['spheroid_hh'] / mydct[myr]['spheroid_vv']
            #~ mydct[myr]['spheroid_Ldr']   = mydct[myr]['spheroid_hv'] / mydct[myr]['spheroid_vv']
        
        mydct['De Wolf (1990), spheroid, dBZdr'] = 10. * np.log10( (
            mydct['De Wolf (1990), spheroid, $\sigma_{hh}$']/
            mydct['De Wolf (1990), spheroid, $\sigma_{vv}$']) )
            
        mydct['Mishchenko (2000), spheroid, dBZdr'] = 10. * np.log10( (
            mydct['Mishchenko (2000), spheroid, $\sigma_{hh}$'] / 
            mydct['Mishchenko (2000), spheroid, $\sigma_{vv}$'] ))

        
        
        if plot == 'spherical':
            nrs = [2,29,6,12]
        if plot == 'spheroid_hh':
            nrs = [2,29, 18,24]
        if plot == 'spheroid_hv':
            nrs = [2,19,25]
        if plot == 'spheroid_vv':
            nrs = [2,29, 20,26]
                    
        if plot == 'spheroid_Zdr':
            for mykey in [              
                'De Wolf (1990), spheroid, dBZdr',
                'Mishchenko (2000), spheroid, dBZdr',
                ]:
                ax.plot(diameter_mm, mydct[mykey], label=mykey, alpha=0.6)
                    

                       
        else:
            for mynr in nrs:
                mytmpname = names[mynr]
                ax.plot(diameter_mm, mydct[mytmpname], label=names[mynr], alpha=0.6)

        ax.legend(frameon=False, loc='upper left')


        ax.set_xlabel('equivolumetric diameter [mm]', fontsize=fontsize0)
        ax.set_ylabel('cross section [dB] (S.I.)', fontsize=fontsize0)
        #ax.set_xscale('log')

        ax.set_title(mywavel[sel_abc])
        
        plt.tight_layout()
        
        plt.savefig('output3_'+myname[sel_abc]+'_'+plot+'.png')
        plt.close(fig)        
        plt.clf()
