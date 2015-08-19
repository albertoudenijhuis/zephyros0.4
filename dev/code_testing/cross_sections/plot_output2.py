#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib



data = []
for rawdata in open("output2.txt",'r').readlines():
    rawdata = map(float, rawdata.split())
    data.append(rawdata)   

data = np.array(data)



#nice, we obtained the data. Now we are going to make contour plots.
n = 100
wavelength = data[:,0]

names = {}
names[1] = 'De Wolf (1990), $\sigma_{hh}$'
names[2] = 'De Wolf (1990), $\sigma_{hv}$'
names[3] = 'De Wolf (1990), $\sigma_{vv}$'
names[4] = 'Mishchenko (2000), $\sigma_{hh}$'
names[5] = 'Mishchenko (2000), $\sigma_{hv}$'
names[6] = 'Mishchenko (2000), $\sigma_{vv}$'

names[14] 	= 'spherical Rayleigh backscattering'
names[16] 	= 'spherical Mie backscattering (Bohren-Huffman)'


plotnames = {}
plotnames[1] = 'hh'
plotnames[2] = 'hv'
plotnames[3] = 'vv'


for mynr in [1,2,3]:
    mynr2 = mynr + 3

    fontsize0 = 20
    matplotlib.rc('xtick', labelsize=fontsize0) 
    matplotlib.rc('ytick', labelsize=fontsize0) 
    fig = plt.figure(figsize=(8,8))
    ax=plt.subplot(111)


    ax.plot(wavelength, data[:,14], label=names[14], alpha=0.8)
    ax.plot(wavelength, data[:,16], label=names[16], alpha=0.8)
    
    ax.plot(wavelength, data[:,mynr], label=names[mynr], alpha=0.8)
    ax.plot(wavelength, data[:,mynr2], label=names[mynr2], alpha=0.8)
    
    ax.legend(frameon=False, loc='lower right')


    plt.gca().invert_xaxis()


    ax.set_xlabel('wavelength [m]', fontsize=fontsize0)
    ax.set_ylabel('cross section [dB] (S.I.)', fontsize=fontsize0)
    ax.set_xscale('log')
    

    plt.tight_layout()
    plt.savefig('output2_'+plotnames[mynr]+'.png')
    plt.close(fig)        
    plt.clf()
