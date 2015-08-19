#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib



data = []
for rawdata in open("output.txt",'r').readlines():
    rawdata = map(float, rawdata.split())
    data.append(rawdata)   

data = np.array(data)



#nice, we obtained the data. Now we are going to make contour plots.
n = 100
AZ = np.rad2deg(np.reshape(data[:,0], (n,n)))
EL = np.rad2deg(np.reshape(data[:,1], (n,n)))

#update data
mymax = np.max(np.absolute(np.hstack(data[:,2:])))
data[:,2:] = np.where(np.abs(data[:,2:]) <= (mymax * 1.e-20), np.nan, data[:,2:])

names = {}
names[2] = 'De Wolf (1990), $\sigma_{hh}$'
names[3] = 'De Wolf (1990), $\sigma_{hv}$'
names[4] = 'De Wolf (1990), $\sigma_{vv}$'
names[5] = 'Mishchenko (2000), $\sigma_{hh}$'
names[6] = 'Mishchenko (2000), $\sigma_{hv}$'
names[7] = 'Mishchenko (2000), $\sigma_{vv}$'

plotnames = {}
plotnames[2] = 'dewolf_hh'
plotnames[3] = 'dewolf_hv'
plotnames[4] = 'dewolf_vv'
plotnames[5] = 'mishchenko_hh'
plotnames[6] = 'mishchenko_hv'
plotnames[7] = 'mishchenko_vv'


mymin_hh_vv = np.nanmin(10. * np.log10(np.absolute(np.hstack(data[:,[2,4,5,7]]))))
mymax_hh_vv = np.nanmax(10. * np.log10(np.absolute(np.hstack(data[:,[2,4,5,7]]))))

mymin_hv = np.nanmin(10. * np.log10(np.absolute(np.hstack(data[:,[3,6]]))))
mymax_hv = np.nanmax(10. * np.log10(np.absolute(np.hstack(data[:,[3,6]]))))

plot_max = {}
plot_max[2] = mymax_hh_vv
plot_max[3] = mymax_hv
plot_max[4] = mymax_hh_vv
plot_max[5] = mymax_hh_vv
plot_max[6] = mymax_hv
plot_max[7] = mymax_hh_vv

plot_min = {}
for i in range(2,8):
	plot_min[i] = plot_max[i] - 20.


#plot_min = {}
#plot_min[2] = mymin_hh_vv
#plot_min[3] = mymin_hv
#plot_min[4] = mymin_hh_vv
#plot_min[5] = mymin_hh_vv
#plot_min[6] = mymin_hv
#plot_min[7] = mymin_hh_vv



for mynr in [2,3,4,5,6,7]:
    fontsize0 = 20
    matplotlib.rc('xtick', labelsize=fontsize0) 
    matplotlib.rc('ytick', labelsize=fontsize0) 
    fig = plt.figure(figsize=(8,8))
    ax=plt.subplot(111)

    #Z = np.log(np.reshape(data[:,mynr], (n, n)))
    Z = np.reshape(data[:,mynr], (n, n))
    ZdB = 10. * np.log10(Z)
    
    v = np.linspace(plot_min[mynr], plot_max[mynr], 100, endpoint=True)
    CS = ax.contourf(AZ, EL, ZdB, v, cmap=plt.cm.rainbow, vmin=plot_min[mynr], vmax=plot_max[mynr])
    #CS = ax.contourf(AZ, EL, ZdB, 100, cmap=plt.cm.rainbow)
    CB = plt.colorbar(CS)
    CB.set_label('cross section [dB] S.I.', fontsize=fontsize0)


    ax.set_xlabel('radar azimuth angle [deg]', fontsize=fontsize0)
    ax.set_ylabel('radar elevation angle [deg]', fontsize=fontsize0)

    plt.title(names[mynr])
    plt.tight_layout()
    plt.savefig('output_'+plotnames[mynr]+'.png')
    plt.close(fig)        
    plt.clf()
