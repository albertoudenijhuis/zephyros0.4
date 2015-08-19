#!/usr/bin/env python

import warnings; warnings.simplefilter("ignore")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
fontsize0 = 20
fontsize1 = 16
matplotlib.rc('xtick', labelsize=fontsize0) 
matplotlib.rc('ytick', labelsize=fontsize0) 

def main():
        
    data = []
    for rawdata in open("output.txt",'r').readlines():
        rawdata = map(float, rawdata.split())
        data.append(rawdata)   

    data = np.array(data)


    for mystyle in ['linear', 'log']:
        #make a nice plot.
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(1,1,1)

        plt2 = ax.plot(data[:,0], data[:,4] ,   label='z-dir, KH05, sheroid', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,6] ,   label='xy-dir, KH05, sheroid', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,9] ,   label='z-dir, KH05, spherical', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,11] ,   label='xy-dir, KH05, spherical', linewidth=2)

        if mystyle == 'log':
            ax.set_xscale('log')
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
        else:
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
            
        ax.set_xlabel(r'$D_{eqvol}$ [mm]')
        ax.set_ylabel(r'$\tau_I$ [m]')

        #ax.set_ybound(-2, 16.)
        
        myname = "inertialtime_"

        plt.tight_layout()
        
        plt.savefig(myname+mystyle+".png")
        plt.close(fig)


    for mystyle in ['linear', 'log']:
        #make a nice plot.
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(1,1,1)

        plt2 = ax.plot(data[:,0], data[:,3] ,   label='z-dir, KH05, sheroid', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,5] ,   label='xy-dir, KH05, sheroid', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,8] ,   label='z-dir, KH05, spherical', linewidth=2)
        plt2 = ax.plot(data[:,0], data[:,10] ,   label='xy-dir, KH05, spherical', linewidth=2)

        if mystyle == 'log':
            ax.set_xscale('log')
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
        else:
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
            
        ax.set_xlabel(r'$D_{eqvol}$ [mm]')
        ax.set_ylabel(r'$d_I$ [m]')

        #ax.set_ybound(-2, 16.)
        
        myname = "inertialdistance_"

        plt.tight_layout()
        
        plt.savefig(myname+mystyle+".png")
        plt.close(fig)
    

main();
