#!/usr/bin/env python

import warnings; warnings.simplefilter("ignore")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
fontsize0 = 18
fontsize1 = 12
matplotlib.rc('xtick', labelsize=fontsize0) 
matplotlib.rc('ytick', labelsize=fontsize0) 

#also label fonts should be enlarged !!!

linestyles = {
	0: {'color': 'black', 'ls': '-'},
	1: {'color': 'blue', 'ls': '--'},
	2: {'color': 'green', 'ls': '-.'},
	3: {'color': 'red', 'ls': ':'},
	4: {'color': 'magenta', 'ls': '-', 'dashes':[4, 1, 1, 1, 1, 1]},
	}



def main():
        
    data = []
    for rawdata in open("output.txt",'r').readlines():
        rawdata = map(float, rawdata.split())
        data.append(rawdata)   

    data = np.array(data)


    for mystyle in ['linear', 'log']:
        #make a nice plot.
        fig = plt.figure(figsize=(6,4))
        ax = fig.add_subplot(1,1,1)


        plt2 = ax.plot(data[:,0], data[:,4] ,   label=r'$z$-dir, spheroid  ($v_t$ small)', linewidth=2, **linestyles[0])
        plt2 = ax.plot(data[:,0], data[:,6] ,   label=r'$z$-dir, spheroid ($v_t$ large)', linewidth=2, **linestyles[1])
        plt2 = ax.plot(data[:,0], data[:,8] ,   label=r'$x$/$y$-dir, spheroid', linewidth=2, **linestyles[2])
        plt2 = ax.plot(data[:,0], data[:,13] ,   label=r'$z$-dir, spherical ($v_t$ large)', linewidth=2, **linestyles[3])
        plt2 = ax.plot(data[:,0], data[:,15] ,   label=r'$x$/$y$-dir, spherical', linewidth=2, **linestyles[4])

        #plt2 = ax.plot(data[:,0], data[:,11] ,   label=r'$z$-dir, $v_t$ small, spherical', linewidth=2)


        if mystyle == 'log':
            ax.set_xscale('log')
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
        else:
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
            
        ax.set_xlabel(r'$D$ [mm]', fontsize=fontsize0)
        ax.set_ylabel(r'$\tau_I$ [s]', fontsize=fontsize0)

        #ax.set_ybound(-2, 16.)
        
        myname = "inertialtime_"

        plt.tight_layout()
        
        plt.savefig(myname+mystyle+".eps")
        plt.close(fig)


    for mystyle in ['linear', 'log']:
        #make a nice plot.
        fig = plt.figure(figsize=(6,4))
        ax = fig.add_subplot(1,1,1)

        

        plt2 = ax.plot(data[:,0], data[:,3] ,   label=r'$z$-dir, spheroid ($v_t \ll \mathcal{O}(v^{\prime}_p)$)', linewidth=2, **linestyles[0])
        plt2 = ax.plot(data[:,0], data[:,5] ,   label=r'$z$-dir, spheroid ($v_t \gg \mathcal{O}(v^{\prime}_p)$)', linewidth=2, **linestyles[1])
        plt2 = ax.plot(data[:,0], data[:,7] ,   label=r'$x$/$y$-dir, spheroid', linewidth=2, **linestyles[2])
        #plt2 = ax.plot(data[:,0], data[:,12] ,   label=r'$z$-dir, spherical ($v_t$ large)', linewidth=2, linestyle='--', color=c1)
        plt2 = ax.plot(data[:,0], data[:,12] ,   label=r'$z$-dir, spherical ($v_t \gg \mathcal{O}(v^{\prime}_p)$)', linewidth=2, **linestyles[3])
        plt2 = ax.plot(data[:,0], data[:,14] ,   label=r'$x$/$y$-dir, spherical', linewidth=2, **linestyles[4])

        #plt2 = ax.plot(data[:,0], data[:,10] ,   label=r'$z$-dir, $v_t$ small, spherical', linewidth=2)


        if mystyle == 'log':
            ax.set_xscale('log')
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
        else:
            ax.legend(loc='upper left', frameon=False, fontsize=fontsize1)
            
        ax.set_xlabel(r'$D$ [mm]', fontsize=fontsize0)
        ax.set_ylabel(r'$d_I$ [m]', fontsize=fontsize0)

        #ax.set_ybound(-2, 16.)
        
        myname = "inertialdistance_"

        plt.tight_layout()
        
        plt.savefig(myname+mystyle+".eps")
        plt.close(fig)
    

main();
