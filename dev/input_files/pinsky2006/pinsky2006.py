#!/usr/bin/env python

__author        = 'Albert Oude Nijhuis <albertoudenijhuis@gmail.com>'
__institute     = 'Delft University of Technology'
__date          = '2014, October 29'
__version       = '0.1'
__project       = 'EU FP7 program, the UFO project.'
__dissemination = 'Confidential, only for members of the UFO project. Potentially public in the future.'
__description   = \
"""
Repetition of generation of homogeneous isotropic turbulence, pinsky 2006.
Time integration is not applied.

"""

import numpy as np

from pprint import pprint 
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import matplotlib

from scipy.special import jn_zeros, jn
from scipy import integrate
          
from copy import deepcopy


def pinsky2006_R(dct, r):
    #dct['eps']
    #dct['nu']
    
    C = 2.0
    G = dct['eps'] * ( (30. * dct['nu']) ** -1.)
    F = ((15. * C * dct['nu']) ** (-3. / 2.)) * (dct['eps'] ** (1./2.))
    Y = 1. + F * (r**2.)
        
    R = (
            (2. * dct['E']) - 
            (   ((4./3.) * G / F) 
                * 
                (   (2. * (Y ** (1./3.)))
                    -
                    (Y ** (2./3.))
                    - 
                    (Y**(5./3.))
                )
            )
        )
    return R
    
#input
#dct['r0'], dct['K'], dct['M']
def pinsky2006_muk_lambdak(dct):
    dct['R'] = lambda r: pinsky2006_R(dct, r)

    #calculations
    dct['r0'] = 40. #m, value in the article
    
    dct['muk'] = jn_zeros(0, dct['K'])
    dct['lambdak'] = dct['muk'] / dct['r0']
    
    #obtain coefficients alphak
    dct['alphaksq'] = np.zeros(dct['K'])
    for k in range(dct['K']):
        dct['alphaksq'][k] = (2. / dct['M']) * ((dct['r0'] / (dct['muk'][k] * jn(1, dct['muk'][k]))) ** 2.)
        myfun = lambda s: s * dct['R'](dct['r0'] * s) * jn(0, dct['muk'][k] * s)
        tmp = integrate.quad(myfun, 0., 1.)[0]
        dct['alphaksq'][k] *= tmp
    dct['alphak'] = np.where(dct['alphaksq'] <= 0., 0., np.sqrt(dct['alphaksq']))



def obtain_velocity(dct):
    pinsky2006_muk_lambdak(dct)

    #grid
    dct['NN']       = dct['N'] * dct['N']
    dct['x'] = range(dct['N'])
    dct['z'] = range(dct['N'])
    dct['xx'], dct['zz'] = np.meshgrid(dct['x'], dct['z'])  
    
    #a and b are random numbers with zero mean, and variance of 1.
    dct['a'] = np.random.normal(loc=0, scale=1,  size=dct['K']*dct['M']).reshape( (dct['K'],dct['M']))
    dct['b'] = np.random.normal(loc=0, scale=1,  size=dct['K']*dct['M']).reshape( (dct['K'],dct['M']))
    
    
    dct['Vx'] = np.zeros((dct['N'],dct['N']))
    dct['Vz'] = np.zeros((dct['N'],dct['N']))
    
    
    phim = (2. * np.pi / dct['M']) * np.arange(1, dct['M'] + 0.1)
    for k in range(dct['K']):
        print "{}/{}".format(k, dct['K'])
        for i in range(dct['N']):
            for j in range(dct['N']):

                tmparr = (  (dct['a'][k,:] * np.cos(dct['lambdak'][k] * ((dct['xx'][i,j] * np.cos(phim)) + (dct['zz'][i,j] * np.sin(phim))))) +
                            (dct['b'][k,:] * np.sin(dct['lambdak'][k] * ((dct['xx'][i,j] * np.cos(phim)) + (dct['zz'][i,j] * np.sin(phim)))))
                            )

                dct['Vx'][i,j] += dct['alphak'][k] * np.sum(dct['lambdak'][k] * np.sin(phim) * tmparr)
                dct['Vz'][i,j] += dct['alphak'][k] * np.sum(dct['lambdak'][k] * np.cos(phim) * tmparr)
                
                
def repeat_article_cases():


    for plot in [
        'pinsky2006_test'
        ]:

        dct = {}
        dct['eps']  = 1.e-4     #typical value of EDR, S.I. units
        dct['nu']   = 1.e-5     #typical value of kinematic viscosity, S.I. units
        dct['E']    = 1.        #mean square value of each flow velocity component
        dct['N']    = 32
        dct['K']    = 150
        dct['M']    = 180

        if plot == 'pinsky2006_test':
            dct['R'] = lambda r: pinsky2006_R(dct, r)

        obtain_velocity(dct)


        
        #plot   
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}
        matplotlib.rc('font', **font)

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1,1,1)
        ax.grid(True)
        skip_i = 1
        Q = ax.quiver(dct['xx'][::skip_i,::skip_i],dct['zz'][::skip_i,::skip_i],dct['Vx'][::skip_i,::skip_i],dct['Vz'][::skip_i,::skip_i],
            )#pivot='mid', scale=1,headwidth=5)

        ax.set_xbound(0, dct['N'])
        ax.set_ybound(0, dct['N'])
        
        plt.tight_layout()

        myname=plot+'.png'
        plt.savefig(myname)

        plt.close(fig)

if __name__ == '__main__':
    repeat_article_cases()

    mydct = {'eps': 1.e-2, 'nu': 1.e-5, 'E': 1., 'r0': 1.e2, 'K': 100, 'M': 100}
    pinsky2006_muk_lambdak(mydct)
    f = open("./pinsky2006_M{}_K{}.z".format(mydct['M'], mydct['K']),'w')
    for i in range(mydct['K']):
        f.write("{:.20e}\n".format(mydct['lambdak'][i]))
    for i in range(mydct['K']):
        f.write("{:.20e}\n".format(mydct['muk'][i]))
    f.close()
