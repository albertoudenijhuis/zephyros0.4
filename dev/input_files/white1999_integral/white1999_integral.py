#!/usr/bin/env python2.7

import numpy as np
from math import gamma    
from scipy.integrate import dblquad

def white1999_I(a, b, L):    
    def myf(phi, theta):
        return (
            12. * gamma(2./3.) * 
            (np.sin(theta) ** 3.) 
            * ((
            ((b ** 2.) * (np.cos(theta)**2.)) +
            ((a ** 2.) * (np.sin(theta)**2.)) +
            (((L**2.)/12.) * (np.sin(theta)**2.) * (np.cos(phi)**2.))
            )**(1./3.))
        )

    return dblquad(myf, 0, np.pi/2., lambda theta: 0., lambda theta: np.pi/2.)[0]
    
if __name__ == '__main__':

	tmpa = 1.e10
	tmpb = 1.e10
	tmpL = 1.e10
	tmpI = white1999_I(tmpa, tmpb, tmpL)
	
	print "tmpa = {:.2e}".format(tmpa)
	print "tmpb = {:.2e}".format(tmpb)
	print "tmpL = {:.2e}".format(tmpL)
	print "I = {:.2e}".format(tmpI)
	print "Isqrt = {:.2e}".format(np.sqrt(tmpI))
