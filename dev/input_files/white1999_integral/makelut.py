#!/usr/bin/env python

import numpy as np
import white1999_integral
import sys

def main():
    na = 50
    nb = 50
    nL = 25
    ntot = na * nb * nL
    lst_a = 10. ** np.linspace(-2., 4., na-2, endpoint=True)       #m
    lst_b = 10. ** np.linspace(-2., 4., nb-2, endpoint=True)       #m
    lst_L = 10. ** np.linspace(-1., 3., nL-2, endpoint=True)       #s

    #add some extreme value points, to avoid weird extrapolation
    lst_a = np.hstack((1.e-10, lst_a, 1.e10))
    lst_b = np.hstack((1.e-10, lst_b, 1.e10))
    lst_L = np.hstack((1.e-10, lst_L, 1.e10))

    txt0 = "{:<15}".format(na)  
    txt1 = "{:<15}".format(nb)
    txt2 = "{:<15}".format(nL)
    txt3 = "{:<15}".format(ntot)

    for a in lst_a:
        txt0 += "{:<15.3e}".format(np.log(a))
    for b in lst_b:
        txt1 += "{:<15.3e}".format(np.log(b))
    for L in lst_L:            
        txt2 += "{:<15.3e}".format(np.log(L))

    teller = -1

    for L in lst_L:            
        for b in lst_b:
            for a in lst_a:            

                print "{}/{}\n".format(teller, ntot); sys.stdout.flush()

            
                teller += 1;
                 
                #calculate I
                I = white1999_integral.white1999_I(a, b, L)
                Isqrt = np.sqrt(I)

                txt3 += "{:<15.3e}".format(Isqrt)

    text = \
"""
# configuration file for Zephyros
# '# ' number sign + blank space at the beginning of a line means a comment and these lines are skipped
# comments can also be put in parentheses: (comment)
# the parsing module can not read booleans => use integers with 0 = .false. and any other number is true

section general
subsection overall
version_number       0.4    ( version number of the tool that corresponds to this configuration file)

subsection white1999_integral
vec_ln_a            """+txt0+"""
vec_ln_b            """+txt1+"""
vec_ln_L            """+txt2+"""
integral_sqrt       """+txt3+"""
"""

    f = open("white1999_integral.cfg", "w")
    f.write(text)
    f.close()



main()
