#!/usr/bin/env python

import os, sys; sys.path.append(os.path.expanduser("~/mylib/python")); import grab
import os, sys; sys.path.append(os.path.expanduser("~/mylib/python")); import scinetcd2 as sc

import numpy as np


rawdata = {}
myfile = "stratocumulus_lwc.nc"
root    = sc.load(myfile)

rawdata['x'] = sc.readvar(root, 'x')
rawdata['y'] = sc.readvar(root, 'y')
rawdata['z'] = np.array(sc.readvar(root, 'z')) - 0.

rawdata['lwc'] = np.nan_to_num(sc.readvar(root, 'lwc'))
root.close()

x_slice = slice(None, 40, 4)
y_slice = slice(None, 40, 4)
z_slice = slice(None, None, 4)

x_vec = rawdata['x'][x_slice] - 15.e3
y_vec = rawdata['y'][y_slice] - 15.e3
z_vec = rawdata['z'][z_slice]
lwc = rawdata['lwc'][z_slice, y_slice, x_slice]

lwc = np.where(lwc == 0., 1.e-10, lwc)

n = np.product(lwc.shape)
txt_x   = "{:<15}".format(len(x_vec))
txt_y   = "{:<15}".format(len(y_vec))
txt_z   = "{:<15}".format(len(z_vec))
txt_t   = "{:<15}".format(1)
txt_N0  = "{:<15}".format(n)
txt_lwc = "{:<15}".format(n)

for i in range(len(x_vec)):
    txt_x   += "{:<15.3e}".format(x_vec[i])

for j in range(len(y_vec)):
    txt_y   += "{:<15.3e}".format(y_vec[j])

for k in range(len(z_vec)):
    txt_z   += "{:<15.3e}".format(z_vec[k])

txt_t += "{:<15.3e}".format(0.)


for k in range(len(z_vec)):
    for j in range(len(y_vec)):
        for i in range(len(x_vec)):
            txt_N0  += "{:<15.3e}".format(1.e5)
            
            #txt_lwc += "{:<15.3e}".format(np.abs(x_vec[i])/1.e5 + 1.e-10 )
            txt_lwc += "{:<15.3e}".format(lwc[k,j,i])
            
            
text = \
"""
# Hogan cloud 
# - stratocumulus
# - domain is 15 km x 15 km (xy-plane)
# - vertical domain is 0 - 1 km (z - plane)
#
#
# configuration file for Zephyros
# '# ' number sign + blank space at the beginning of a line means a comment and these lines are skipped
# comments can also be put in parentheses: (comment)
# the parsing module can not read booleans => use integers with 0 = .false. and any other number is true

section general
subsection overall
version_number       0.4    ( version number of the tool that corresponds to this configuration file)

section simulation
subsection scattererfield
psd                                 0       (counter from 0 ... max 100)
psd_distribution_type               1       (0 = discrete, 1 = gamma-distribution)
psd_particle_type                   2       (1 = spherical droplet, 2 = speroid droplet)
psd_vec_x                          """+txt_x+"""
psd_vec_y                          """+txt_y+"""
psd_vec_z                          """+txt_z+"""
psd_vec_t                          """+txt_t+"""
psd_grid_lwc_gm3                         """+txt_lwc+"""        (gridded liquid water content in g per m3 (optional). Number density will be rescaled accordingly)  
psd_grid_gammadistribution_N0            """+txt_N0+"""     (gridded number density parameter of the gamma distribution)
psd_gammadistribution_mu            5.
psd_gammadistribution_D0_mm         .1
psd_gammadistribution_dmin_mm       0.01
psd_gammadistribution_dmax_mm       8.5
psd_n_diameters                     10


"""

f = open("hogan_statocumulus_cloud.cfg", "w")
f.write(text)
f.close()
