#!/usr/bin/env python

import os, sys; sys.path.append(os.path.expanduser("~/mylib/python")); import grab
import numpy as np

rawdata = np.array(grab.grablines("Segelstein.txt", skip=4))

n = len(rawdata)
wavelenght = 1.e-6 * rawdata[:,0]
realindex = rawdata[:,1]
imagindex = rawdata[:,2]

txt0 = "{:<15}".format(n)
txt1 = "{:<15}".format(n)
txt2 = "{:<15}".format(n)

for i in range(n):
	txt0 += "{:<15.3e}".format(wavelenght[i])
	txt1 += "{:<15.6e}".format(realindex[i])
	txt2 += "{:<15.6e}".format(imagindex[i])

text = \
"""
# configuration file for Zephyros
# '# ' number sign + blank space at the beginning of a line means a comment and these lines are skipped
# comments can also be put in parentheses: (comment)
# the parsing module can not read booleans => use integers with 0 = .false. and any other number is true

section general
subsection overall
version_number       0.4    ( version number of the tool that corresponds to this configuration file)

subsection water_refractive_index
wavelength_m   """+txt0+"""
realindex      """+txt1+"""
imagindex      """+txt2+"""
"""

f = open("water_refractive_index_segelstein.cfg", "w")
f.write(text)
f.close()
