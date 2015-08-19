#!/usr/bin/env python

import numpy as np
import white1999_integral
import sys

import os, sys; sys.path.append(os.path.expanduser('~/mylib/python/')); import grab

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc; rc('text',usetex=True)

def main():
	
	myfile = "white1999_integral.cfg"
	
	rawdata = grab.grablines(myfile, [12,13,14,15])

	vec_ln_a = rawdata[0][1:]
	vec_ln_b = rawdata[1][1:]
	vec_ln_L = rawdata[2][1:]
	integral_sqrt = np.array(rawdata[3][1:])
	#integral_sqrt = np.arange(len(rawdata[3][1:]))
	integral_sqrt = integral_sqrt.reshape((len(vec_ln_a), len(vec_ln_b), len(vec_ln_L)), order='F')
	
	#strip extreme values
	vec_ln_a = vec_ln_a[1:-1]
	vec_ln_b = vec_ln_b[1:-1]
	vec_ln_L = vec_ln_L[1:-1]
	integral_sqrt = integral_sqrt[1:-1,1:-1,1:-1]
	
	#for i2 in range(len(vec_ln_L)):
	#	for i1 in range(len(vec_ln_b)):
	#		for i0 in range(len(vec_ln_a)):
	#			print "integral_sqrt[{},{},{}] = {}".format(i0, i1, i2, integral_sqrt[i0,i1,i2])
	
	
	#TBD...
	for myplot in ['avsb', 'avsL', 'bvsL']:
		fig = plt.figure(figsize=(8,8))
		fontsize0 = 20
		matplotlib.rc('xtick', labelsize=fontsize0) 
		matplotlib.rc('ytick', labelsize=fontsize0) 
		ax = fig.add_subplot(1,1,1)
		plt_cmap = plt.cm.rainbow

		if myplot == 'avsb':
			tmpdata = integral_sqrt[:,:,0]
			ax.set_xlabel('a')
			ax.set_ylabel('b')
		if myplot == 'avsL':
			tmpdata = integral_sqrt[:,0,:]
			ax.set_xlabel('a')
			ax.set_ylabel('L')
		if myplot == 'bvsL':
			tmpdata = integral_sqrt[0,:,:]
			ax.set_xlabel('b')
			ax.set_ylabel('L')
			
		print tmpdata
		CS = ax.contourf(tmpdata)
		
		
		myname='test_'+myplot
		plt.savefig(myname+".png", dpi=100, bbox_inches='tight')
		plt.clf()
		
main()
