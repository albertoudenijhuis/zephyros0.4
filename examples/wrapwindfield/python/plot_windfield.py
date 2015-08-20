#!/usr/bin/env python2.7

import sys, os; sys.path.append(os.path.expanduser("~/tools/zephyros0.4/wrapwindfield")); import windfield

import numpy as np
from copy import deepcopy

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text',usetex=True)

from matplotlib.mlab import griddata

myplotdcts = {}
#myplotdcts['vector'] 							= {'windfield': 'vector', 'plane':'xz', 'scale':10.e3, 'unit':'m'}
#myplotdcts['grid'] 								= {'windfield': 'grid', 'plane':'xz', 'scale':10.e3, 'unit':'m'}
#myplotdcts['wave'] 								= {'windfield': 'wave', 'plane':'xz', 'scale':10.e3, 'unit':'m'}
#myplotdcts['rankine_vortex'] 					= {'windfield': 'rankine_vortex', 'plane':'xy', 'scale':10.e3, 'unit':'m'}
#myplotdcts['lamb_oseen_vortex'] 				= {'windfield': 'lamb_oseen_vortex', 'plane':'xy', 'scale':10.e3, 'unit':'m'}
myplotdcts['turbulence_mann1998_small_scales'] 	= {'windfield': 'turbulence_mann1998_small_scales', 'plane':'xz', 'scale':40., 'unit':'m'}
#myplotdcts['turbulence_pinsky2006_xyplane'] 	= {'windfield': 'turbulence_pinsky2006_xyplane', 'plane':'xy', 'scale':100., 'unit':'m'}
#myplotdcts['profile_turbulence_stochastic'] 	= {'windfield': 'profile_turbulence_stochastic', 'plane':'xz', 'scale':50., 'unit':'m'}


for myplotname in myplotdcts.keys():
	myplotdct = myplotdcts[myplotname]
	print myplotname

	wind = {}
	wind['cfg_filename']                =  "../../../input_files/general/standard_output.cfg;"
	wind['cfg_filename']                +=  "../../../input_files/general/water_refractive_index_segelstein.cfg;"
	wind['cfg_filename']                +=  "../../../input_files/general/white1999_integral.cfg;"

	wind['cfg_filename']                +=  "../../../input_files/general/atmosphere/US1976.cfg;"
	
	wind['cfg_filename']                += "../../../input_files/general/instruments/idra.cfg;"
	wind['cfg_filename']                += "../../../input_files/simulation/radarfilter/all.cfg;"
	wind['cfg_filename']                += "../../../input_files/simulation/scatterers/rain.cfg;"
	wind['cfg_filename']                +=  "../../../input_files/simulation/wind/"+myplotdct['windfield']+".cfg;"
	wind['additional_output_filename']  =  "../data/"+myplotname+"_additional_output.zout"

	n = 15

	if myplotdct['plane'] == 'xy':
		xmin = -.5 * myplotdct['scale']; xmax = .5 * myplotdct['scale'] + 1.e-10
		ymin = -.5 * myplotdct['scale']; ymax = .5 * myplotdct['scale'] + 1.e-10
		wind['x_slice'] = slice(xmin, xmax, myplotdct['scale'] / n)
		wind['y_slice'] = slice(ymin, ymax, myplotdct['scale'] / n)
		wind['z_slice'] = slice(0., 0.1, 1.)     
		wind['t_slice'] = slice(0., 0.1, 1.)     #0
	if myplotdct['plane'] == 'xz':
		xmin = -.5 * myplotdct['scale']; xmax = .5 * myplotdct['scale'] + 1.e-10
		zmin = 0.; zmax = 1. * myplotdct['scale'] + 1.e-10
		wind['x_slice'] = slice(xmin, xmax, myplotdct['scale'] / n)
		wind['y_slice'] = slice(0., 0.1, 1.)     
		wind['z_slice'] = slice(zmin, zmax, myplotdct['scale'] / n)
		wind['t_slice'] = slice(0., 0.1, 1.)     
	if myplotdct['plane'] == 'yz':
		ymin = -.5 * myplotdct['scale']; ymax = .5 * myplotdct['scale'] + 1.e-10
		zmin = 0.; zmax = 1. * myplotdct['scale'] + 1.e-10
		wind['x_slice'] = slice(0., 0.1, 1.)     
		wind['y_slice'] = slice(ymin, ymax, myplotdct['scale'] / n)
		wind['z_slice'] = slice(zmin, zmax, myplotdct['scale'] / n)
		wind['t_slice'] = slice(0., 0.1, 1.)     

	x, y, z, t = np.mgrid[wind['x_slice'], wind['y_slice'], wind['z_slice'], wind['t_slice']]
	wind['x'] = deepcopy(np.ndarray.flatten(x))     
	wind['y'] = deepcopy(np.ndarray.flatten(y))     
	wind['z'] = deepcopy(np.ndarray.flatten(z))     
	wind['t'] = deepcopy(np.ndarray.flatten(t))     
	wind['n'] = len(wind['x'])
	del x,y,z, t

	windfield.windfield(wind)
			
			
	fontsize0 = 14
	matplotlib.rc('xtick', labelsize=fontsize0) 
	matplotlib.rc('ytick', labelsize=fontsize0) 

	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(1,1,1)
	if myplotdct['unit'] == 'km':
		sc = 1.e-3
	if myplotdct['unit'] == 'm':
		sc = 1.

	mydelta = 1.e-25
	wind['u'] = np.where(np.abs(wind['u']) < mydelta, 0., wind['u'])
	wind['v'] = np.where(np.abs(wind['v']) < mydelta, 0., wind['v'])
	wind['w'] = np.where(np.abs(wind['w']) < mydelta, 0., wind['w'])

	
	vel_abs = np.sqrt(wind['u']**2. + wind['v']**2. + wind['w']**2.)
	
	if myplotdct['plane'] == 'xy':
		Q = ax.quiver(
			sc * wind['x'],
			sc * wind['y'],
			wind['u'],
			wind['v'],)
		ax.set_xlabel('x ['+myplotdct['unit']+']', fontsize=fontsize0)
		ax.set_ylabel('y ['+myplotdct['unit']+']', fontsize=fontsize0)
		extralab = 'z=0'
		ax.set_xbound(sc * xmin, sc * xmax)
		ax.set_ybound(sc * ymin, sc * ymax)
		
	if myplotdct['plane'] == 'xz':
		Q = ax.quiver(
			sc * wind['x'],
			sc * wind['z'],
			wind['u'],
			wind['w'],)
		ax.set_xlabel('x ['+myplotdct['unit']+']', fontsize=fontsize0)
		ax.set_ylabel('z ['+myplotdct['unit']+']', fontsize=fontsize0)
		extralab = '@y=0'
		ax.set_xbound(sc * xmin, sc * xmax)
		ax.set_ybound(sc * zmin, sc * zmax)
		
	if myplotdct['plane'] == 'yz':
		Q = ax.quiver(
			sc * wind['y'],
			sc * wind['z'],
			wind['v'],
			wind['w'],)
		ax.set_xlabel('y ['+myplotdct['unit']+']', fontsize=fontsize0)
		ax.set_ylabel('z ['+myplotdct['unit']+']', fontsize=fontsize0)
		extralab = 'x=0'
		ax.set_xbound(sc * ymin, sc * ymax)
		ax.set_ybound(sc * zmin, sc * zmax)	  
	  
	qk = ax.quiverkey(Q, 0.9, 1.05, 1, r'1 ms$^{-1}$',
						labelpos='W',
						fontproperties={'weight': 'bold', 'size':fontsize0})
		
	ax.set_title('wind field '+extralab, fontsize=fontsize0)

	plt.tight_layout()
	
	plt.savefig('plots/windfield_'+myplotname+'.png')
	plt.clf()


	print 'u'
	print 'min:', np.min(wind['u'])
	print 'max:', np.max(wind['u'])
	print 'mean:', np.mean(wind['u'])
