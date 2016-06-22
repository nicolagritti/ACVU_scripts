# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:31:00 2015

@author: kienle
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from generalFunctions import *
from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

path = 'X:\\Nicola\\160606_noiseTest_beads'
worms = ['C04']#,'C02']

def getPiezo( path, worm, tRow ):

	with open( os.path.join( path,worm,tRow.fName+'.txt' ), 'r') as f:

		line = ''
		while 'Time [ms]' not in line:
		    line = f.readline()

		lines = f.readlines()


	t = np.array( [ np.float(i.strip().split('\t')[0]) for i in lines ] ) * 1000
	_in = np.array( [ np.float(i.strip().split('\t')[1]) for i in lines ] )
	_out = np.array( [ np.float(i.strip().split('\t')[2]) for i in lines ] )

	# plt.plot(t,_in)
	# plt.plot(t,_out)
	# plt.show()
	return ( t, _in, _out )


def plotSingleWormData( path, worm, ax1 ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )



	for idx, tRow in timesDF.iterrows():

		print(tRow.fName)
		(t,_in,_out) = getPiezo( path, worm, tRow )

	ax1.plot( t, _in, '-b', lw=2 )
	ax1.plot( t, _out, '-g', lw=2 )


### setup figure for the timeseries
fig1 = plt.figure(figsize=(5.8,3.8))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(18)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(18)
ax1.set_ylim((-17.2,0.5))
ax1.set_xlim((-5,500))

plotSingleWormData(path,worms[0],ax1)

plt.show()