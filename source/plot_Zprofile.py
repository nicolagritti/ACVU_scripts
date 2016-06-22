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

def plotFluorescence( path, worm, ax1, channel = '488nm', cname = '1.p' ):

	print(path, worm)

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gonadPosDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )
	flatField = load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )

	### plot the timeseries
	for idx, trow in timesDF.iterrows():
		print(trow.fName)
		tidx = trow.tidxRel

		currentCells = extract_current_cell_fluo( cellPosDF, cellOutDF, cellFluoDF, tidx )
		currentCell = currentCells.ix[ currentCells.cname == cname ].squeeze()

		gonadPos = extract_pos( gonadPosDF.ix[ gonadPosDF.tidx == tidx ].squeeze() )

		imgpxl = currentCell.imgPxl
		medianCorrection = np.median( ( flatField - darkField ) )
		darkF = crop_image( darkField, gonadPos, 512 )
		flatF = crop_image( flatField, gonadPos, 512 )
		den = np.clip( flatF.astype(np.float) - darkF.astype(np.float), 0, None )

		cellPos = extract_3Dpos( currentCell )
		cellOut = extract_out( currentCell ) * imgpxl / 1000.
		drift = np.array([currentCell.ix['X'+channel]-currentCell.X,currentCell.ix['Y'+channel]-currentCell.Y])

		#load (and correct) image	
		fileName = os.path.join( path, w + '_analyzedImages', trow.fName + channel + '.tif')
		stack = load_stack( fileName )

		data = []
		for img in stack:
			imgCorr = img.astype(np.float) #np.clip( img.astype(np.float) - darkF.astype(np.float), 0, None ) * medianCorrection / den
			imgCorrSmall = imgCorr[cellPos[1]-imgpxl/2:cellPos[1]+imgpxl/2+1,cellPos[0]-imgpxl/2:cellPos[0]+imgpxl/2+1]

			### calculate the signal
			signal = calculate_fluo_intensity( imgCorrSmall, drift, cellOut )
			# area = calculate_area( imgCorrSmall, drift, cellOut )
			data.append(signal)

		ax1.plot( data, '-', lw=2 )


if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,3.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)


	ax1.plot([0,0],[0,2000],'--', color = 'black', lw=1)


	path = 'X:\\Nicola\\160606_noiseTest_beads'
	worms = ['C04']#,'C02']

	for idx, w in enumerate( worms ):
		
		plotFluorescence( path, w, ax1, channel = '561nm' )

	plt.show()
