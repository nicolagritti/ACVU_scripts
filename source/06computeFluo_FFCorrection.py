import glob
from tifffile import *
from generalFunctions import *
import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib.path import Path
import pickle
import scipy.interpolate as ip
from scipy.stats import gaussian_kde
from scipy import interpolate
import shutil

import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def computeFluorescence( path, worm, channels = ['488nm'] ):

	print( path, worm )
	rawImgsPath = os.path.join( path, worm + '_analyzedImages' )

	pickleFile = worm + '_06cellFluo.pickle'
	imgpxl = 50

	# load pickle files
	paramsDF = load_data_frame( path, worm + '_01params.pickle' )
	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gpDF = load_data_frame( path, worm + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )

	# if the cellFluo pickle file already exists, load it, otherwise create a blank one
	if os.path.isfile( os.path.join( path, pickleFile ) ):
	    cellFluoDF = load_data_frame( path, pickleFile )

	else:
	    cellFluoDF = create_cell_fluo( timesDF, list( cellPosDF.keys() ) )

	# load darkField
	darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )

	for channel in channels:

		print(channel)

		# load flatField
		flatField = load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )
		medianCorrection = np.median( ( flatField - darkField ) )

		# for idx, trow in timesDF.ix[ timesDF.tidxRel == 58 ].iterrows(): ### THIS IS TO TEST A PARTICULAR TIMEPOINT!!!
		for idx, trow in timesDF.iterrows():

			# if there is the cropped image
			if os.path.isfile( os.path.join( rawImgsPath, trow.fName + channel + '.tif' ) ):

				print('tidx: ' + str( trow.tidxRel ), trow.fName, trow.timesRel )

				# find all cells labeled in the timepoint
				currentCells = extract_current_cell_out( cellPosDF, cellOutDF, trow.tidxRel )
				gonadPos = extract_pos( gpDF.ix[ gpDF.tidx == trow.tidxRel ].squeeze() )

				# load the stack and perform FF correction
				imgs = load_stack( os.path.join( rawImgsPath, trow.fName + channel + '.tif' ) )
				imgsCorr = flat_field_correction( imgs, darkField, flatField, gonadPos )

				# print('current cells:')
				# print(currentCells)

				### for each labeled cell, calculate the signal
				for jdx, cell in currentCells.iterrows():

					### if there is an outline
					if is_outline_cell( cell ):

						print('detected cell/background with outline: ', cell.cname)

						cellPos = extract_3Dpos( cell )
						cellOut = extract_out( cell )

						### find the XYpos with maximum intensity

						if not cell.cname[:2] == 'b_':

							## this is for the cells - WITH DRIFT CORRECTION
							_range = 10
							signals = np.zeros((2*_range+1,2*_range+1))

							for i in np.arange(2*_range+1):
								for j in np.arange(2*_range+1):
									signals[i,j] = calculate_fluo_intensity( imgsCorr, imgpxl, cellPos+np.array([i-_range,j-_range,0]), cellOut )

							# print( np.where(signals == np.max(signals)) )
							drift = np.array( [ np.where(signals == np.max(signals))[0][0], np.where(signals == np.max(signals))[1][0] ] ) - _range
							# print(drift)
							signal = np.max( signals )

							# ## this is for the cells - NO DRIFT CORRECTION
							# signal = calculate_fluo_intensity( imgsCorr, imgpxl, cellPos, cellOut )
							# drift = [0,0]


						else:

							### this is for backgrounds in which an outline has been drawn
							signal = calculate_fluo_intensity( imgsCorr, imgpxl, cellPos, cellOut )

							drift = [ 0, 0 ]

						# plt.figure()
						# plt.imshow(signals)

						# plt.figure()
						# plt.imshow( imgsCorr[ cell.Z, cell.Y-25:cell.Y+25, cell.X-25:cell.X+25 ], cmap = 'gray', interpolation = 'nearest', clim = [0,4000] )
						# plt.plot(np.append(cell.Xout,cell.Xout[0]),np.append(cell.Yout,cell.Yout[0]),'-o',color='yellow')
						# print(signal)
						# plt.figure()
						# plt.imshow( imgsCorr[ cell.Z, cell.Y+drift[1]-25:cell.Y+drift[1]+25, cell.X+drift[0]-25:cell.X+drift[0]+25 ], cmap = 'gray', interpolation = 'nearest' )
						# plt.plot(np.append(cell.Xout,cell.Xout[0]),np.append(cell.Yout,cell.Yout[0]),'-o',color='red')

					else:

						### if the outline of the background has not been drawn, just take a small area arund its center
						if cell.cname[:2] == 'b_':

							print('detected background without outline: ', cell.cname)
		
							cellPos = extract_3Dpos( cell )

							signal = calculate_fluo_intensity_bckg( imgsCorr, 6, cellPos )

							drift = [ 0, 0 ]

					### update the full dictionary
					newCellFluoDF = update_cell_fluo_DF( cellFluoDF, trow, cell, channel, drift, signal )
					cellFluoDF = newCellFluoDF

				# currentCells = extract_current_cell_fluo( cellFluoDF, trow.tidxRel )
				# print(currentCells)

	save_data_frame( cellFluoDF, path, pickleFile )
	plt.show()


if __name__ == '__main__':


	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

	worms = ['C19']

	for w in worms:
	    computeFluorescence( path = path, worm = w, channels = ['488nm'] )

