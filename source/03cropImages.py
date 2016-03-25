import glob
from tifffile import *
from generalFunctions import *
import numpy as np
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.stats import gaussian_kde
from scipy import interpolate
import shutil

import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

### GUI TO BE IMPLEMENTED

def cropImages( path, worm, channels = ['488nm','CoolLED','561nm'], size = 512 ):

	print( path, worm )
	rawImgsPath = os.path.join( path, worm )

	# load pickle files
	paramsDF = load_data_frame( path, worm + '_01params.pickle' )
	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	gpDF = load_data_frame( path, worm + '_02gonadPos.pickle' )

	# build the movie for each of the input channels
	outpath = os.path.join( path, worm + '_analyzedImages' )

	for channel in channels:
		
		for idx, row in timesDF.iterrows():

			# extract gonad position
			gp = extract_pos( gpDF.ix[ gpDF.tidx == row.tidxRel ].squeeze() )

			# if the gonad is marked
			if not np.sum( np.isnan( gp ) ):

				print( row.tidxRel, row.timesRel, row.fName + channel + '.tif' )

				# copy metadatafile
				if not os.path.isfile( os.path.join( outpath, row.fName + '.txt' ) ):
					shutil.copyfile( os.path.join( rawImgsPath, row.fName + '.txt' ), os.path.join( outpath, row.fName + '.txt' ) )

				if not os.path.isfile( os.path.join( outpath, row.fName + channel + '.tif' ) ):
					# load the Z-stack
					imgs = load_stack( os.path.join( rawImgsPath, row.fName + channel + '.tif') )

					# crop the Z-stack
					cropstack = crop_image( imgs, gp, size )
								
					# save Z-stack
					imsave( os.path.join( outpath, row.fName + channel + '.tif' ), cropstack.astype( np.uint16 ) )

if __name__ == '__main__':

	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

	worms = ['C14','C15','C16','C17','C18','C19']
	wprms = ['C16']

	for w in worms:
	    cropImages( path = path, worm = w , channels = ['488nm','CoolLED','561nm'] )
