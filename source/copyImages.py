
import glob
from generalFunctions import *
import numpy as np
import shutil
import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

path  = 'X:\\Simone\\160420_LIN12_GFP_hist_mCherry'
worms = ['C17','C19','C20','C21','C22','C23','C24','C26','C27','C28','C29','C30']
channels = ['CoolLED','488nm','561nm']
for worm in worms:

	print( path, worm )
	rawImgsPath = os.path.join( path, worm )

	# build the movie for each of the input channels
	outpath = os.path.join( path, worm )

	fName = 'z089_'
	outfName = 'z090_'

		# copy metadatafile
	shutil.copyfile( os.path.join( rawImgsPath, fName + '.txt' ), os.path.join( outpath, outfName + '.txt' ) )

	for channel in channels:

		# load the Z-stack
		imgs = load_stack( os.path.join( rawImgsPath, fName + channel + '.tif') )

		# save Z-stack
		imsave( os.path.join( outpath, outfName + channel + '.tif' ), imgs.astype( np.uint16 ) )
