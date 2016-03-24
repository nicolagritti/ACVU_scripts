import glob
from generalFunctions import *
from tifffile import *
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pickle
import os


def makeMovie( path, worm, hatchingTidx, magnification = 60, channels = ['CoolLED'], scaleFactor = 4 ):
	'''
		this function builds a movie from downsized images of a certain worm.
		
		Inputs:
		
			- path, worm: folder with the (original) images
		
			- hatching time: index of the first image after hatching

			- magnification: objective used (default: 60X)

			- channels: 'LED' for making a moving of transmission images, 
				'488nm' or '561nm' for fluorescence. default: only 'LED'
				NB: for 'LED' each frame is the mean projection of each 3D stack, 
					for '488nm' and '561nm' it's the maximum projection

			- scaleFactor: downsizing parameter (default: 4)
		
	'''


	rawImgsPath = os.path.join( path, worm )
	print(rawImgsPath)

	### CREATE PARAM AND TIME PICKLE FILE
	# create pickle file if it doesn't exist yet
	if not os.path.isfile( os.path.join( path, worm + '_01params.pickle' ) ):

		parDF = create_params( path, worm, magnification, scaleFactor, hatchingTidx )
		save_data_frame( parDF, path, worm + '_01params.pickle' )

	else:

		parDF = load_data_frame( path, worm + '_01params.pickle' )

	# print(parDF)

	# extract times in absolute values relative to hatching
	if not os.path.isfile( os.path.join( path, worm + '_01times.pickle' ) ):

		tDF = create_times( rawImgsPath, zero=hatchingTidx )
		save_data_frame( tDF, path, worm + '_01times.pickle' )

	else:

		tDF = load_data_frame( path, worm + '_01times.pickle' )

	# print(tDF)

	# create new folder for analyzed images
	outpath = os.path.join( path, worm + '_analyzedImages' )
	os.makedirs( outpath, exist_ok = True)	

	### BUILD THE MOVIE FOR EACH OF THE INPUT CHANNELS

	for channel in channels:
		
		### CREATE THE MOVIE WITHOUT TIMESTAMP

		continueMakingMovie = not os.path.isfile( os.path.join( outpath, channel + '_movie.tif' ) )
				
		if continueMakingMovie:

			# load the file list with all the raw images
			flist = glob.glob( os.path.join( rawImgsPath,  '*' + channel + '.tif')	)
			flist.sort()
			# print(len(flist))

			# create a blank movie list
			movieFinal = []

			# minimum and maximum values to optimize the dynamic range of the timeframes
			_min = 2**16
			_max = 0

			for idx, f in enumerate( flist ):
				print(f)

				# load the Z-stack
				imgs = load_stack( f )

				# downsized Z-stack
				smallimgs = downsizeStack( imgs, scaleFactor )

				# compute the mean or maximum projection depending on the channel input
				if channel == 'CoolLED':
					frame = np.mean(smallimgs,0).astype(np.uint16)
				else:
					frame = np.max(smallimgs,0).astype(np.uint16)

				# append the projection to the movie
				movieFinal.append( frame )

				# update the minimum and maximum value for brightness and contrast based on previous and current frames of the movie
				_min = np.min(frame) * ( np.min(frame) < _min ) + _min * ( np.min(frame) >= _min )
				_max = np.max(frame) * ( np.max(frame) > _max ) + _max * ( np.max(frame) <= _max )

			# change b&c
			movieFinal = ( np.array( movieFinal ) - _min ) / ( _max - _min )

			# convert images to 8 bit
			movieFinal = ( ( 2**8 - 1. ) * movieFinal  ).astype( np.uint8 )

			# save the movie without timestamps
			imsave( os.path.join( outpath, channel + '_movie.tif' ), np.array( movieFinal ) )

		### CREATE AND SAVE THE MOVIE WITH TIMESTAMPS
		
		continueMakingMovie = not os.path.isfile( os.path.join( outpath, channel + '_movieWithTime.tif' ) )
				
		if continueMakingMovie:

			if os.path.isfile( os.path.join( outpath, channel + '_movie.tif' ) ):
				movieFinal = load_stack( os.path.join( outpath, channel + '_movie.tif' ) )

			movieFinalWithTime = np.zeros( ( len(movieFinal), movieFinal[0].shape[0], movieFinal[0].shape[1] ) ).astype( np.uint8 )

			for idx, img in enumerate( movieFinal ):

				# create the PythonImageLibrary object of the frame to write the time text
				imgpil = Image.fromarray( movieFinal[idx], 'L' )

				# write the time
				font = ImageFont.truetype( "calibri.ttf", 60 )
				draw = ImageDraw.Draw(imgpil)
				draw.text((0,0),'%d h' % np.floor( tDF.timesRel[idx] ),fill='white',font=font)

				movieFinalWithTime[idx] = np.asarray( imgpil )

			imsave( os.path.join( outpath, channel + '_movieWithTime.tif' ), np.array( movieFinalWithTime ) )



if __name__ == '__main__':

	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

	worms = [ 'C15', 'C16', 'C17', 'C18', 'C19' ]

	hatchingTidx = [3,1,5,4,8]

	magnification = 60

	for val in zip( worms, hatchingTidx ):
	    makeMovie( path = path, worm = val[0], hatchingTidx = val[1], 
	    	magnification = magnification, channels = [ 'CoolLED' ], scaleFactor = 4 )
