import glob
from generalFunctions import *
# from tifffile import *
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import pickle
import os
from PyQt4 import QtGui, QtCore
import sys


class GUI(QtGui.QWidget):
    
	def __init__(self):

		super(GUI, self).__init__()

		self.setWindowTitle( 'Movie making' )

		self.pathDial = 'Y:\\Images'

		self.initUI()

	#-----------------------------------------------------------------------------------------------
	# INITIALIZATION OF THE WINDOW - DEFINE AND PLACE ALL THE WIDGETS
	#-----------------------------------------------------------------------------------------------

	def initUI(self):

		# SET THE GEOMETRY

		mainWindow = QtGui.QGridLayout()
		mainWindow.setSpacing(15)

		self.setLayout(mainWindow)

		# DEFINE ALL WIDGETS AND BUTTONS

		loadExpBtn = QtGui.QPushButton('Load Experiment Folder')
		self.expLbl = QtGui.QLineEdit('')

		loadWormsBtn = QtGui.QLabel('Selected Worms')
		self.wormsLbl = QtGui.QLineEdit('')

		hatchingLbl = QtGui.QLabel('Insert hatching timepoint')
		self.hatchingTp = QtGui.QLineEdit('')

		magnificationLbl = QtGui.QLabel('Insert Magnification')
		self.magnification = QtGui.QComboBox(self)
		self.magnification.addItem('40')
		self.magnification.addItem('60')
		self.magnification.setCurrentIndex(1)

		compressionLbl = QtGui.QLabel('Insert Desired Compression')
		self.compression = QtGui.QComboBox(self)
		for i in [2,4,8,16]:
		    self.compression.addItem(str(i))
		self.compression.setCurrentIndex(1)

		movieLbl = QtGui.QLabel('Select Desired Channels')
		self.channelLED = QtGui.QCheckBox('CoolLED')
		self.channel488 = QtGui.QCheckBox('488nm')
		self.channel561 = QtGui.QCheckBox('561nm')
		self.channelLED.setChecked(True)

		makeMovieBtn = QtGui.QPushButton('Create the Movie!')

		# PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

		mainWindow.addWidget(loadExpBtn,0,0)
		mainWindow.addWidget(self.expLbl,0,1)

		mainWindow.addWidget(loadWormsBtn,1,0)
		mainWindow.addWidget(self.wormsLbl,1,1)

		mainWindow.addWidget(hatchingLbl,2,0)
		mainWindow.addWidget(self.hatchingTp,2,1)

		mainWindow.addWidget(magnificationLbl,3,0)
		mainWindow.addWidget(self.magnification,3,1)

		mainWindow.addWidget(compressionLbl,4,0)
		mainWindow.addWidget(self.compression,4,1)

		mainWindow.addWidget(movieLbl,5,0)
		mainWindow.addWidget(self.channelLED,5,1)
		mainWindow.addWidget(self.channel488,6,1)
		mainWindow.addWidget(self.channel561,7,1)

		mainWindow.addWidget(makeMovieBtn,8,0,1,2)

		self.show()

		# BIND BUTTONS TO FUNCTIONS

		loadExpBtn.clicked.connect(self.selectExp)
		makeMovieBtn.clicked.connect(self.createMovie)

	#-----------------------------------------------------------------------------------------------
	# BUTTON FUNCTIONS
	#-----------------------------------------------------------------------------------------------

	def selectExp(self):

		### store the folder
		self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'Y:\\Images')
		self.exp = self.pathDial.split('\\')[-1]
		self.expLbl.setText(self.exp)

		### update the worms
		wormNames = [ i for i in os.listdir(self.pathDial) if len(i) == 3]
		wormNames.sort()
		wormList = ''
		for i in wormNames:
			wormList += i + ', '
		self.wormsLbl.setText(wormList)

	def createMovie(self):
		
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
		
		if ( len(self.expLbl.text()) == 0 ) or ( len(self.wormsLbl.text()) == 0 ) or ( len(self.hatchingTp.text()) == 0 ):
			QtGui.QMessageBox.about(self,'Error!','No experiment or worms selected!')
			return

		# extract info for creating the movies

		path = self.pathDial

		magnification = int( self.magnification.currentText() )

		scaleFactor = int( self.compression.currentText() )

		channels = []
		if self.channelLED.isChecked():
			channels.append('CoolLED')
		if self.channel488.isChecked():
			channels.append('488nm')
		if self.channel561.isChecked():
			channels.append('561nm')

		worms = self.wormsLbl.text().replace(' ', '')
		if worms[-1] == ',':
			worms = worms[:-1]
		worms = worms.split(',')

		hatchingTidxs = self.hatchingTp.text().replace(' ','')
		if hatchingTidxs[-1] == ',':
			hatchingTidxs = hatchingTidxs[:-1]
		hatchingTidxs = [ int(i) for i in hatchingTidxs.split(',') ]

		# print(path,worms,hatchingTidxs,magnification,scaleFactor,channels)

		# check if there is a mistake in the wormlist
		if len(worms) != len(hatchingTidxs):
			QtGui.QMessageBox.about(self,'Error!','Lengths of hatching times and worms do not match!')
			return

		# create the movie for each worm
		for val in zip( worms, hatchingTidxs ):

			worm = val[0]
			hatchingTidx = val[1]

			rawImgsPath = os.path.join( path, worm )

			### CREATE PARAM AND TIME PICKLE FILE
			# create pickle file if it doesn't exist yet
			if not os.path.isfile( os.path.join( path, worm + '_01params.pickle' ) ):

				parDF = create_params( path, worm, magnification, scaleFactor, hatchingTidx )
				save_data_frame( parDF, path, worm + '_01params.pickle' )

			else:

				parDF = load_data_frame( path, worm + '_01params.pickle' )

			# extract times in absolute values relative to hatching
			if not os.path.isfile( os.path.join( path, worm + '_01times.pickle' ) ):

				tDF = create_times( rawImgsPath, zero=hatchingTidx )
				save_data_frame( tDF, path, worm + '_01times.pickle' )

			else:

				tDF = load_data_frame( path, worm + '_01times.pickle' )

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
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication( sys.argv )
    
    gui = GUI()
    app.setStyle( "plastique" )
    # app.installEventFilter(gui)
    sys.exit( app.exec_() )

