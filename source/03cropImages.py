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
from PyQt4 import QtGui, QtCore
import sys
import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

### GUI TO BE IMPLEMENTED!!!

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

		sizeLbl = QtGui.QLabel('Insert Desired Crop Size')
		self.size = QtGui.QLineEdit('512')

		movieLbl = QtGui.QLabel('Select Desired Channels')
		self.channelLED = QtGui.QCheckBox('CoolLED')
		self.channel488 = QtGui.QCheckBox('488nm')
		self.channel561 = QtGui.QCheckBox('561nm')
		self.channelLED.setChecked(True)
		self.channel488.setChecked(True)
		self.channel561.setChecked(True)

		cropBtn = QtGui.QPushButton('Crop the Images!')

		# PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

		mainWindow.addWidget(loadExpBtn,0,0)
		mainWindow.addWidget(self.expLbl,0,1)

		mainWindow.addWidget(loadWormsBtn,1,0)
		mainWindow.addWidget(self.wormsLbl,1,1)

		mainWindow.addWidget(sizeLbl,2,0)
		mainWindow.addWidget(self.size,2,1)

		mainWindow.addWidget(movieLbl,3,0)
		mainWindow.addWidget(self.channelLED,3,1)
		mainWindow.addWidget(self.channel488,4,1)
		mainWindow.addWidget(self.channel561,5,1)

		mainWindow.addWidget(cropBtn,6,0,1,2)

		self.show()

		# BIND BUTTONS TO FUNCTIONS

		loadExpBtn.clicked.connect(self.selectExp)
		cropBtn.clicked.connect(self.cropImages)

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
			if os.path.isfile(os.path.join( self.pathDial, i + '_02gonadPos.pickle' )):
				wormList += i + ', '
		self.wormsLbl.setText(wormList)

	def cropImages( self ):

		'''
			this function crop the images with marked gonad of a certain worm.
			
			Inputs:
			
				- path, worm: folder with the (original) images
			
				- channels: 'LED' for making a moving of transmission images, 
					'488nm' or '561nm' for fluorescence. default: all available channels
					NB: for 'LED' each frame is the mean projection of each 3D stack, 
						for '488nm' and '561nm' it's the maximum projection

				- size: size of the cropped area (default: 512)
			
		'''
		
		if ( len(self.expLbl.text()) == 0 ) or ( len(self.wormsLbl.text()) == 0 ):
			QtGui.QMessageBox.about(self,'Error!','No experiment or worms selected!')
			return

		# extract info for image cropping
		path = self.pathDial

		size = int( self.size.text() )

		channels = []
		if self.channelLED.isChecked():
			channels.append('CoolLED')
		if self.channel488.isChecked():
			channels.append('488nm')
		if self.channel561.isChecked():
			channels.append('488nm')

		worms = self.wormsLbl.text().replace(' ', '')
		if worms[-1] == ',':
			worms = worms[:-1]
		worms = worms.split(',')

		# print(path,worms,channels,size)

		# check if there is a mistake in the wormlist
		for worm in worms:
			if not os.path.isfile(os.path.join( self.pathDial, worm + '_02gonadPos.pickle' )):
				QtGui.QMessageBox.about(self,'Error!','Worm ' + worm + ' has no marked gonad!')
				return

		# crop the images for each worm
		for worm in worms:

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
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication( sys.argv )
    
    gui = GUI()
    app.setStyle( "plastique" )
    # app.installEventFilter(gui)
    sys.exit( app.exec_() )
