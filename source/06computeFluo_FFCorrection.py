

"""
PyQt seam cells analysis GUI

NB: the python package tifffile (Gohlke) needs to be installed.

author: Nicola Gritti
last edited: June 2015
"""

import sys
from tifffile import *
from generalFunctions import *
import pickle
import os
from PyQt4 import QtGui, QtCore
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends import qt_compat
import glob
import pandas as pd
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE

class GUI(QtGui.QWidget):
    
	def __init__(self):

		super(GUI, self).__init__()

		self.setWindowTitle( 'Outline Cells' )
		self.cellNames = ['1.p','4.a','1.pp','4.aa','1.ppa','1.ppp','4.aaa','4.aap','b_1','b_4']
		self.initUI()
    
    #-----------------------------------------------------------------------------------------------
    # INITIALIZATION OF THE WINDOW - DEFINE AND PLACE ALL THE WIDGETS
    #-----------------------------------------------------------------------------------------------

	def initUI(self):
		# SET THE GEOMETRY

		mainWindow = QtGui.QVBoxLayout()
		mainWindow.setSpacing(15)

		fileBox = QtGui.QHBoxLayout()
		spaceBox1 = QtGui.QHBoxLayout()
		rawDataBox = QtGui.QHBoxLayout()

		mainWindow.addLayout(fileBox)
		mainWindow.addLayout(spaceBox1)
		mainWindow.addLayout(rawDataBox)

		Col1 = QtGui.QGridLayout()
		Col2 = QtGui.QHBoxLayout()
		Col3 = QtGui.QVBoxLayout()

		rawDataBox.addLayout(Col1)
		rawDataBox.addLayout(Col2)
		rawDataBox.addLayout(Col3)

		self.setLayout(mainWindow)

		# DEFINE ALL WIDGETS AND BUTTONS

		loadBtn = QtGui.QPushButton('Load DataSet')
		computeBtn = QtGui.QPushButton('Compute 488 signal (with autmatic drift correction)')
		saveBtn = QtGui.QPushButton('Save data (F12)')

		tpLbl = QtGui.QLabel('Relative Tp:')
		slLbl = QtGui.QLabel('Slice:')
		fNameLbl = QtGui.QLabel('File name:')

		self.tp = QtGui.QSpinBox(self)
		self.tp.setValue(-5)
		self.tp.setMaximum(100000)

		self.sl = QtGui.QSpinBox(self)
		self.sl.setValue(0)
		self.sl.setMaximum(100000)

		self.fName = QtGui.QLabel('')

		self._488nmBtn = QtGui.QRadioButton('488nm')
		self._561nmBtn = QtGui.QRadioButton('561nm')
		self.CoolLEDBtn = QtGui.QRadioButton('CoolLED')

		self.sld1 = QtGui.QSlider(QtCore.Qt.Vertical, self)
		self.sld1.setMaximum(2**16-1)
		self.sld1.setValue(0)
		self.sld2 = QtGui.QSlider(QtCore.Qt.Vertical, self)
		self.sld2.setMaximum(2**16)
		self.sld2.setValue(2**16-1)

		self.fig1 = Figure((8.0, 8.0), dpi=100)
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax1 = self.fig1.add_subplot(111)
		self.canvas1 = FigureCanvas(self.fig1)
		self.canvas1.setFocusPolicy( QtCore.Qt.ClickFocus )
		self.canvas1.setFocus()
		self.canvas1.setFixedSize(QtCore.QSize(600,600))
		self.canvas1.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

		self.fig2 = Figure((4.0, 4.0), dpi=100)
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2 = FigureCanvas(self.fig2)
		self.canvas2.setFixedSize(QtCore.QSize(300,300))
		self.canvas2.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )

		# PLACE ALL THE WIDGET ACCORDING TO THE GRIDS

		fileBox.addWidget(loadBtn)
		fileBox.addWidget(saveBtn)

		spaceBox1.addWidget(self.HLine())

		Col1.addWidget(tpLbl, 0, 0)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(self.tp, 0, 1)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(slLbl, 1, 0)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(self.sl, 1, 1)#, 1, 1, Qt.AlignTop)
		Col1.addWidget(fNameLbl, 2, 0)
		Col1.addWidget(self.fName, 2, 1)
		Col1.addWidget(self._488nmBtn, 3, 0 )
		Col1.addWidget(self._561nmBtn, 4, 0 )
		Col1.addWidget(self.CoolLEDBtn, 5, 0 )

		Col2.addWidget(self.sld1)
		Col2.addWidget(self.sld2)
		Col2.addWidget(self.canvas1)

		Col3.addWidget(self.canvas2)
		Col3.addWidget(computeBtn)

		self.setFocus()
		self.show()

		# BIND BUTTONS TO FUNCTIONS

		loadBtn.clicked.connect(self.selectWorm)
		saveBtn.clicked.connect(self.saveData)

		self.tp.valueChanged.connect(self.loadNewStack)
		self.sl.valueChanged.connect(self.updateAllCanvas)
		self.sld1.valueChanged.connect(self.updateBC)
		self.sld2.valueChanged.connect(self.updateBC)

		self._488nmBtn.toggled.connect(self.radio488Clicked)
		self._561nmBtn.toggled.connect(self.radio561Clicked)
		self.CoolLEDBtn.toggled.connect(self.radioCoolLEDClicked)

		computeBtn.clicked.connect(self.computeFluo)

		self.fig1.canvas.mpl_connect('button_press_event',self.onMouseClickOnCanvas1)        
		self.fig1.canvas.mpl_connect('scroll_event',self.wheelEvent)        

	#-----------------------------------------------------------------------------------------------
	# FORMATTING THE WINDOW
	#-----------------------------------------------------------------------------------------------

	def center(self):
        
		qr = self.frameGeometry()
		cp = QtGui.QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())
        
	def HLine(self):

		toto = QtGui.QFrame()
		toto.setFrameShape(QtGui.QFrame.HLine)
		toto.setFrameShadow(QtGui.QFrame.Sunken)
		return toto

	def VLine(self):
	    
	    toto = QtGui.QFrame()
	    toto.setFrameShape(QtGui.QFrame.VLine)
	    toto.setFrameShadow(QtGui.QFrame.Sunken)
	    return toto

	def heightForWidth(self, width):

		return width

	#-----------------------------------------------------------------------------------------------
	# BUTTON FUNCTIONS
	#-----------------------------------------------------------------------------------------------

	def selectWorm(self):

		### store the folders
		self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101')#'Y:\\Images')
		self.worm = self.pathDial.split("\\")[-1].split('_')[0]
		self.path = os.path.dirname( self.pathDial )
		self.setWindowTitle('Outline Cells - ' + self.pathDial)

		### give error message if there is no CoolLED movie in the selected folder
		if not os.path.isfile( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) ):
			QtGui.QMessageBox.about(self,'Warning!','There is no movie in this folder! Create a movie first!')
			return

		### load parameters and times dataframes
		self.paramsDF = load_data_frame( self.path, self.worm + '_01params.pickle' )
		self.timesDF = load_data_frame( self.path, self.worm + '_01times.pickle' )
		self.gpDF = load_data_frame( self.path, self.worm + '_02gonadPos.pickle' )
		self.cellPosDF = load_data_frame( self.path, self.worm + '_04cellPos.pickle' )
		self.cellOutDF = load_data_frame( self.path, self.worm + '_05cellOut.pickle' )

		# extract some info
		self.compression = self.paramsDF.compression
		self.hatchingtidx = int( self.paramsDF.tidxHatch )

		### if the cellOutline pickle file already exists, load it, otherwise create a blank one
		if os.path.isfile( os.path.join(self.path, self.worm + '_06cellFluo.pickle' ) ):
			self.cellFluoDF = load_data_frame( self.path, self.worm + '_06cellFluo.pickle' )

		else:
			self.cellFluoDF = create_cell_fluo( self.timesDF, self.cellNames )

		# detect available channels
		self.channels = []
		chns = ['CoolLED','488nm','561nm']
		for c in chns:

			if os.path.isfile( os.path.join( self.pathDial, c + '_movie.tif' ) ):

				self.channels.append(c)
		self.currentChannel = 'CoolLED'

		### detect size of the cropped images
		tp = first_tidx_pos_all_cells( self.cellPosDF )
		tRow = self.timesDF.ix[ self.timesDF.tidxRel == tp ].squeeze()
		fileName = os.path.join( self.pathDial, tRow.fName + 'CoolLED.tif')
		firststack = load_stack( fileName )
		self.cropsize = firststack.shape[1]
		self.nslices = firststack.shape[0]

		### intialize canvases
		self.initializeCanvas1()
		self.initializeCanvas2()

		### extract current cells already labeled
		self.currentCells = extract_current_cell_out( self.cellPosDF, self.cellOutDF, self.tp.value() )
		self.analyzedCell = '---'

		### update the text of the fileName
		self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == tp, 'fName' ].values[0])

		### set the timepoint to the hatching time
		self.tp.setMinimum(np.min(self.timesDF.tidxRel))
		self.tp.setMaximum(np.max(self.timesDF.tidxRel))

		### set the max slice number
		self.sl.setMaximum( self.nslices-1 )

		if tp != self.tp.value():
			self.tp.setValue( tp )
		else:
			self.loadNewStack()

		self.CoolLEDBtn.setChecked(True)    # this uppdates the canvases once more

		# self.pathDial.show()
		self.setFocus()

	def loadNewStack(self):

		# before changing timepoint, print labeled cells and check if they are OK
		print( 'cells labeled:\n ', self.currentCells )

		# print(self.fList['gfp'][self.tp.value()])
		tRow = self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value() ].squeeze()

		### update the text of the fileName
		self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])

		print( 'Loading... ', self.pathDial, tRow.fName )

		# calculate the max value of the previous stack
		try:
			prevmax = np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] )
		# if it's the first time a stack is to be loaded (so if there is no previous stack), set it to zero
		except:
			prevmax = 0

		# load all the available stacks - this is the slowest part of the code!!!
		self.stacks = {}
		for ch in self.channels:
			fileName = os.path.join( self.pathDial, tRow.fName + ch + '.tif')
			if os.path.isfile( fileName ):
				# print(MultiImage('X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101\\C02_analyzedImages\\Z003_488nm.tif'))
				# print(fileName, MultiImage( fileName ))
				# self.stacks[ch] = MultiImage( fileName )
				self.stacks[ch] = load_stack( fileName )
			# if there are no files for the timepoint, create a blank image
			else:
				self.stacks[ch] = prevmax*np.ones((self.nslices,self.cropsize,self.cropsize))

		### extract current cells already labeled
		self.currentCells = extract_current_cell_out( self.cellPosDF, self.cellOutDF, self.tp.value() )
		self.currentDriftedCells = extract_current_cell_fluo( self.cellFluoDF, self.tp.value() )
		print(self.currentCells)

		# if there are cells labeled and if the previously currently analyzed cell is present, set it as the currently labeled cell and select the right slice
		if len(self.currentCells) > 0:

			### update currently analyzed cell
			if self.analyzedCell not in list( self.currentCells.cname ):
				self.analyzedCell = self.currentCells.cname[0]

			### update slice
			self.sl.setValue( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'Z' ] )

		# if the BC bound are different, the BCsliderMinMax will automatically update all canvases. Otherwise, manually update it!
		newmax = np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] )
		if prevmax != newmax:
			self.setBCslidersMinMax()            
		else:
			self.updateAllCanvas()

	def saveData(self):

		save_data_frame( self.cellFluoDF, self.path, self.worm + '_06cellFluo.pickle' )
		self.setFocus()
	    
	def updateAllCanvas(self):
		self.updateCanvas1()
		self.updateCanvas2()
        
	def radio488Clicked(self, enabled):
		print('radio 488 clicked')

		if enabled:
			if '488nm' in self.channels:
				self.currentChannel = '488nm'
				self.setFocus()
				self.updateAllCanvas()
			else:
				self.CoolLEDBtn.setChecked(True)
				QtGui.QMessageBox.about(self, 'Warning', 'No 488nm channel!')

	def radio561Clicked(self, enabled):
	    print('radio 561 clicked')

	    if enabled:
	        if '561nm' in self.channels:
	            self.currentChannel = '561nm'
	            self.setFocus()
	            self.updateAllCanvas()
	        else:
	            self.CoolLEDBtn.setChecked(True)
	            QtGui.QMessageBox.about(self, 'Warning', 'No 561nm channel!')

	def radioCoolLEDClicked(self, enabled):
		print('radio LED clicked')

		if enabled:
			self.currentChannel = 'CoolLED'
			self.setFocus()
			self.updateAllCanvas()

	def updateBC(self):
		# change brightness and contrast
		self.imgplot1.set_clim( self.sld1.value(), self.sld2.value() )  
		self.imgplot2.set_clim( self.sld1.value(), self.sld2.value() )  
		self.canvas1.draw()
		self.canvas2.draw()

	#-----------------------------------------------------------------------------------------------
	# DEFAULT FUNCTION FOR KEY AND MOUSE PRESS ON WINDOW
	#-----------------------------------------------------------------------------------------------

	def keyPressEvent(self, event):

		# print(event.key())

		# change timepoint
		if event.key() == QtCore.Qt.Key_Right:
			self.changeSpaceTime( 'time', +1 )

		elif event.key() == QtCore.Qt.Key_Left:
			self.changeSpaceTime( 'time', -1 )

		# change slice
		elif event.key() == QtCore.Qt.Key_Up:
			self.changeSpaceTime( 'space', +1 )
		    
		elif event.key() == QtCore.Qt.Key_Down:
			self.changeSpaceTime( 'space', -1 )

		elif event.key() == QtCore.Qt.Key_Space:

			if len( self.currentCells ) > 0:
				idx = np.mod( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].index + 1, len(self.currentCells) )
				self.analyzedCell = self.currentCells.cname.values[idx][0]
				self.sl.setValue( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'Z' ] )

				self.updateAllCanvas()

		self.setFocus()

	def wheelEvent(self,event):
		if self.canvas1.underMouse():
			step = event.step
		else:          
			step = event.delta()/abs(event.delta())
		self.sl.setValue( self.sl.value() + step) 

	#-----------------------------------------------------------------------------------------------
	# ADDITIONAL FUNCTIONS FOR KEY AND MOUSE PRESS ON CANVASES
	#-----------------------------------------------------------------------------------------------

	def onMouseClickOnCanvas1(self, event):
		print('move center of the cell - TO BE IMPLEMENTED')

		self.updateCanvas1()
		self.setFocus()
		# print(event.button,event.xdata,event.ydata)

	#-----------------------------------------------------------------------------------------------
	# UTILS
	#-----------------------------------------------------------------------------------------------

	def setBCslidersMinMax(self):
		self.sld1.setMaximum( np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] ) )
		self.sld1.setMinimum(0)
		self.sld2.setMaximum (np.max( [ np.max(self.stacks[ch]) for ch in self.channels ] ) )
		self.sld2.setMinimum(0)

	def initializeCanvas1(self):
		print('initializing canvas1')

		self.fig1.clf()
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax1 = self.fig1.add_subplot(111)
		self.canvas1.draw()

		# plot the first blank image with the right size
		self.ax1.cla()
		self.imgplot1 = self.ax1.imshow( np.zeros((50,50)), cmap = 'gray', interpolation = 'nearest' )

		# remove the white borders
		self.ax1.autoscale(False)
		self.ax1.axis('Off')
		self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

		# plot cell pos and name
		self.text_c1 = []
		self.out_c1 = []
		self.outDrift_c1 = []
		self.centerDrift_c1 = []

		# redraw the canvas
		self.canvas1.draw()
		self.setFocus()

	def updateCanvas1(self):
		print('updating canvas1')

		# clear cell text and points
		# print(self.text1,self.points1)
		for text in self.text_c1:
			text.remove()
		self.text_c1 = []

		for points in self.out_c1:
			self.ax1.lines.remove(points)
		self.out_c1 = []

		for points in self.outDrift_c1:
			self.ax1.lines.remove(points)
		self.outDrift_c1 = []

		for points in self.centerDrift_c1:
			self.ax1.lines.remove(points)
		self.centerDrift_c1 = []

		# if no cells labeled, leave the image blank
		if len( self.currentCells ) == 0:
			self.imgplot1.set_data( np.ones((50,50)) )
			self.canvas1.draw()
			return

		# extract current cell data
		pos = extract_3Dpos( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

		# plot the image
		imgpxl = 50
		self.imgplot1.set_data( self.stacks[self.currentChannel][self.sl.value(),pos[1]-imgpxl/2:pos[1]+imgpxl/2+1,pos[0]-imgpxl/2:pos[0]+imgpxl/2+1] )

		# change brightness and contrast
		self.imgplot1.set_clim( self.sld1.value(), self.sld2.value() )  

		# print cell name
		if pos[2] == self.sl.value():
			self.text_c1.append( self.ax1.text( 1, 2, self.analyzedCell, color='yellow', size='medium', alpha=.8,
					rotation=0, fontsize = 20 ) )

		### draw outline
		outline = extract_out( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

		# print(self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze())

		if len( outline ) > 1:
			outline = np.vstack( [ outline, outline[0] ] )

		self.plot_c1.append( self.ax1.plot( outline[:,0], outline[:,1], '-x', color='yellow', ms=2, alpha=1., lw = .5 )[0] )

		# draw drifted outline
		drift


		# # redraw the canvas
		self.canvas1.draw()
		self.setFocus()

	def initializeCanvas2(self):
		print('initializing canvas2')

		self.fig2.clf()
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2.draw()

		# plot the first blank image with the right size
		self.ax2.cla()
		self.imgplot2 = self.ax2.imshow( np.zeros((self.cropsize,self.cropsize)), cmap = 'gray' )

		# remove the white borders
		self.ax2.autoscale(False)
		self.ax2.axis('Off')
		self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)

		# plot cell pos and name
		self.text_c2 = []
		self.plot_c2 = []

		# redraw the canvas
		self.canvas2.draw()
		self.setFocus()

	def updateCanvas2(self):
		print('updating canvas2')

		# plot the image
		self.imgplot2.set_data( self.stacks[self.currentChannel][self.sl.value()] )

		# change brightness and contrast
		self.imgplot2.set_clim( self.sld1.value(), self.sld2.value() )  

		# clear cell text and points
		# print(self.text1,self.points1)
		for text in self.text_c2:
			text.remove()
		self.text_c2 = []

		for points in self.plot_c2:
			self.ax2.lines.remove(points)
		self.plot_c2 = []

		# extract current cell data
		for idx, cell in self.currentCells.iterrows():

			if cell.Z == self.sl.value():

				color = 'red'
				if cell.cname == self.analyzedCell:
					color = 'yellow'

				self.text_c2.append( self.ax2.text( cell.X, cell.Y + 18, cell.cname, color=color, size='medium', alpha=.8,
						rotation=0 ) )
				self.plot_c2.append( self.ax2.plot( cell.X, cell.Y, 'o', color=color, alpha = .8, mew = 0 )[0] )


		# redraw the canvas
		self.canvas2.draw()
		self.setFocus()

	def changeSpaceTime(self, whatToChange, increment):


		if whatToChange == 'time':

			newCellOutDF = update_cell_out_DF( self.currentCells, self.cellOutDF, self.tp.value() )
			self.cellOutDF = newCellOutDF

			# if they are OK (and not going to negative times), change timepoint
			self.tp.setValue( self.tp.value() + increment )

		if whatToChange == 'space':
			self.sl.setValue( self.sl.value() + increment )

	def computeFluo( self ):#path, worm, channels = ['488nm'] ):

		path = self.path
		worm = self.worm

		print( path, worm )
		rawImgsPath = os.path.join( path, worm + '_analyzedImages' )

		imgpxl = 50

		# load pickle files
		paramsDF = self.paramsDF
		timesDF = self.timesDF
		gpDF = self.gpDF
		cellPosDF = self.cellPosDF
		cellOutDF = self.cellOutDF
		cellFluoDF = self.cellFluoDF

		currentCells = self.currentCells

		# load darkField
		darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )

		channels = [ '488nm' ]
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
						self.cellFluoDF = newCellFluoDF

					# currentCells = extract_current_cell_fluo( cellFluoDF, trow.tidxRel )
					# print(currentCells)
		self.updateCanvas1()


	# if __name__ == '__main__':


	# 	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

	# 	worms = ['C19']

	# 	for w in worms:
	# 	    computeFluorescence( path = path, worm = w, channels = ['488nm'] )



if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication(sys.argv)
    
    gui = GUI()
    app.setStyle("plastique")
    # app.installEventFilter(gui)
    sys.exit(app.exec_())
    





# import glob
# from tifffile import *
# from generalFunctions import *
# import numpy as np
# import os.path
# import matplotlib.pyplot as plt
# from matplotlib.path import Path
# import pickle
# import scipy.interpolate as ip
# from scipy.stats import gaussian_kde
# from scipy import interpolate
# import shutil

# import os
# import matplotlib as mpl
# mpl.rcParams['pdf.fonttype'] = 42

# if __name__ == '__main__':


# 	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

# 	worms = ['C19']

# 	for w in worms:
# 	    computeFluorescence( path = path, worm = w, channels = ['488nm'] )

