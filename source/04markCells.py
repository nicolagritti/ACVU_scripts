

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
        
        self.setWindowTitle( 'Label Cells' )
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
        saveBtn = QtGui.QPushButton('Save data (F12)')
        
        tpLbl = QtGui.QLabel('Relative Tp:')
        slLbl = QtGui.QLabel('Slice:')
        fNameLbl = QtGui.QLabel('File name:')
        
        self.tp = QtGui.QSpinBox(self)
        self.tp.setValue(0)
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

        self.cellTbl = QtGui.QTableWidget()

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

        Col3.addWidget(self.cellTbl)
        Col3.addWidget(self.canvas2)
        
        self.setFocus()
        self.show()
        
        # BIND BUTTONS TO FUNCTIONS
        
        loadBtn.clicked.connect(self.selectWorm)
        saveBtn.clicked.connect(self.saveData)

        self.tp.valueChanged.connect(self.loadNewStack)
        self.sl.valueChanged.connect(self.updateAllCanvas)
        self.sld1.valueChanged.connect(self.updateAllCanvas)
        self.sld2.valueChanged.connect(self.updateAllCanvas)

        self._488nmBtn.toggled.connect(self.radioClicked)
        self._561nmBtn.toggled.connect(self.radioClicked)
        self.CoolLEDBtn.toggled.connect(self.radioClicked)

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
        self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'Y:\\Images')
        self.worm = self.pathDial.split("\\")[-1].split('_')[0]
        self.path = os.path.dirname( self.pathDial )
        self.setWindowTitle('Mark Cells - ' + self.pathDial)
        
        ### give error message if there is no CoolLED movie in the selected folder
        if not os.path.isfile( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) ):
            QtGui.QMessageBox.about(self,'Warning!','There is no movie in this folder! Create a movie first!')
            return

        ### load parameters and times dataframes
        self.paramsDF = load_data_frame( self.path, self.worm + '_01params.pickle' )
        self.timesDF = load_data_frame( self.path, self.worm + '_01times.pickle' )
        self.gpDF = load_data_frame( self.path, self.worm + '_02gonadPos.pickle' )

        # extract some info
        self.compression = self.paramsDF.compression
        self.hatchingtidx = int( self.paramsDF.tidxHatch )

        ### if the cellPos pickle file already exists, load it, otherwise create a blank one
        if os.path.isfile( os.path.join(self.path, self.worm + '_04cellPos.pickle' ) ):
            self.cellPosDF = load_data_frame( self.path, self.worm + '_04cellPos.pickle' )
        
        else:
            self.cellPosDF = create_cell_pos( self.timesDF, self.cellNames )

        ### load all movies (without timestamps, we will add it later on)
        self.LEDmovie = load_stack( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) )

        ### set the timepoint to the hatching time
        self.tp.setMinimum(np.min(self.timesDF.tidxRel))
        self.tp.setMaximum(np.max(self.timesDF.tidxRel))
        self.tp.setValue( 0 )

        ### extract current cells already labeled
        self.currentCells = extract_current_cell_pos( self.cellPosDF, self.tp.value() )

        # detect available channels
        self.channels = []
        chns = ['CoolLED','488nm','561nm']
        for c in chns:

            if os.path.isfile( os.path.join( self.pathDial, c + '_movie.tif' ) ):

                self.channels.append(c)
        self.currentChannel = self.channels[0]

        ### update the text of the fileName
        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])
        
        self.loadNewStack()

        # self.pathDial.show()
        self.updateAllCanvas()
        self.setFocus()

    def loadNewStack(self):

        # print(self.fList['gfp'][self.tp.value()])
        tRow = self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value() ].squeeze()

        print( 'Loading... ', self.pathDial, tRow.fName )

        # load all the available stacks
        self.stacks = {}
        for ch in self.channels:
            fileName = os.path.join( self.pathDial, tRow.fName + ch + '.tif')
            if os.path.isfile( fileName ):
                self.stacks[ch] = load_stack( fileName )

        if len( self.stacks.keys() ) > 0:
            # print(self.stacks.keys(), self.stacksStraight)
            self.sl.setMaximum(self.stacks[self.currentChannel].shape[0]-1)

            self.setBCslidersMinMax()
        
        ### update the text of the fileName
        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])

        ### extract current cells already labeled
        self.currentCells = extract_current_cell_pos( self.cellPosDF, self.tp.value() )

        # self.updateTable()
        self.updateAllCanvas()

    def saveData(self):
        
        if self.checkConsistencyCellNames():
            save_data_frame( self.cellPosDF, self.path, self.worm + '_04cellPos.pickle' )        
        else:
            QtGui.QMessageBox.about(self,'Warning!','There is a mistake in the cell labels!')
        self.setFocus()
        
    def updateAllCanvas(self):
        self.updateRadioBtn()
        self.updateCanvas1()
        self.updateCanvas2()
        
    def radioClicked(self):
        if self._488nmBtn.isChecked():
            if '488nm' in self.channels:
                self.currentChannel = '488nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 488nm channel!')
        elif self._561nmBtn.isChecked():
            if '561nm' in self.channels:
                self.currentChannel = '561nm'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No 561nm channel!')
        elif self.CoolLEDBtn.isChecked():
            if 'CoolLED' in self.channels:
                self.currentChannel = 'CoolLED'
            else:
                QtGui.QMessageBox.about(self, 'Warning', 'No CoolLED channel!')
        self.setBCslidersMinMax()
        self.resetBC()
        self.setFocus()
        self.updateAllCanvas()

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

        # key press on cropped image
        if self.canvas1.underMouse():
            self.onKeyPressOnCanvas1(event)
            
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

    def onKeyPressOnCanvas1(self, event):
        
        motherCells = [ QtCore.Qt.Key_1, QtCore.Qt.Key_4, QtCore.Qt.Key_B ]
        daughterCells = [ QtCore.Qt.Key_A, QtCore.Qt.Key_P ]

        # find the position of the cursor relative to the image in pixel
        imgshape = self.stacks[self.currentChannel][self.sl.value()].shape
        canshape = self.canvas1.size()
        cf = imgshape[0]/canshape.width()
        refpos = self.canvas1.mapFromGlobal(QtGui.QCursor.pos())
        refpos = np.array([ int( refpos.x() * cf ), int( refpos.y() * cf )])
        refpos = np.append(refpos,self.sl.value())

        ### find the closest cell to the cursor
        idx = closer_cell( refpos.astype( np.uint16 ), self.currentCells )

        ### assign the name to the cell
        if any( [ event.key() == cn for cn in motherCells ] ):
            # if labeling bckg, add the 1 or 2 to the name
            if self.currentCells.ix[ idx, 'cname' ] == 'b_':
                self.currentCells.ix[ idx, 'cname' ] += QtGui.QKeySequence(event.key()).toString().lower()
            else:
                # if not, rename the cell from scratch
                if event.key() == QtCore.Qt.Key_B:
                    self.currentCells.ix[ idx, 'cname' ] = QtGui.QKeySequence(event.key()).toString().lower() + '_'
                else:
                    self.currentCells.ix[ idx, 'cname' ] = QtGui.QKeySequence(event.key()).toString().lower() + '.'

        # add the anterior/posterior to the cell name
        elif any( [ event.key() == cp for cp in daughterCells ] ):
            # don't do it if labeling background (bckg doesn't have a/p!)
            if self.currentCells.ix[ idx, 'cname' ] == 'b_':
                return
            else:
                self.currentCells.ix[ idx, 'cname' ] += QtGui.QKeySequence(event.key()).toString().lower()
        
        # remove the last entry of the name with backspace                
        elif event.key() == QtCore.Qt.Key_Backspace:
            self.currentCells.ix[ idx, 'cname' ] = self.currentCells.ix[ idx, 'cname' ][:-1]

        self.updateCanvas1()
        self.setFocus()

    def onMouseClickOnCanvas1(self, event):

        refpos = np.array( [ event.xdata, event.ydata, self.sl.value() ] )  

        if event.button == 1:

            # create an empty cell in the currentCells df: the only entries are tidx, xyzpos and cname
            newcell = create_single_cell_pos( refpos.astype(np.uint16), self.tp.value() )
            self.currentCells = pd.concat( [ self.currentCells, newcell ] )
            
        elif event.button == 3:

            # remove a cell (the closest to the cursor at the moment of right-click)
            idx = closer_cell( refpos.astype(np.uint16), self.currentCells )
            self.currentCells = self.currentCells.drop( [ idx ] )

        self.currentCells = self.currentCells.reset_index(drop=True)
        
        self.updateCanvas1()
        self.setFocus()

    #-----------------------------------------------------------------------------------------------
    # UTILS
    #-----------------------------------------------------------------------------------------------

    def updateRadioBtn(self):
        if self.currentChannel == '488nm':
            self._488nmBtn.setChecked(True)
        elif self.currentChannel == '561nm':
            self._561nmBtn.setChecked(True)
        elif self.currentChannel == 'CoolLED':
            self.CoolLEDBtn.setChecked(True)
        self.setFocus()

    def setBCslidersMinMax(self):
        self.sld1.setMaximum(np.max(self.stacks[self.currentChannel]))
        self.sld1.setMinimum(np.min(self.stacks[self.currentChannel]))
        self.sld2.setMaximum(np.max(self.stacks[self.currentChannel]))
        self.sld2.setMinimum(np.min(self.stacks[self.currentChannel]))

    def resetBC(self):
        self.sld1.setValue(np.min(self.stacks[self.currentChannel]))
        self.sld2.setValue(np.max(self.stacks[self.currentChannel]))
        
    def updateCanvas1(self):

        self.fig1.clf()
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax1 = self.fig1.add_subplot(111)
        self.canvas1.draw()

        if len( self.stacks.keys() ) == 0:
            # if no images are found, leave the canvas empty
            return

        # plot the image
        self.ax1.cla()
        imgplot = self.ax1.imshow( self.stacks[self.currentChannel][self.sl.value()], cmap = 'gray' )
        
        # remove the white borders and plot outline and spline
        self.ax1.autoscale(False)
        self.ax1.axis('Off')
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # cell text on the image
        for idx, cell in self.currentCells.iterrows():

            if cell.Z == self.sl.value():

                self.ax1.text( cell.X, cell.Y + 18, cell.cname, color='red', size='medium', alpha=.8,
                        rotation=0 )
                self.ax1.plot( cell.X, cell.Y, 'o', color='red', alpha = .8, mew = 0 )

        # change brightness and contrast
        self.sld1.setValue(np.min([self.sld1.value(),self.sld2.value()]))
        self.sld2.setValue(np.max([self.sld1.value(),self.sld2.value()]))
        imgplot.set_clim(self.sld1.value(), self.sld2.value())  

        # redraw the canvas
        self.canvas1.draw()
        self.setFocus()

    def updateCanvas2(self):
        
        # plot the image
        self.ax2.cla()
        imgplot = self.ax2.imshow( self.LEDmovie[ self.tp.value() + self.hatchingtidx ], cmap = 'gray' )
        
        # remove the white borders and plot outline and spline
        self.ax2.autoscale(False)
        self.ax2.axis('Off')
        self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # print gonad position
        gonadPos = extract_pos( self.gpDF.ix[ self.gpDF.tidx == self.tp.value() ].squeeze() ) / self.compression
        self.ax2.plot( gonadPos[0], gonadPos[1], 'o', color='red', ms=10, mew=0, alpha=.5, lw = 0 )      

        # print time
        # print(self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'timesRel' ])
        self.ax2.text( 5, 25, '%.2f' % self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'timesRel' ].values[0], color = 'red' )     

        # redraw the canvas
        self.canvas2.draw()
        self.setFocus()

    def checkConsistencyCellNames( self ):
        ### check consistency of cell names
        tp = self.tp.value()

        # if no cells are labeled (currentCells df is empty), remove all labeled cells in the cellPosDF and return
        if len( self.currentCells ) == 0:
            newCellPosDF = update_cell_pos_DF( self.currentCells, self.cellPosDF, tp )
            correctCellNames = True

        if len( self.currentCells ) > 0:

            correctCellNames = check_cell_names( self.currentCells, self.cellNames )

            # if cells are not properly labeled, give a Warning
            if not correctCellNames:
                QtGui.QMessageBox.about(self,'Warning!','There is a mistake in the cell labels!')

            # else, update final cellPosDF and return OK
            else:
                newCellPosDF = update_cell_pos_DF( self.currentCells, self.cellPosDF, tp )
                self.cellPosDF = newCellPosDF

        return correctCellNames

    def changeSpaceTime(self, whatToChange, increment):


        if whatToChange == 'time':

            # before changinf timepoint, print labeled cells and check if they are OK
            print( self.currentCells )

            # if they are OK (and not going to negative times), change timepoint
            if self.checkConsistencyCellNames() and ( self.tp.value() + increment ) >= 0 :
                self.tp.setValue( self.tp.value() + increment )

        if whatToChange == 'space':
            self.sl.setValue( self.sl.value() + increment )


if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication(sys.argv)
    
    gui = GUI()
    app.setStyle( "plastique" )
    app.installEventFilter( gui )
    sys.exit( app.exec_() )
    


