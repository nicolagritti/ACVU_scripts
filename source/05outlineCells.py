

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
        self.imgpxl = 50
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

        self.initializeCanvas1()
        self.initializeCanvas2()
        self.show()
        
        # BIND BUTTONS TO FUNCTIONS
        
        loadBtn.clicked.connect(self.selectWorm)
        saveBtn.clicked.connect(self.saveData)

        self.tp.valueChanged.connect(self.loadNewStack)
        self.sl.valueChanged.connect(self.updateAllCanvas)
        self.sld1.valueChanged.connect(self.changeBC)
        self.sld2.valueChanged.connect(self.changeBC)

        self._488nmBtn.toggled.connect(self.radioClicked)
        self._561nmBtn.toggled.connect(self.radioClicked)
        self.CoolLEDBtn.toggled.connect(self.radioClicked)

        self.fig1.canvas.mpl_connect('button_press_event',self.onMouseClickOnCanvas1)        
        self.fig1.canvas.mpl_connect('scroll_event',self.wheelEvent)        

    #-----------------------------------------------------------------------------------------------
    # FORMATTING THE WINDOW
    #-----------------------------------------------------------------------------------------------

    def HLine(self):
        
        toto = QtGui.QFrame()
        toto.setFrameShape(QtGui.QFrame.HLine)
        toto.setFrameShadow(QtGui.QFrame.Sunken)
        return toto

    #-----------------------------------------------------------------------------------------------
    # BUTTON FUNCTIONS
    #-----------------------------------------------------------------------------------------------

    def selectWorm(self):

        ### store the folders
        self.pathDial = QtGui.QFileDialog.getExistingDirectory(self, 'Select a folder', 'Y:\\Images')
        self.worm = self.pathDial.split("\\")[-1].split('_')[0]
        self.path = os.path.dirname( self.pathDial )
        self.setWindowTitle('Outline Cells - ' + self.pathDial)
        
        ### give error message if there is no CoolLED movie in the selected folder
        if not os.path.isfile( os.path.join( self.pathDial, 'CoolLED_movie.tif' ) ):
            QtGui.QMessageBox.about(self,'Warning!','There is no movie in this folder! Create a movie first!')
            return

        # detect available channels
        self.channels = []
        chns = ['CoolLED','488nm','561nm']
        for c in chns:

            if os.path.isfile( os.path.join( self.pathDial, c + '_movie.tif' ) ):

                self.channels.append(c)

        self.currentChannel = self.channels[0]

        ### load parameters and times dataframes
        self.paramsDF = load_data_frame( self.path, self.worm + '_01params.pickle' )
        self.timesDF = load_data_frame( self.path, self.worm + '_01times.pickle' )
        self.gpDF = load_data_frame( self.path, self.worm + '_02gonadPos.pickle' )
        self.cellPosDF = load_data_frame( self.path, self.worm + '_04cellPos.pickle' )

        # extract some info
        self.compression = self.paramsDF.compression
        self.hatchingtidx = int( self.paramsDF.tidxHatch )

        ### if the cellOutline pickle file already exists, load it, otherwise create a blank one
        if os.path.isfile( os.path.join(self.path, self.worm + '_05cellOut.pickle' ) ):
            self.cellOutDF = load_data_frame( self.path, self.worm + '_05cellOut.pickle' )
        
        else:
            self.cellOutDF = create_cell_out( self.timesDF, self.cellNames )

        self.analyzedCell = '---'

        ### set the timepoint to the hatching time
        self.tp.setMinimum(np.min(self.timesDF.tidxRel))
        self.tp.setMaximum(np.max(self.timesDF.tidxRel))
        self.tp.setValue( first_tidx_pos_all_cells( self.cellPosDF ) )
        
        # self.pathDial.show()
        self.updateAllCanvas()
        self.setFocus()

    def loadNewStack(self):
        
        # print(self.fList['gfp'][self.tp.value()])
        tRow = self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value() ].squeeze()

        print( 'Loading... ', self.pathDial, tRow.fName )

        ### update the text of the fileName
        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])

        ### extract current cells already labeled
        self.currentCells = extract_current_cell_out( self.cellPosDF, self.cellOutDF, self.tp.value() )

        if len(self.currentCells) > 0:
            ### update current analyzed cell
            if self.analyzedCell not in list( self.currentCells.cname ):
                self.analyzedCell = self.currentCells.cname[0]

        ### update the text of the fileName
        self.fName.setText( self.timesDF.ix[ self.timesDF.tidxRel == self.tp.value(), 'fName' ].values[0])

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

        if len( self.currentCells ) > 0:
            ### update slice
            self.sl.setValue( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell, 'Z' ] )
        
        # self.updateTable()
        self.updateAllCanvas()

    def saveData(self):
        
        save_data_frame( self.cellOutDF, self.path, self.worm + '_05cellOut.pickle' )       
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

        pos = np.array( [ int(event.xdata), int(event.ydata) ] )

        outline = extract_out( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

        if event.button == 1:

            if np.isnan( outline[0,0] ):
                outline = np.array( [ pos ] )
            else:
                outline = np.vstack( [ outline, pos ] )

        elif event.button == 3:

            if len( outline[ :, 0 ] ) == 1:
                outline = np.array( [ [ np.nan, np.nan ] ] )
            else:
                outline = outline[:-1,:]

        idx = self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].index

        self.currentCells.Xout.values[ idx ] = [ outline[:,0] ]
        self.currentCells.Yout.values[ idx ] = [ outline[:,1] ]

        self.updateCanvas1()
        self.setFocus()
        # print(event.button,event.xdata,event.ydata)

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

    def initializeCanvas1(self):

        self.fig1.clf()
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax1 = self.fig1.add_subplot(111)
        self.canvas1.draw()

        # plot the image
        self.ax1.cla()
        self.imgplot1 = self.ax1.imshow( np.zeros((self.imgpxl,self.imgpxl)), cmap = 'gray', interpolation = 'nearest' )
        
        # remove the white borders
        self.ax1.autoscale(False)
        self.ax1.axis('Off')
        self.fig1.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # plot cell pos and name
        self.text1 = []
        self.outline1 = []

        # redraw the canvas
        self.canvas1.draw()
        self.setFocus()

    def updateCanvas1(self):

        if ( len( self.stacks.keys() ) == 0 ) or ( len( self.currentCells ) == 0 ):
            # if no images are found, leave the canvas empty
            self.initializeCanvas1()
            return

        # extract current cell data
        pos = extract_3Dpos( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

        # replot the image
        imgpxl = self.imgpxl
        self.imgplot1.set_data( self.stacks[self.currentChannel][self.sl.value(),pos[1]-imgpxl/2:pos[1]+imgpxl/2+1,pos[0]-imgpxl/2:pos[0]+imgpxl/2+1] )
        
        # clear cell text and points
        for text in self.text1:
            text.remove()
        self.text1 = []

        for outline in self.outline1:
            self.ax1.lines.remove(outline)
        self.outline1 = []

        # print cell name
        if pos[2] == self.sl.value():
            self.text1.append( self.ax1.text( 1, 2, self.analyzedCell, color='yellow', size='medium', alpha=.8,
                    rotation=0, fontsize = 20 ) )

        ### draw outline
        outline = extract_out( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

        # print(self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze())

        if len( outline ) > 1:
            outline = np.vstack( [ outline, outline[0] ] )

        self.outline1.append( self.ax1.plot( outline[:,0], outline[:,1], '-x', color='yellow', ms=2, alpha=1., lw = .5 )[0] )

        # change brightness and contrast
        self.sld1.setValue(np.min([self.sld1.value(),self.sld2.value()]))
        self.sld2.setValue(np.max([self.sld1.value(),self.sld2.value()]))

        # # redraw the canvas
        self.canvas1.draw()
        self.setFocus()

    def changeBC(self):

        self.imgplot1.set_clim(self.sld1.value(), self.sld2.value())  

        self.imgplot2.set_clim(self.sld1.value(), self.sld2.value()) 

    def initializeCanvas2(self):

        self.fig2.clf()
        self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        self.ax2 = self.fig2.add_subplot(111)
        self.canvas2.draw()

        # plot the image
        self.ax2.cla()
        size = 2048 / self.compression
        self.imgplot2 = self.ax2.imshow( np.zeros((size,size)), cmap = 'gray' )
        
        # remove the white borders
        self.ax2.autoscale(False)
        self.ax2.axis('Off')
        self.fig2.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        # plot cell pos and name
        self.text2 = []
        self.points2 = []

        # redraw the canvas
        self.canvas2.draw()
        self.setFocus()

    def updateCanvas2(self):
        
        # plot the image
        self.imgplot2.set_data( self.stacks[self.currentChannel][self.sl.value()] )
        
        # clear cell text and points
        for text in self.text2:
            text.remove()
        self.text2 = []

        for points in self.points2:
            self.ax2.lines.remove(points)
        self.points2 = []

        # draw cell text and point
        if len( self.currentCells ) > 0:
            pos = extract_3Dpos( self.currentCells.ix[ self.currentCells.cname == self.analyzedCell ].squeeze() )

            for idx, cell in self.currentCells.iterrows():

                if cell.Z == self.sl.value():

                    color = 'red'
                    if cell.cname == self.analyzedCell:
                        color = 'yellow'

                    self.text2.append( self.ax2.text( cell.X, cell.Y + 18, cell.cname, color=color, size='medium', alpha=.8,
                            rotation=0 ) )
                    self.points2.append( self.ax2.plot( cell.X, cell.Y, 'o', color=color, alpha = .8, mew = 0 )[0] )

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

if __name__ == '__main__':
    
    app = QtGui.QApplication.instance() # checks if QApplication already exists 
    if not app: # create QApplication if it doesnt exist 
        app = QtGui.QApplication(sys.argv)
    
    gui = GUI()
    app.setStyle("plastique")
    # app.installEventFilter(gui)
    sys.exit(app.exec_())
    


