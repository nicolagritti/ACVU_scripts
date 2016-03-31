# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 08:57:15 2015

@author: gritti

NB: this file is also in Miniconda3/Lib/site-packages. This way the file doesn't need to be in the script folder.
I have a copy here just as a backup!

"""

import glob
from tifffile import *
import numpy as np
#import PIL.Image as Image
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
import re
from datetime import datetime
from skimage.filters import threshold_otsu, rank
from skimage import measure
from skimage.morphology import remove_small_objects, disk
from scipy.ndimage import morphology, filters
import pandas as pd
from matplotlib.path import Path

# from uiget_dir import *
#import psutil




##############################################################################
# LOADING STACK, downsize and crop images
##############################################################################

def load_stack(filename):

    '''
    This function reads a (multi)tif file and returns the images in a numpy 
    array using tifffile package from Gohlke.
    '''

    with TiffFile(filename) as f:
        return f.asarray()

def downsizeStack( stack, scaleFactor ):

	if len(stack.shape) == 3:
	    
	    smallstack = []
	    for img in stack:
	        Nbig = img.shape[0]
	        Nsmall = img.shape[0]/scaleFactor
	        smallimg = ( img.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1) ).astype(np.uint16)
	        smallstack.append( smallimg )
	    smallstack = np.array(smallstack)

	    return smallstack

	if len(stack.shape) == 2:

		Nbig = stack.shape[0]
		Nsmall = stack.shape[0]/scaleFactor
		smallimg = ( stack.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1) ).astype(np.uint16)

		return smallimg


def crop_image( imgs, c, size ):

    # print(imgs.shape, c,size)
    dim = imgs.shape

    if len(imgs.shape) == 3:
        cropstack = np.zeros( ( imgs.shape[0], size, size ) )

        cropstack[ : , 
                    -np.min( [ c[1]-size/2, 0 ] ) : size-np.max( [ c[1]+size/2-dim[1]+1, 0 ] ) , 
                    -np.min( [ c[0]-size/2, 0 ] ) : size-np.max( [ c[0]+size/2-dim[2]+1, 0 ] ) ] = imgs[ :,
                    np.max( [ c[1]-size/2, 0 ] ) : np.min( [ c[1]+size/2, dim[1]-1 ] ) , 
                    np.max( [ c[0]-size/2, 0 ] ) : np.min( [ c[0]+size/2, dim[2]-1 ] ) ]
    if len(imgs.shape) == 2:
        cropstack = np.zeros( ( size, size ) )

        cropstack[
                    -np.min( [ c[1]-size/2, 0 ] ) : size-np.max( [ c[1]+size/2-dim[0]+1, 0 ] ) , 
                    -np.min( [ c[0]-size/2, 0 ] ) : size-np.max( [ c[0]+size/2-dim[1]+1, 0 ] ) ] = imgs[
                    np.max( [ c[1]-size/2, 0 ] ) : np.min( [ c[1]+size/2, dim[0]-1 ] ) , 
                    np.max( [ c[0]-size/2, 0 ] ) : np.min( [ c[0]+size/2, dim[1]-1 ] ) ]

    return cropstack



##############################################################################
# UTILS
##############################################################################

def distance( pos1, pos2 ):
    return np.sqrt(np.sum((pos1-pos2)**2))

def save_data_frame( df, path, fileName ):

    pickleName = os.path.join( path, fileName )
    pickle.dump( df, open( pickleName, 'wb' ), protocol=2 )

def load_data_frame( path, fileName ):

    pickleName = os.path.join( path, fileName )
    return pickle.load( open( pickleName, 'rb' ) )

def mag2pxlsize(magnification):
    return 6.5/magnification

def extract_pos( series ):
    return np.array( [ series.X, series.Y ] )

def extract_pos488( series ):
    return np.array( [ series.X488nm, series.Y488nm ] )

def extract_pos561( series ):
    return np.array( [ series.X561nm, series.Y561nm ] )

def extract_488( series ):
    return series._488nm

def extract_561( series ):
    return series._561nm

def extract_3Dpos( series ):
    return np.array( [ series.X, series.Y, series.Z ] )

def extract_out( series ):
    return np.transpose( np.array( [ series.Xout, series.Yout ] ) )

def closer_cell( pos, cells ):

    index = np.nan

    minDist = 100000

    for idx, c in cells.iterrows():
        cpos = extract_3Dpos( c )

        if distance( pos, cpos ) < minDist:

            minDist = distance( pos, cpos )
            index = idx

    return index

def closer_2Dpos( refpos, pos ):

    index = np.nan

    minDist = 100000

    for idx, c in pos.iterrows():
        cpos = extract_pos( c )

        if distance( refpos, cpos ) < minDist:

            minDist = distance( refpos, cpos )
            index = idx

    return index

### functions for parameters and times pickle files

def create_params( path, wormName, magnification, scaleFactor, firsttp):

    worm = {
        'exp': '',
        'worm': '',
        'magnification': np.nan,
        'pxlSize': np.nan,
        'tidxHatch': np.nan,
        'compression': np.nan
        }
        
    df = pd.Series( worm )

    # insert parameters
    df.exp = path.split('\\')[-1]
    df.worm = wormName
    df.magnification = magnification
    df.pxlSize = mag2pxlsize(magnification)
    df.tidxHatch = firsttp
    df.compression = scaleFactor
    
    return df 

def create_times(path, zero=0):

    '''
    Reads all the metadata files in the specified folder and returns a DataFrame
    containing the absolute and relative (to hatching) time indexes, as well as the
    actual times (in minute) relative to hatching time.
    '''

    zero = int(zero)
    flist = glob.glob( os.path.join( path, 'z*.txt' ) )
    flist.sort()
    
    # print(flist,path)
    times = []
    tidx = []
    ftimezero = flist[zero]

    timesDF = { 'fName': [ f.split('\\')[-1][:5] for f in flist ],
              'tidxRel': [np.nan for i in flist],
              'timesRel': [np.nan for i in flist] }

    timesDF = pd.DataFrame( timesDF )

    with open(ftimezero, 'r') as f:

        line = ''
        while 'Date/Time' not in line:
            line = f.readline()

        timezero = datetime.strptime(line.strip().split(': ')[1:][0], '%Y-%m-%d %H:%M:%S')

    for idx, fname in enumerate( flist ):

        with open(fname, 'r') as f:

            line = ''
            while 'Date/Time' not in line:
                line = f.readline()

            date_time = datetime.strptime(line.strip().split(': ')[1:][0], '%Y-%m-%d %H:%M:%S')

            tidx.append( idx-zero )
            times.append( ( date_time - timezero ).total_seconds() / 60. / 60. )

    timesDF.tidxRel = tidx
    timesDF.timesRel = times

    return timesDF

### FUNCTIONS FOR GONAD POSITION

def create_gonad_pos( tDF ):

    gpDF = pd.DataFrame( { 'tidx': tDF.tidxRel,
            'X': np.nan,
            'Y': np.nan } )

    return gpDF

### CELL POSITION FUNCTIONS

def is_position_cell( cell ):

    if not np.isnan( cell.X ):
        return True
    
    return False

def first_tidx_pos_single_cell( cell ):

    cell.isThere = np.isnan( cell.X.values ) == False

    if np.sum( cell.isThere ) > 0:
        firsttidx = np.min( cell.ix[ cell.isThere == True, 'tidx' ] )

    else:
        firsttidx = np.nan

    return firsttidx

def last_tidx_pos_single_cell( cell ):

    cell.isThere = np.isnan( cell.X.values ) == False

    if np.sum( cell.isThere ) > 0:
        lasttidx = np.max( cell.ix[ cell.isThere == True, 'tidx' ] )

    else:
        lasttidx = np.nan

    return lasttidx

def create_cell_pos( tDF, cNames ):

    singleCell = pd.DataFrame( { 'tidx': tDF.tidxRel,
            'X': np.nan,
            'Y': np.nan,
            'Z': np.nan } )

    cellPosDF = { cn: singleCell.copy() for cn in cNames }

    return cellPosDF

def extract_current_cell_pos( dictionary, tp ):

    info = []
    for key in dictionary.keys():
        cell = dictionary[ key ].ix[ dictionary[key].tidx == tp ].squeeze()
        if is_position_cell( cell ):
            info.append([key,cell.X,cell.Y,cell.Z])

    cellsDF = pd.DataFrame( { 'cname': [ row[0] for row in info ],
                'X': [ row[1] for row in info ],
                'Y': [ row[2] for row in info ],
                'Z': [ row[3] for row in info ],
                'tidx': tp  } )

    return cellsDF

def create_single_cell_pos( pos, tp ):
    cell = pd.DataFrame( { 'cname': '?',
                            'X': pos[0],
                            'Y': pos[1],
                            'Z': pos[2],
                            'tidx': tp }, index = [0] )
    return cell

def check_cell_names( currentCells, availableNames ):

    # check if all cells are labeled
    if any( [ i == '?' for i in currentCells.cname ] ):
        return False

    # check if backgrounds are labeled
    for bckg in [i for i in availableNames if i[0]=='b']:
        if bckg not in list( currentCells.cname ):
            return False
    
    for cname in list( currentCells.cname ):
        # check if cells are labeled only once
        if list( currentCells.cname ).count( cname )>1:
            return False
        # check if all cells have available names
        if not cname in availableNames:
            return False

    # check that all the cells belonging to the same lineage have the same length names (no mother/daughter in the same tp)
    for mothers in ['1.','4.','b_']:
        cnames = list( currentCells.ix[ currentCells.cname.str.startswith( mothers ), 'cname' ] )
        if not all(len(x) == len(cnames[0]) for x in cnames):
            return False

    # if everything is correct,return True
    return True

def update_cell_pos_DF( currentCells, cellPosDF, tp ):

    newCellPosDF = cellPosDF.copy()

    for cname in newCellPosDF.keys():
        newCellPosDF[ cname ].ix[ newCellPosDF[ cname ].tidx == tp, 'X' ] = np.nan
        newCellPosDF[ cname ].ix[ newCellPosDF[ cname ].tidx == tp, 'Y' ] = np.nan
        newCellPosDF[ cname ].ix[ newCellPosDF[ cname ].tidx == tp, 'Z' ] = np.nan

    for idx, cell in currentCells.iterrows():
        pos = extract_3Dpos( cell )
        newCellPosDF[ cell.cname ].ix[ newCellPosDF[ cell.cname ].tidx == tp, 'X' ] = pos[0]
        newCellPosDF[ cell.cname ].ix[ newCellPosDF[ cell.cname ].tidx == tp, 'Y' ] = pos[1]
        newCellPosDF[ cell.cname ].ix[ newCellPosDF[ cell.cname ].tidx == tp, 'Z' ] = pos[2]

    return newCellPosDF

### CELL OUTLINE FUNCTIONS

def first_tidx_pos_all_cells( dictionary ):

    tps = []
    
    for cname in dictionary.keys():
        tps.append( first_tidx_pos_single_cell( dictionary[ cname ] ) )

    firsttidx = np.nanmin( tps )

    if np.isnan( firsttidx ):
        firsttidx = 0

    return firsttidx

def last_tidx_pos_all_cells( dictionary ):

    tps = []
    
    for cname in dictionary.keys():
        tps.append( last_tidx_pos_single_cell( dictionary[ cname ] ) )

    lasttidx = np.nanmax( tps )

    if np.isnan( lasttidx ):
        lasttidx = 0

    return lasttidx

def is_outline_cell( cell ):

    if len( cell.Xout ) > 1:
        return True
    
    return False

def create_cell_out( tDF, cNames ):

    singleCell = pd.DataFrame( { 'tidx': tDF.tidxRel,
            'Xout': np.nan,
            'Yout': np.nan } )

    singleCell.Xout = singleCell.Xout.astype(object)
    singleCell.Yout = singleCell.Yout.astype(object)

    for idx, row in singleCell.iterrows():
        singleCell.Xout.values[ idx ] = np.array( [ np.nan ] )
        singleCell.Yout.values[ idx ] = np.array( [ np.nan ] )

    cellOutDF = { cn: singleCell.copy() for cn in cNames }

    return cellOutDF

def extract_current_cell_out( dict_pos, dict_out, tp ):

    cellsDF = extract_current_cell_pos( dict_pos, tp )

    if len(cellsDF) > 0:
        cellsDF.ix[:,'Xout'] = np.nan
        cellsDF.ix[:,'Yout'] = np.nan

        cellsDF.Xout = cellsDF.Xout.astype(object)
        cellsDF.Yout = cellsDF.Yout.astype(object)

        for idx, row in cellsDF.iterrows():

            outline = extract_out( dict_out[ row.cname ].ix[ dict_out[ row.cname ].tidx == tp ].squeeze() )

            cellsDF.Xout.values[ idx ] = outline[:,0]
            cellsDF.Yout.values[ idx ] = outline[:,1]
    # print(cellsDF)

    return cellsDF

def update_cell_out_DF( currentCells, cellOutDF, tp ):

    newCellOutDF = { cn: np.nan for cn in cellOutDF.keys() }

    for cname in cellOutDF.keys():

        newCell = cellOutDF[cname].copy()

        index = newCell.ix[ newCell.tidx == tp ].index

        newCell.Xout.values[index] = [ np.array( [ np.nan ] ) ]
        newCell.Yout.values[index] = [ np.array( [ np.nan ] ) ]

        newCellOutDF[ cname ] = newCell

    for idx, cell in currentCells.iterrows():
        pos = extract_out( cell )

        index = newCellOutDF[ cell.cname ].ix[ newCellOutDF[ cell.cname ].tidx == tp ].index

        newCellOutDF[ cell.cname ].Xout.values[index] = [ pos[:,0] ]
        newCellOutDF[ cell.cname ].Yout.values[index] = [ pos[:,1] ]

    return newCellOutDF

### CELL FLUO FUNCTIONS

def create_cell_fluo( tDF, cNames ):

    singleCell = pd.DataFrame( { 'tidx': tDF.tidxRel,
            '_488nm': np.nan,
            '_561nm': np.nan,
            'X488nm': np.nan,
            'Y488nm': np.nan,
            'X561nm': np.nan,
            'Y561nm': np.nan } )

    cellFluoDF = { cn: singleCell.copy() for cn in cNames }

    return cellFluoDF

### FLUORESCENCE QUANTIFICATION FUNCTIONS

def flat_field_correction( imgs, darkField, flatField, gp ):

    size = imgs[0].shape[0]
    medianCorrection = np.median( ( flatField - darkField ) )

    darkF = crop_image( darkField, gp, size )

    flatF = crop_image( flatField, gp, size )

    num = np.clip( imgs.astype(np.float) - darkF.astype(np.float), 0, None )
    den = np.clip( flatF.astype(np.float) - darkF.astype(np.float), 0, None )

    return ( medianCorrection * num / den ).astype(np.uint16)

def extract_current_cell_fluo( dict_pos, dict_out, dict_fluo, tp ):

    cellsDF = extract_current_cell_out( dict_pos, dict_out, tp )

    if len(cellsDF) > 0:
        cellsDF.ix[:,'_488nm'] = np.nan
        cellsDF.ix[:,'_561nm'] = np.nan
        cellsDF.ix[:,'X488nm'] = np.nan
        cellsDF.ix[:,'X561nm'] = np.nan
        cellsDF.ix[:,'Y488nm'] = np.nan
        cellsDF.ix[:,'Y561nm'] = np.nan

        for cname in dict_fluo.keys():
            cell = dict_fluo[ cname ].ix[ dict_fluo[ cname ].tidx == tp ].squeeze()

            pos488 = extract_pos488( cell )
            pos561 = extract_pos561( cell )
            _488 = extract_488( cell )
            _561 = extract_561( cell )

            cellsDF.ix[ cellsDF.cname == cname, '_488nm' ] = _488
            cellsDF.ix[ cellsDF.cname == cname, '_561nm' ] = _561
            cellsDF.ix[ cellsDF.cname == cname, 'X488nm' ] = pos488[0]
            cellsDF.ix[ cellsDF.cname == cname, 'X561nm' ] = pos561[0]
            cellsDF.ix[ cellsDF.cname == cname, 'Y488nm' ] = pos488[1]
            cellsDF.ix[ cellsDF.cname == cname, 'Y561nm' ] = pos561[1]

    return cellsDF

def update_cell_fluo_DF( currentCells, cellFluoDF, tp ):

    newCellFluoDF = { cn: np.nan for cn in cellFluoDF.keys() }

    for cname in cellFluoDF.keys():

        newCell = cellFluoDF[cname].copy()

        index = newCell.ix[ newCell.tidx == tp ].index

        newCell._488nm.values[index] = np.nan
        newCell._561nm.values[index] = np.nan
        newCell.X488nm.values[index] = np.nan
        newCell.Y488nm.values[index] = np.nan
        newCell.X561nm.values[index] = np.nan
        newCell.Y561nm.values[index] = np.nan

        newCellFluoDF[ cname ] = newCell

    for idx, cell in currentCells.iterrows():
        pos488 = extract_pos488( cell )
        pos561 = extract_pos561( cell )
        _488 = extract_488( cell )
        _561 = extract_561( cell )

        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, 'X488nm' ] = pos488[0]
        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, 'Y488nm' ] = pos488[1]
        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, 'X561nm' ] = pos561[0]
        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, 'Y561nm' ] = pos561[1]
        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, '_488nm' ] = _488
        newCellFluoDF[ cell.cname ].ix[ newCellFluoDF[ cell.cname ].tidx == tp, '_561nm' ] = _561

    return newCellFluoDF

def calculate_fluo_intensity( img, center, outline ):

    imgpxl = img.shape[0]

    vertices = np.array( [ center[0] + np.append(outline[:,0],outline[0,0]), center[1] + np.append(outline[:,1],outline[0,1]) ] ).T

    p = Path(vertices)
    
    # create the mask (image full of 0 and 1, the 1s are wherethe cell is)
    points = [ (i,j) for i in np.arange(imgpxl) for j in np.arange(imgpxl) ]
    mask = p.contains_points(points).reshape(imgpxl,imgpxl).T

    # plt.figure()
    # plt.imshow(mask)
    # plt.show()
    # plt.figure()
    # plt.imshow(img)
    # plt.show()

    return np.sum( mask * img ) / np.sum( mask )

def calculate_fluo_intensity_bckg( img, imgpxl, center ):

    img = img[ center[1]-imgpxl/2:center[1]+imgpxl/2, center[0]-imgpxl/2:center[0]+imgpxl/2 ]
    size = len( img.flatten() )

    return np.sum( img ) / size

def update_current_cell_fluo( currentCells, cell, channel, drift, signal ):

    newCurrentCells = currentCells.copy()

    newCurrentCells.ix[ newCurrentCells.cname == cell.cname, '_' + channel ] = signal
    newCurrentCells.ix[ newCurrentCells.cname == cell.cname, 'X' + channel ] = int( cell.X + drift[0] )
    newCurrentCells.ix[ newCurrentCells.cname == cell.cname, 'Y' + channel ] = int( cell.Y + drift[1] )

    return newCurrentCells


def find_interesting_cells( currentCells ):
    cnames = list( currentCells.cname )

    if '1.p' in cnames:
        cell1 = currentCells.ix[ currentCells.cname == '1.p' ]
        if '4.a' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.a' ]
        if '4.aa' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.aa' ]

    if '1.pp' in cnames:
        cell1 = currentCells.ix[ currentCells.cname == '1.pp' ]
        if '4.a' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.a' ]
        if '4.aa' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.aa' ]
        if '4.aaa' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.aaa' ]

    if '1.ppp' in cnames:
        cell1 = currentCells.ix[ currentCells.cname == '1.ppp' ]
        if '4.aa' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.aa' ]
        if '4.aaa' in cnames:
            cell2 = currentCells.ix[ currentCells.cname == '4.aaa' ]

    cell1 = cell1.squeeze()
    cell2 = cell2.squeeze()

    bckg1 = currentCells.ix[ currentCells.cname == 'b_'+cell1.cname[0] ].squeeze()
    bckg2 = currentCells.ix[ currentCells.cname == 'b_'+cell2.cname[0] ].squeeze()

    # print(cell1,cell2)
    return ( cell1, cell2, bckg1, bckg2 )

### FUNCTIONS FOR AP-DV IDENTIFICATION

def create_apdv_pos( tDF ):

    singlePos = pd.DataFrame( { 'tidx': tDF.tidxRel,
            'X': np.nan,
            'Y': np.nan } )

    apdvPosDF = { cn: singlePos.copy() for cn in [ 'a', 'p', 'd' ] }

    return apdvPosDF

def extract_current_apdv_pos( dictionary, tp ):

    info = []
    for key in dictionary.keys():
        pos = dictionary[ key ].ix[ dictionary[key].tidx == tp ].squeeze()
        if is_position_cell( pos ):
            info.append([key,pos.X,pos.Y])

    apdvDF = pd.DataFrame( { 'pname': [ row[0] for row in info ],
                'X': [ row[1] for row in info ],
                'Y': [ row[2] for row in info ],
                'tidx': tp  } )

    return apdvDF

def create_single_apdv_pos( pos, tp ):
    pos = pd.DataFrame( { 'pname': '?',
                            'X': pos[0],
                            'Y': pos[1],
                            'tidx': tp }, index = [0] )
    return pos

def update_apdv_pos_DF( currentPos, posDF, tp ):

    newPosDF = posDF.copy()

    for key in newPosDF.keys():
        newPosDF[ key ].ix[ newPosDF[ key ].tidx == tp, 'X' ] = np.nan
        newPosDF[ key ].ix[ newPosDF[ key ].tidx == tp, 'Y' ] = np.nan

    for idx, pos in currentPos.iterrows():
        p = extract_pos( pos )
        newPosDF[ pos.pname ].ix[ newPosDF[ pos.pname ].tidx == tp, 'X' ] = p[0]
        newPosDF[ pos.pname ].ix[ newPosDF[ pos.pname ].tidx == tp, 'Y' ] = p[1]

    return newPosDF
