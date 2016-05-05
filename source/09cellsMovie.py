# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:00:56 2015

@author: gritti


NB: CHANGE outpath ACCORDING TO THE COMPUTER YOU ARE USING!
"""

import glob
from tifffile import *
from generalFunctions import *
import numpy as np
import PIL
from PIL import Image, ImageDraw, ImageFont
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
from scipy.ndimage import interpolation
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def isDaughter(cname, cdata):
    return ( ( cdata.cname.str.startswith( cname ) ) & ( cdata.cname.str.len() == (len(cname)+1) ) ).any()

def coordTransformO2C(r,size):
    return np.array( [ size/2+r[0], size/2-r[1] ] )

def coordTransformC2O(r,size):
    return np.array( [ r[0]-size/2, size/2-r[1] ] )

def rotateCoords(r, alpha):
    return np.array( [ r[0]*np.cos(alpha) - r[1]*np.sin(alpha), r[0]*np.sin(alpha) + r[1]*np.cos(alpha) ] )

def Yreflect(r):
	return np.array( [ r[0], -r[1] ] )

def inverse(pos1,pos2,posn):
	sign = (pos2[0]-pos1[0])*(posn[1]-pos1[1])-(pos2[1]-pos1[0])*(posn[0]-pos1[0])
	return np.sign(sign) == 1

def extractCoords(cell):
	# print(cell)
	# return np.array([cell.X,cell.Y])
	if np.isnan( cell.X488nm ):
		cell.X488nm = cell.X
	if np.isnan( cell.Y488nm ):
		cell.Y488nm = cell.Y
	return np.array([cell.X488nm,cell.Y488nm])

def extractCoordsRot(cell):
	return np.array([cell.cXposRot,cell.cYposRot])

def extractCoordsO(celltp,cname):
	return np.array([celltp.ix[celltp.cname==cname,'cXposO'].values[0],celltp.ix[celltp.cname==cname,'cYposO'].values[0]])

def extractRefPointsCoordsO(postp,pname):
	return np.array([postp.ix[postp.pname==pname,'cXposO'].values[0],postp.ix[postp.pname==pname,'cYposO'].values[0]])

def extractCoordsORot(celltp,cname):
	return np.array([celltp.ix[celltp.cname==cname,'cXposORot'].values[0],celltp.ix[celltp.cname==cname,'cYposORot'].values[0]])

def extractRefPointsCoordsORot(postp,pname):
	return np.array([postp.ix[postp.pname==pname,'cXposORot'].values[0],postp.ix[postp.pname==pname,'cYposORot'].values[0]])

def flatFieldCorrection( imgs, darkField, flatField, bodyData, tidx ):
	size = imgs[0].shape[0]
	medianCorrection = np.median( ( flatField - darkField ) )

	gp = np.floor( bodyData.ix[ bodyData.tidx == tidx, 'gonadPos' ].values[0] )

	darkF = np.zeros( ( size, size ) )
	flatF = np.zeros( ( size, size ) )

	darkF[ -np.min( [ gp[1]-size/2, 0 ] ) : size-np.max( [ gp[1]+size/2-2047, 0 ] ) , 
							-np.min( [ gp[0]-size/2, 0 ] ) : size-np.max( [ gp[0]+size/2-2047, 0 ] ) ] = darkField[
							np.max( [ gp[1]-size/2, 0 ] ) : np.min( [ gp[1]+size/2, 2047 ] ) , 
							np.max( [ gp[0]-size/2, 0 ] ) : np.min( [ gp[0]+size/2, 2047 ] ) ]
	flatF[ -np.min( [ gp[1]-size/2, 0 ] ) : size-np.max( [ gp[1]+size/2-2047, 0 ] ) , 
							-np.min( [ gp[0]-size/2, 0 ] ) : size-np.max( [ gp[0]+size/2-2047, 0 ] ) ] = flatField[
							np.max( [ gp[1]-size/2, 0 ] ) : np.min( [ gp[1]+size/2, 2047 ] ) , 
							np.max( [ gp[0]-size/2, 0 ] ) : np.min( [ gp[0]+size/2, 2047 ] ) ]

	# print(flatF,medianCorrection)

	# plt.figure()
	# plt.imshow(darkF,cmap='gray')
	# plt.figure()
	# plt.imshow(flatF,cmap='gray')

	return ((imgs - darkF)/((flatF-darkF)/medianCorrection)).astype(np.uint16)

def rotateImage(imgs, alpha):
	return np.array( [ interpolation.rotate(img,alpha*180/np.pi,reshape=False) for img in imgs ] )

def flipImage(imgs,refPoints):
	if refPoints.ix[refPoints.pname=='d','cYposRot'].values[0] > refPoints.ix[refPoints.pname=='a','cYposRot'].values[0]:
		return np.array( [ img[::-1, :] for img in imgsRot ] ), True
	else:
		return imgs, False

# path = 'X:\\Simone\\16_01_25_MCHERRY_LAG2YFP\\'
path = 'X:\\Simone\\160420_LIN12_GFP_hist_mCherry'

worms = ['C01']
channel = '561nm'

size=40

### INCLUDE DRIFT CORRECTION

for w in worms:
	# w = 'C01'

	### Nicola computer
	outpath = os.path.join( path, w + '_analyzedImages' )
	# ### Simone computer
	# outpath = 'Y:\\Presentation\\Figures\\160211_ACVU\\movie_' + path.split('\\')[-2] + '_%s' %w

	print(outpath)

	### load parameters and times dataframes
	paramsDF = load_data_frame( path, w + '_01params.pickle' )
	timesDF = load_data_frame( path, w + '_01times.pickle' )
	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, w + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
	apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )
	
	### load darkField and flatField
	darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )
	flatField = load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )
	
	### create the movies
	movie = [[],[]]

	# for idx, trow in timesDF.ix[timesDF.tidxRel==42].iterrows(): # test single timepoints
	for idx, trow in timesDF.iterrows():

		### skip problematic timepoints
		# if trow.tidxRel == 42:
		# 	continue

		currentCells = extract_current_cell_fluo( cellPosDF, cellOutDF, cellFluoDF, trow.tidxRel )
		currentCellsOriginalPos = extract_current_cell_pos( cellPosDF, trow.tidxRel )

		### if there are only 2 cells (one is the bckg) then skip the timepoint
		if len( currentCells.cname ) > 2:

			### find out which cells we want to show
			(c1,c2,b1,b2) = find_interesting_cells( currentCells )
			interestingCells = pd.concat([c1, c2],axis=1).transpose().reset_index()
			refPoints = extract_current_apdv_pos( apdvPosDF, trow.tidxRel )
			print(trow.tidxRel,list(interestingCells.cname))

			### load stack
			imgs = load_stack( os.path.join( outpath, trow.fName + channel + '.tif') )
			# print(imgfile)
			imgSize = imgs[0].shape[0]

			### correct for XY
			gonadPos = extract_pos( gpDF.ix[ gpDF.tidx == trow.tidxRel ].squeeze() )
			imgsXYCorr = flat_field_correction( imgs, darkField, flatField, gonadPos )
			# print(imgsXYCorr[22])

			# ### plot raw data
			# plt.figure()
			# plt.xlim(0,imgs[21].shape[0])
			# plt.ylim(0,imgs[21].shape[1])
			# plt.gca().invert_yaxis()
			# plt.imshow(imgs[21],cmap='gray')
			# colors = ['r','g']
			# for kdx, cell in interestingCells.iterrows():
			# 	pos = extractCoords(cell)
			# 	plt.plot( pos[0], pos[1], 'o',color = colors[kdx])
			# for kdx, cell in refPoints.iterrows():
			# 	pos = extract_pos(cell) - gonadPos + 256
			# 	plt.plot( pos[0], pos[1], 'o',color = 'b')

			### get coordinates with respect to center of the image (origin)
			for jdx, cell in interestingCells.iterrows():

				posC = extractCoords(cell)
				posO = coordTransformC2O(posC, imgSize)

				interestingCells.ix[ interestingCells.cname == cell.cname, 'cXposO' ] = posO[0]
				interestingCells.ix[ interestingCells.cname == cell.cname, 'cYposO' ] = posO[1]

			for jdx, pos in refPoints.iterrows():
				# print(pos)
				posC = extract_pos(pos) - gonadPos + 256
				posO = coordTransformC2O(posC, imgSize)

				refPoints.ix[ refPoints.pname == pos.pname, 'cXposO' ] = posO[0]
				refPoints.ix[ refPoints.pname == pos.pname, 'cYposO' ] = posO[1]

			### rotate image

			vect = extractRefPointsCoordsO(refPoints,'p') - extractRefPointsCoordsO(refPoints,'a')
			alpha = -np.arctan2(vect[1],vect[0])
			# print(alpha)
			imgsRot = rotateImage(imgsXYCorr, alpha)

			### rotate coordinates
			for jdx, cell in interestingCells.iterrows():

				posO = extractCoordsO(interestingCells, cell.cname)
				posORot = rotateCoords(posO,alpha)

				interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
				interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = posORot[1]

				interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
				interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]

			for jdx, pos in refPoints.iterrows():

				posO = extractRefPointsCoordsO(refPoints,pos.pname)
				posORot = rotateCoords(posO,alpha)

				refPoints.ix[refPoints.pname==pos.pname,'cXposORot'] = posORot[0]
				refPoints.ix[refPoints.pname==pos.pname,'cYposORot'] = posORot[1]

				refPoints.ix[refPoints.pname==pos.pname,'cXposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[0]
				refPoints.ix[refPoints.pname==pos.pname,'cYposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[1]

			### flip image
			imgsFinal, flipTF = flipImage( imgsRot, refPoints )

			### flip coordinates
			if flipTF:
				# print('flipping...')
				for jdx, cell in interestingCells.iterrows():
					posORot = extractCoordsORot(interestingCells, cell.cname)

					interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
					interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = -posORot[1]

					interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
					interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]

			# ### plot rotated and flipped image
			# plt.figure()
			# plt.xlim(0,imgsFinal[21].shape[0])
			# plt.ylim(0,imgsFinal[21].shape[1])
			# plt.gca().invert_yaxis()
			# plt.imshow(imgsFinal[21],cmap='gray')
			# colors = ['r','g']
			# for kdx, cell in interestingCells.iterrows():
			# 	pos = extractCoordsRot(cell)
			# 	plt.plot( pos[0], pos[1], 'o',color = colors[kdx])
			# for kdx, cell in refPoints.iterrows():
			# 	pos = extractCoordsRot(cell)
			# 	plt.plot( pos[0], pos[1], 'o',color = 'b')

			### crop images
			for jdx, cell in interestingCells.iterrows():

				pos = extractCoordsRot(cell)
				small = imgsFinal[ currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cell.cname,'Z'], pos[1]-size:pos[1]+size, pos[0]-size:pos[0]+size ]

				# print(cell, small.shape)
				# plt.figure()
				# plt.imshow(small,cmap='gray')

				movie[jdx].append(small)
				# print(np.array(movie[jdx]).shape)

	# print(np.array(movie[0]).shape)
	## save movies
	imsave( os.path.join( outpath, 'cell_'+channel+'_1.tif' ), np.array(movie[0]) )
	imsave( os.path.join( outpath, 'cell_'+channel+'_4.tif' ), np.array(movie[1]) )

plt.show()

### DO NOT INCLUDE DRIFT CORRECTION

# for w in worms:
# 	# w = 'C01'

# 	### Nicola computer
# 	outpath = os.path.join( path, w + '_analyzedImages' )
# 	# ### Simone computer
# 	# outpath = 'Y:\\Presentation\\Figures\\160211_ACVU\\movie_' + path.split('\\')[-2] + '_%s' %w

# 	print(outpath)

# 	### load parameters and times dataframes
# 	paramsDF = load_data_frame( path, w + '_01params.pickle' )
# 	timesDF = load_data_frame( path, w + '_01times.pickle' )
# 	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
# 	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
# 	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
# 	apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )
	
# 	### load darkField and flatField
# 	darkField = load_stack( 'X:\\Orca_calibration\\AVG_darkField.tif' )
# 	flatField = load_stack( os.path.join( path, 'AVG_flatField_'+channel+'.tif' ) )
	
# 	### create the movies
# 	movie = [[],[]]

# 	# for idx, trow in timesDF.ix[timesDF.tidxRel<60].iterrows():
# 	for idx, trow in timesDF.iterrows():
# 		# tidx = 37 ### example
# 		currentCells = extract_current_cell_pos( cellPosDF, trow.tidxRel )
# 		currentCellsOriginalPos = extract_current_cell_pos( cellPosDF, trow.tidxRel )

# 		### if there are only 2 cells (one is the bckg) then skip the timepoint
# 		if len( currentCells.cname ) > 2:

# 			### find out which cells we want to show
# 			(c1,c2,b1,b2) = find_interesting_cells( currentCells )
# 			interestingCells = pd.concat([c1, c2],axis=1).transpose().reset_index()
# 			refPoints = extract_current_apdv_pos( apdvPosDF, trow.tidxRel )
# 			print(trow.tidxRel,list(interestingCells.cname))

# 			### load stack
# 			imgs = load_stack( os.path.join( outpath, trow.fName + channel + '.tif') )
# 			# print(imgfile)
# 			imgSize = imgs[0].shape[0]

# 			# ### plot raw data
# 			# plt.figure()
# 			# plt.xlim(0,imgs[22].shape[0])
# 			# plt.ylim(0,imgs[22].shape[1])
# 			# plt.gca().invert_yaxis()
# 			# plt.imshow(imgs[22],cmap='gray')
# 			# for cname in [c1,c2,'n.']:
# 			# 	cell = cdatatp.ix[ cdatatp.cname == cname ].squeeze()
# 			# 	pos = extractCoords(cell)
# 			# 	plt.plot( pos[0], pos[1], 'o')
# 			# print(imgs[22])

# 			### correct for XY
# 			gonadPos = extract_pos( gpDF.ix[ gpDF.tidx == trow.tidxRel ].squeeze() )
# 			imgsXYCorr = flat_field_correction( imgs, darkField, flatField, gonadPos )
# 			# print(imgsXYCorr[22])

# 			### get coordinates with respect to center of the image (origin)
# 			for jdx, cell in interestingCells.iterrows():

# 				posC = extractCoords(cell)
# 				posO = coordTransformC2O(posC, imgSize)

# 				interestingCells.ix[ interestingCells.cname == cell.cname, 'cXposO' ] = posO[0]
# 				interestingCells.ix[ interestingCells.cname == cell.cname, 'cYposO' ] = posO[1]

# 			for jdx, pos in refPoints.iterrows():

# 				posC = extract_pos(pos) - gonadPos + 256
# 				posO = coordTransformC2O(posC, imgSize)

# 				refPoints.ix[ refPoints.pname == pos.pname, 'cXposO' ] = posO[0]
# 				refPoints.ix[ refPoints.pname == pos.pname, 'cYposO' ] = posO[1]

# 			### rotate image
# 			vect = extractRefPointsCoordsO(refPoints,'a') - extractRefPointsCoordsO(refPoints,'p')
# 			alpha = -np.arctan2(vect[1],vect[0])
# 			imgsRot = rotateImage(imgsXYCorr, alpha)

# 			### rotate coordinates
# 			for jdx, cell in interestingCells.iterrows():

# 				posO = extractCoordsO(interestingCells, cell.cname)
# 				posORot = rotateCoords(posO,alpha)

# 				interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
# 				interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = posORot[1]

# 				interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
# 				interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]

# 			for jdx, pos in refPoints.iterrows():

# 				posO = extractRefPointsCoordsO(refPoints,pos.pname)
# 				posORot = rotateCoords(posO,alpha)

# 				refPoints.ix[refPoints.pname==pos.pname,'cXposORot'] = posORot[0]
# 				refPoints.ix[refPoints.pname==pos.pname,'cYposORot'] = posORot[1]

# 				refPoints.ix[refPoints.pname==pos.pname,'cXposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[0]
# 				refPoints.ix[refPoints.pname==pos.pname,'cYposRot'] = coordTransformO2C( extractRefPointsCoordsORot(refPoints,pos.pname), imgSize )[1]

# 			### flip image
# 			imgsFinal, flipTF = flipImage( imgsRot, refPoints )

# 			### flip coordinates
# 			if flipTF:
# 				for jdx, cell in interestingCells.iterrows():
# 					posORot = extractCoordsORot(interestingCells, cell.cname)

# 					interestingCells.ix[interestingCells.cname==cell.cname,'cXposORot'] = posORot[0]
# 					interestingCells.ix[interestingCells.cname==cell.cname,'cYposORot'] = -posORot[1]

# 					interestingCells.ix[interestingCells.cname==cell.cname,'cXposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[0]
# 					interestingCells.ix[interestingCells.cname==cell.cname,'cYposRot'] = coordTransformO2C( extractCoordsORot(interestingCells,cell.cname), imgSize )[1]

# 			# ### plot rotated and flipped image
# 			# plt.figure()
# 			# plt.xlim(0,imgsFinal[22].shape[0])
# 			# plt.ylim(0,imgsFinal[22].shape[1])
# 			# plt.gca().invert_yaxis()
# 			# plt.imshow(imgsFinal[22],cmap='gray')
# 			# for cname in [c1,c2,'n.']:
# 			# 	cell = cdatatp.ix[ cdatatp.cname == cname ].squeeze()
# 			# 	pos = extractCoordsRot(cell)
# 			# 	plt.plot( pos[0], pos[1], 'o')

# 			### crop images
# 			for jdx, cell in interestingCells.iterrows():

# 				pos = extractCoordsRot(cell)
# 				small = imgsFinal[ currentCellsOriginalPos.ix[currentCellsOriginalPos.cname == cell.cname,'Z'], pos[1]-25:pos[1]+25, pos[0]-25:pos[0]+25 ]
# 				# print(small.shape)
# 				# plt.figure()
# 				# plt.imshow(small,cmap='gray')
# 				movie[jdx].append(small)

# 	print(np.array(movie[0]).shape)
# 	### save movies
# 	imsave( os.path.join( outpath, 'cell_'+channel+'_1_noshifted.tif' ), np.array(movie[0]) )
# 	imsave( os.path.join( outpath, 'cell_'+channel+'_4_noshifted.tif' ), np.array(movie[1]) )

# plt.show()
