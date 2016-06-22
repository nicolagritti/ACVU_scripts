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
from matplotlib.colors import LinearSegmentedColormap
mpl.rcParams['pdf.fonttype'] = 42

# path = 'X:\\Simone\\16_01_25_MCHERRY_LAG2YFP\\'
path = 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry'
worms = ['C02','C05','C07','C08','C11','C12','C15','C19','C20','C21','C25','C27','C28','C33']

path = 'X:\\Simone\\160516_lag2_YFP_hist_mCherry'
worms = ['C02','C07','C09','C14','C16','C17','C19','C20','C21','C22','C24','C25']
firsttidx = [65, 84,   85,   56,   78,   87,   92,   25,   44,   70,   42,   32]
channel = '488nm'


def plotPhaseFluorescence( path, worm, ax1, color = 'k', channel = '488nm' ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	cellPosDF = load_data_frame( path, worm + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, worm + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']]
	colors = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
				np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )

	### find ecdysis timepoint

	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >= 0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - mintp
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]

	### set relative time

	# # relative to L2 start
	# tpRel = tpL2

	# # relative to hatching time
	# tpRel=0

	# relative to first cell division
	tdiv1 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].timesRel )
	tdiv2 = np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].timesRel )
	tpRel = np.min([tdiv1,tdiv2])

	print(tpRel)

	### plot the timeseries
	for idx, lin in enumerate( lineages ):

		for jdx, key in enumerate( lin ):

			cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
			bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]
			times = timesDF.ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]

			# with background correction
			# bckg = cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF[ key ]._488nm ) ]	# use same bckg for both cells
			# print(bckg)
			ax1.plot( times.timesRel-tpRel, ( cell._488nm - bckg._488nm ) / bckg._488nm, '-', color = colors[idx,jdx], lw=2 )

			# #without background correction
			# ax1.plot( times.timesRel-tpL2, cell.fluo488nm, '-', color = colors[idx,jdx], lw=2 )

	### plot the ratio
	data = pd.DataFrame( { 'tidx': timesDF.tidxRel,
							'times': timesDF.timesRel,
							'ratio': np.nan,
							'lineage': np.nan } )
	# print(data)
	colorsRatio = ['black','blue','magenta']

	for idx, trow in timesDF.iterrows():

		currentCells = extract_current_cell_fluo( cellPosDF, cellOutDF, cellFluoDF, trow.tidxRel )
		# print(currentCells)

		if len(currentCells) > 2:

			(cell1, cell2, b1, b2) = find_interesting_cells( currentCells )

			# with bckg correction
			c1signal = ( cell1._488nm - b1._488nm ) / b1._488nm
			c2signal = ( cell2._488nm - b2._488nm ) / b2._488nm

			# # without bckg correction
			# c1signal = cell1.fluo488nm
			# c2signal = cell2.fluo488nm

			# print(cell1.cname,cell1.fluo488nm,cell2.cname,cell2.fluo488nm)
			data.ix[ data.tidx == trow.tidxRel, 'ratio' ] = ( c1signal - c2signal ) / ( c1signal + c2signal )
			data.ix[ data.tidx == trow.tidxRel, 'lineage' ] = int(np.min([len(cell1.cname),len(cell2.cname)])-3)
			# print(cell1.cname,cell2.cname,int(np.min([len(cell1.cname),len(cell2.cname)])-3),colorsRatio[int(np.min([len(cell1.cname),len(cell2.cname)])-3)])

	data = data.ix[ pd.notnull( data.ratio ) ].reset_index()
	# print(data)

	for idx, d in data.iterrows():
		# print(idx,idx+1)
		ax2.plot( [ d.times-tpRel, data.times.values[np.clip(idx+1,0,len(data)-1)]-tpRel ], [ d.ratio, data.ratio.values[np.clip(idx+1,0,len(data)-1)] ], '-', color = colorsRatio[int(d.lineage)], lw=2 )

	# ### plot ecdysis
	# for idx, tidx in enumerate( lethtidx ):
	# 	if tidx >= 0:
	# 		tp = timesDF.ix[ timesDF.tidxRel == tidx, 'timesRel' ].values[0]
	# 		ax1.plot([tp-tpRel,tp-tpRel],[0,2000],'--', color = color, lw=1)
	# 		ax2.plot([tp-tpRel,tp-tpRel],[-1,1],'--', color = color, lw=1)
			# ax2.plot([tp,tp],[-1,1],'--', color = color, lw=1)
