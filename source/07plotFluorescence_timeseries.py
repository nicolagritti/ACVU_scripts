# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:31:00 2015

@author: kienle
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import cm
from generalFunctions import *


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotFluorescence( path, worm, ax1, ax2, color = 'k', channel = '488nm' ):

	timesDF = load_data_frame( path, worm + '_01times.pickle' )
	cellFluoDF = load_data_frame( path, worm + '_06cellFluo.pickle' )

	lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']]
	colors = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
				np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )

	### find ecdysis timepoint
	ecd = np.loadtxt( open( os.path.join( path, 'skin.txt'), 'rb' ) )
	# load ecdysis data
	index = np.where( ecd[:,0] == float(worm[1:]) )
	mintp = np.min( [ i for i in ecd[index, 1:6][0][0] if i >0 ] )
	lethtidx = ecd[ index, 2:6 ][0][0] - mintp
	tpL2 = timesDF.ix[ timesDF.tidxRel == lethtidx[0], 'timesRel' ].values[0]
	tpL2=0
	print(tpL2)

	### plot the timeseries
	for idx, lin in enumerate( lineages ):

		for jdx, key in enumerate( lin ):

			cell = cellFluoDF[ key ].ix[ pd.notnull( cellFluoDF[ key ].fluo488nm ) ]
			bckg = cellFluoDF[ 'b_'+key[0] ].ix[ pd.notnull( cellFluoDF[ key ].fluo488nm ) ]
			times = timesDF.ix[ pd.notnull( cellFluoDF[ key ].fluo488nm ) ]

			# with background correction
			bckg = cellFluoDF[ 'b_1' ].ix[ pd.notnull( cellFluoDF[ key ].fluo488nm ) ]
			ax1.plot( times.timesRel-tpL2, ( cell.fluo488nm - bckg.fluo488nm ) / bckg.fluo488nm, '-', color = colors[idx,jdx], lw=2 )

			# #without background correction
			# ax1.plot( times.timesRel-tpL2, cell.fluo488nm, '-', color = colors[idx,jdx], lw=2 )

	### plot the ratio
	data = pd.DataFrame( { 'tidx': timesDF.tidxRel,
							'times': timesDF.timesRel,
							'ratio': np.nan,
							'lineage': np.nan } )
	colorsRatio = ['white','black','magenta']

	for idx, trow in timesDF.iterrows():

		currentCells = extract_current_cell_fluo( cellFluoDF, trow.tidxRel )
		# print(currentCells)

		if len(currentCells) > 2:

			(cell1, cell2, b1, b2) = find_interesting_cells( currentCells )

			# with bckg correction
			c1signal = ( cell1.fluo488nm - b1.fluo488nm ) / b1.fluo488nm
			c2signal = ( cell2.fluo488nm - b2.fluo488nm ) / b2.fluo488nm

			# # without bckg correction
			# c1signal = cell1.fluo488nm
			# c2signal = cell2.fluo488nm

			# print(cell1.cname,cell1.fluo488nm,cell2.cname,cell2.fluo488nm)
			data.ix[ data.tidx == trow.tidxRel, 'ratio' ] = ( c1signal - c2signal ) / ( c1signal + c2signal )
			data.ix[ data.tidx == trow.tidxRel, 'lineage' ] = int(np.max([len(cell1.cname),len(cell2.cname)])-3)

	data = data.ix[ pd.notnull( data.ratio ) ].reset_index()
	# print(data)

	for idx, d in data.iterrows():
		# print(idx,idx+1)
		ax2.plot( [ d.times-tpL2, data.times.values[np.clip(idx+1,0,len(data)-1)]-tpL2 ], [ d.ratio, data.ratio.values[np.clip(idx+1,0,len(data)-1)] ], '-', color = colorsRatio[int(d.lineage)], lw=2 )

	### plot ecdysis
	for idx, tidx in enumerate( lethtidx ):
		if tidx >= 0:
			tp = timesDF.ix[ timesDF.tidxRel == tidx, 'timesRel' ].values[0]
			ax1.plot([tp-tpL2,tp-tpL2],[0,2000],'--', color = color, lw=1)
			ax2.plot([tp-tpL2,tp-tpL2],[-1,1],'--', color = color, lw=1)
			# ax2.plot([tp,tp],[-1,1],'--', color = color, lw=1)

if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,3.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	ax1.set_ylim((0,6))
	# ax1.set_xlim(-5,10)

	### setup figure for the ratio
	fig2 = plt.figure(figsize=(7,2.75))
	ax2 = fig2.add_subplot(111)
	fig2.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax2.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax2.get_yticklabels():
		tl.set_fontsize(18)
	ax2.set_ylim(-1,1)
	# ax2.set_xlim(-5,10)

	path = 'X:\\Simone\\160129_MCHERRY_HLH2GFP_onHB101'

	worms = ['C19']
	# worms = ['C10','C11','C14','C15','C17','C18','C19']#'C02','C03','C04','C06']

	for idx, w in enumerate( worms ):
		
		plotFluorescence( path, w, ax1, ax2, channel = '488nm' )
	ax2.plot([-5,25],[0,0],'--k')
		
	plt.show()

