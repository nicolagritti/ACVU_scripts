# import glob
# from tifffile import *
# from generalFunctions import *
# import numpy as np
# import PIL
# from PIL import Image, ImageDraw, ImageFont
# import os.path
# import matplotlib.pyplot as plt
# import pickle
# import scipy.interpolate as ip
# from scipy.ndimage import interpolation
# import matplotlib as mpl
# mpl.rcParams['pdf.fonttype'] = 42

# # path = 'X:\\Simone\\16_01_25_MCHERRY_LAG2YFP\\'
# path = 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry'

# worms = ['C02','C05','C07','C08','C11','C12','C15']
# channel = '488nm'

# data = []

# for w in worms:

# 	### load parameters and times dataframes
# 	paramsDF = load_data_frame( path, w + '_01params.pickle' )
# 	timesDF = load_data_frame( path, w + '_01times.pickle' )
# 	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
# 	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
# 	cellOutDF = load_data_frame( path, w + '_05cellOut.pickle' )
# 	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
# 	# apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )

# 	# find tDiv
# 	tidxdiv = [ np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel ), 
# 				np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].tidxRel ) ]

# 	# extract dt
# 	dt = timesDF.ix[ timesDF.tidxRel == tidxdiv[0], 'timesRel' ].values[0] - timesDF.ix[ timesDF.tidxRel == tidxdiv[1], 'timesRel' ].values[0]

# 	# find last tp
# 	tidxfinal = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel )

# 	print(tidxdiv, tidxfinal, dt)

# 	fluo = []
# 	for cname in ['1.ppp','4.aaa']:
# 		f = cellFluoDF[ cname ].ix[ cellFluoDF[cname].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		b = cellFluoDF[ 'b_'+cname[0] ].ix[ cellFluoDF[ 'b_'+cname[0] ].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		fluo.append((f-b)/b)
# 	initial = ( fluo[0] - fluo[1] ) / ( fluo[0] + fluo[1] )

# 	fluo = []
# 	for cname in ['1.ppp','4.aaa']:
# 		f = cellFluoDF[ cname ].ix[ cellFluoDF[cname].tidx == tidxfinal, '_488nm' ].values[0]
# 		b = cellFluoDF[ 'b_'+cname[0] ].ix[ cellFluoDF[ 'b_'+cname[0] ].tidx == tidxfinal, '_488nm' ].values[0]
# 		fluo.append((f-b)/b)
# 	final = ( fluo[0] - fluo[1] ) / ( fluo[0] + fluo[1] )

# 	data.append([dt,initial,final])

# data = np.array(data)
# print(data)

# fig1 = plt.figure(figsize=(5.8,5.8))
# ax1 = fig1.add_subplot(111)
# fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
# for tl in ax1.get_xticklabels():
# 	tl.set_fontsize(18)
# for tl in ax1.get_yticklabels():
# 	tl.set_fontsize(18)
# ax1.set_ylim((-1,1))

# # ax1.set_xlim((9,30))
# ax1.set_xlim(-1.5,1.5)

# ax1.plot(data[:,0],data[:,1],'or',ms=10)
# ax1.plot(data[:,0],data[:,2],'ob',ms=10)
# for d in data:
# 	ax1.plot([d[0],d[0]],[d[1],d[2]],'-k')

# ax1.plot([-2,2],[0,0],'--k')
# ax1.plot([0,0],[-1,1],'--k')

# plt.show()












# import glob
# from tifffile import *
# from generalFunctions import *
# import numpy as np
# import PIL
# from PIL import Image, ImageDraw, ImageFont
# import os.path
# import matplotlib.pyplot as plt
# import pickle
# import scipy.interpolate as ip
# from scipy.ndimage import interpolation
# import matplotlib as mpl
# mpl.rcParams['pdf.fonttype'] = 42

# # path = 'X:\\Simone\\16_01_25_MCHERRY_LAG2YFP\\'
# path = 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry'

# worms = ['C02','C05','C07','C08','C11','C12','C15']
# channel = '488nm'

# data = []

# fig1 = plt.figure(figsize=(5.8,5.8))
# ax1 = fig1.add_subplot(111)
# fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
# for tl in ax1.get_xticklabels():
# 	tl.set_fontsize(18)
# for tl in ax1.get_yticklabels():
# 	tl.set_fontsize(18)
# ax1.set_ylim((-6,6))

# # ax1.set_xlim((9,30))
# ax1.set_xlim(-.6,1.1)

# for w in worms:

# 	### load parameters and times dataframes
# 	paramsDF = load_data_frame( path, w + '_01params.pickle' )
# 	timesDF = load_data_frame( path, w + '_01times.pickle' )
# 	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
# 	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
# 	cellOutDF = load_data_frame( path, w + '_05cellOut.pickle' )
# 	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
# 	# apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )

# 	# find tDiv
# 	tidxdiv = [ np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel ), 
# 				np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].tidxRel ) ]

# 	# find last tp
# 	tidxmax = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel )

# 	dt = timesDF.ix[ timesDF.tidxRel == tidxdiv[0], 'timesRel' ].values[0] - timesDF.ix[ timesDF.tidxRel == tidxdiv[1], 'timesRel' ].values[0]

# 	data = [dt]
# 	colors = ['black','blue','red','green','magenta']
# 	ms=10

# 	if tidxdiv[0] == tidxdiv[1]:
# 		tidxs = [np.min(tidxdiv)-1,np.min(tidxdiv),tidxmax]

# 		f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.pp'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aa'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == tidxmax, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == tidxmax, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		ax1.plot(data[0],data[1],'o',ms=ms,color=colors[0],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[2],'o',ms=ms,color=colors[3],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[3],'o',ms=ms,color=colors[4],alpha = .5,mew=0)

# 		ax1.plot([data[0],data[0]],[data[1],data[2]],'-k')
# 		ax1.plot([data[0],data[0]],[data[2],data[3]],'-k')

# 	if tidxdiv[0] > tidxdiv[1]:
# 		tidxs = [np.min(tidxdiv)-1,np.min(tidxdiv),np.max(tidxdiv)-1,np.max(tidxdiv),tidxmax]

# 		f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.pp'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aa'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.pp'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.pp'].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == tidxmax, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == tidxmax, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		ax1.plot(data[0],data[1],'o',ms=ms,color=colors[0],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[2],'o',ms=ms,color=colors[1],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[3],'o',ms=ms,color=colors[2],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[4],'o',ms=ms,color=colors[3],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[5],'o',ms=ms,color=colors[4],alpha = .5,mew=0)

# 		ax1.plot([data[0],data[0]],[data[1],data[2]],'-k')
# 		ax1.plot([data[0],data[0]],[data[2],data[3]],'-k')
# 		ax1.plot([data[0],data[0]],[data[3],data[4]],'-k')
# 		ax1.plot([data[0],data[0]],[data[4],data[5]],'-k')

# 	if tidxdiv[0] < tidxdiv[1]:
# 		tidxs = [np.min(tidxdiv)-1,np.min(tidxdiv),np.max(tidxdiv)-1,np.max(tidxdiv),tidxmax]

# 		f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.pp'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aa'].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv )-1, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aa'].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.min( tidxdiv ), '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aa'].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.max( tidxdiv )-1, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == np.max( tidxdiv ), '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == tidxmax, '_488nm' ].values[0]
# 		b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r1 = (f1-b1)/b1
# 		f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == tidxmax, '_488nm' ].values[0]
# 		b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == tidxmax, '_488nm' ].values[0]
# 		r4 = (f4-b4)/b4
# 		data.append((r1-r4))

# 		ax1.plot(data[0],data[1],'o',ms=ms,color=colors[0],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[2],'o',ms=ms,color=colors[1],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[3],'o',ms=ms,color=colors[2],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[4],'o',ms=ms,color=colors[3],alpha = .5,mew=0)
# 		ax1.plot(data[0],data[5],'o',ms=ms,color=colors[4],alpha = .5,mew=0)

# 		ax1.plot([data[0],data[0]],[data[1],data[2]],'-k')
# 		ax1.plot([data[0],data[0]],[data[2],data[3]],'-k')
# 		ax1.plot([data[0],data[0]],[data[3],data[4]],'-k')
# 		ax1.plot([data[0],data[0]],[data[4],data[5]],'-k')




# ax1.plot([-6,6],[0,0],'--k')
# ax1.plot([0,0],[-6,6],'--k')

# plt.show()













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
worms = ['C02','C06','C07','C09','C14','C16','C17','C19','C20','C21','C22','C24','C25']
firsttidx = [65, 90,   84,   85,   56,   78,   87,   92,   25,   44,   70,   42,   32]
firsttidx = [60, 85,   80,   79,   51,   74,   79,   88,   21,   42,   61,   37,   24]
lasttidx  = [105,130, 139,  126,   95,  120,  135,  130,   70,   88,  110,   85,   77]

# ### PART TO BE CHANGED
# path = 'X:\\Simone\\160407_HLH2_GFP_hist_mCherry'
# worms = ['C02','C05','C07','C08','C11','C12','C15','C19','C20','C21','C27']
# firsttidx = [44, 58, 58, 58, 50, 51, 52, 53, 55, 58, 51]
# ###

data = []

fig1 = plt.figure(figsize=(6.8,5.8))
ax1 = fig1.add_subplot(111)
fig1.subplots_adjust(left=0.1, right=.9, top=.9, bottom=0.1)
for tl in ax1.get_xticklabels():
	tl.set_fontsize(18)
for tl in ax1.get_yticklabels():
	tl.set_fontsize(18)
# ax1.set_ylim((-6,6))

# ax1.set_xlim((9,30))
# ax1.set_xlim(-.6,1.1)

data = []
dts = []

for widx, w in enumerate(worms):

	### load parameters and times dataframes
	paramsDF = load_data_frame( path, w + '_01params.pickle' )
	timesDF = load_data_frame( path, w + '_01times.pickle' )
	gpDF = load_data_frame( path, w + '_02gonadPos.pickle' )
	cellPosDF = load_data_frame( path, w + '_04cellPos.pickle' )
	cellOutDF = load_data_frame( path, w + '_05cellOut.pickle' )
	cellFluoDF = load_data_frame( path, w + '_06cellFluo.pickle' )
	# apdvPosDF = load_data_frame( path, w + '_08apdvPos.pickle' )

	### find division time delay

	# find tDiv
	tidxdiv = [ np.min( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel ), 
				np.min( timesDF.ix[ pd.notnull( cellPosDF[ '4.aaa' ].X ) ].tidxRel ) ]
	print(tidxdiv)

	dt = timesDF.ix[ timesDF.tidxRel == tidxdiv[0], 'timesRel' ].values[0] - timesDF.ix[ timesDF.tidxRel == tidxdiv[1], 'timesRel' ].values[0]
	dts.append(dt)

	### find fluorescence intensity difference

	# automatic detection of first tp with Z1.ppp and Z4.aaa
	# tdiv = np.max(tidxdiv)
	# manual detection
	tdiv = firsttidx[widx]
	# autmoatic detection of mother cells
	# tdiv = np.min(tidxdiv)-4

	# automatic detection of last tp with Z1.ppp and Z4.aaa
	tidxmax = np.max( timesDF.ix[ pd.notnull( cellPosDF[ '1.ppp' ].X ) ].tidxRel )
	# manual detection
	# tidxmax = lasttidx[widx]

	f1 = cellFluoDF[ '1.pp' ].ix[ cellFluoDF['1.ppp'].tidx == tdiv, '_488nm' ].values[0]
	b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == tdiv, '_488nm' ].values[0]
	r1 = (f1-b1)/b1
	f4 = cellFluoDF[ '4.aa' ].ix[ cellFluoDF['4.aaa'].tidx == tdiv, '_488nm' ].values[0]
	b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == tdiv, '_488nm' ].values[0]
	r4 = (f4-b4)/b4
	data1 = ( r1 - r4 ) / ( r1 + r4 )

	f1 = cellFluoDF[ '1.ppp' ].ix[ cellFluoDF['1.ppp'].tidx == tidxmax, '_488nm' ].values[0]
	b1 = cellFluoDF[ 'b_1' ].ix[ cellFluoDF[ 'b_1' ].tidx == tidxmax, '_488nm' ].values[0]
	r1 = (f1-b1)/b1
	f4 = cellFluoDF[ '4.aaa' ].ix[ cellFluoDF['4.aaa'].tidx == tidxmax, '_488nm' ].values[0]
	b4 = cellFluoDF[ 'b_4' ].ix[ cellFluoDF[ 'b_4' ].tidx == tidxmax, '_488nm' ].values[0]
	r4 = (f4-b4)/b4
	data2 = ( r1 - r4 ) / ( r1 + r4 )

	ax1.text(dt,data2,w)

	### append data to full list

	data.append([data1,data2])

data = np.array(data)

### define my colormap
cdict1 = {'red':   ((0.0, 0.6, 0.6),
					(0.25,0.0,0.0),
                    (0.5, 0.0, 0.0),
					(0.75,1.0,1.0),
                    (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.6, 0.6),
				   (0.25,0.0,0.0),
		           (0.5, 0.0, 0.0),
				   (0.75,0.0,0.0),
                   (1.0, 0.6, 0.6)),

         'blue':  ((0.0, 1.0, 1.0),
				   (0.25,1.0,1.0),
                   (0.5, 0.0, 0.0),
				   (0.75,0.0,0.0),
                   (1.0, 0.6, 0.6))
         }

mycmap = LinearSegmentedColormap('BlueRed1', cdict1)

print(data)
# sc = ax1.scatter(dts,data[:,1],c = data[:,0],cmap = mycmap,s=60,linewidth=0,vmin=-1,vmax=1)
sc = ax1.scatter(dts,data[:,1],c = data[:,0],cmap = plt.cm.seismic,s=60,linewidth=1,vmin=-.5,vmax=.5)
# sc = ax1.scatter(dts,data[:,1],s=60,linewidth=1,vmin=-.5,vmax=.5,color='black')

ax1.plot([-1.5,1.5],[0,0],'--k')
ax1.plot([0,0],[-5,5],'--k')

ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1,1)

ax1.text(2.5,-3,'initial: 4>1',rotation=90,verticalalignment = 'center',fontsize = 18)
ax1.text(2.5,+3,'initial: 1>4',rotation=90,verticalalignment = 'center',fontsize = 18)

ax1.text(-1.9,+2.5,'final: 1>4',rotation=90,verticalalignment = 'center',fontsize = 18)
ax1.text(-1.9,-2.5,'final: 4>1',rotation=90,verticalalignment = 'center',fontsize = 18)

ax1.text(-.75,-6,'1->4',horizontalalignment = 'center',fontsize = 18)
ax1.text(+.75,-6,'4->1',horizontalalignment = 'center',fontsize = 18)

plt.colorbar(sc)
plt.show()
