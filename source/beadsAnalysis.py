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
from skimage import filters


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def plotFluorescence( path, worm, ax1, ax2 ):
	pos = np.loadtxt(open(os.path.join(path,worm+'_Pos.txt')))
	# pos = pos[28:30]
	print(pos)

	flist = glob.glob(os.path.join(path,worm,'*.tif'))
	flist.sort()
	# print(flist)

	imgs = [ load_stack(fName) for fName in flist ]
	# imgs = load_stack(flist[0])

	size = 25

	mean = []
	std = []
	for idx, p in enumerate( pos ):
		print(idx,p[2],p[1])

		movie = []
		movieMask = []
		data = []
		for img in imgs:
		
			imgCrop = img[p[2]-size:p[2]+size,p[1]-size:p[1]+size]
			movie.append(imgCrop)
			
			# plt.figure()
			# plt.imshow(imgCrop)
			# plt.show()

			thr = filters.threshold_otsu(imgCrop)
			if p[0] == 0:
				thr = 0.
			global_otsu = imgCrop >= thr
			movieMask.append( global_otsu )

			# plt.figure()
			# plt.imshow(global_otsu,interpolation='nearest')
			# # plt.show()
			# plt.figure()
			# plt.imshow(imgCrop*global_otsu,interpolation='nearest')
			# plt.show()

			data.append( np.sum(imgCrop*global_otsu) / np.sum(global_otsu) )

		if p[0] == 1:
			ax1.plot(data, '-b', lw=2)
		if p[0] == 0:
			ax1.plot(data, '-r', lw=2)

		mean.append(np.mean(data))
		std.append(np.std(data))#(np.max(data)-np.min(data)))

		# imsave( os.path.join( path, worm + '_pos'+str(idx)+'.tif' ), np.array(movie) )
		# imsave( os.path.join( path, worm + '_pos'+str(idx)+'_mask.tif' ), np.array(movieMask).astype(np.uint16) )

	img = ax2.scatter( pos[:,1] , pos[:,2],c=np.array(std)/np.array(mean),cmap=plt.cm.seismic,s=100 )
	names = np.arange(len(pos[:,1]))
	names = [ str(i) for i in names ]
	print(names)
	for idx in np.arange(len(pos[:,1])):
		ax2.text(pos[idx,1] , pos[idx,2],names[idx])

	plt.colorbar(img)
	

if __name__ == '__main__':

	### setup figure for the timeseries
	fig1 = plt.figure(figsize=(5.8,3.8))
	ax1 = fig1.add_subplot(111)
	fig1.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax1.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax1.get_yticklabels():
		tl.set_fontsize(18)
	ax1.set_ylim((0,2000))

	### setup figure for the timeseries
	fig2 = plt.figure(figsize=(5.8,3.8))
	ax2 = fig2.add_subplot(111)
	fig2.subplots_adjust(left=0.15, right=.95, top=.95, bottom=0.15)
	for tl in ax2.get_xticklabels():
		tl.set_fontsize(18)
	for tl in ax2.get_yticklabels():
		tl.set_fontsize(18)
	ax2.set_ylim((0,2048))
	ax2.set_xlim((0,2048))

	path = 'X:\\Nicola\\160609_noiseTest_beads'
	worm = 'C01'

	path = 'X:\\Nicola\\160616_testCamera'
	worm = 'beads_561'

	plotFluorescence( path, worm, ax1, ax2 )
		
	ax2.invert_yaxis()
	plt.show()

