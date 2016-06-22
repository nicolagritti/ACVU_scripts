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
	size = 128
	pos = ( 1+np.arange(2048/(2*size)) ) *2*size - size
	# pos=[1,2,3]
	pos = np.array( [ [ i, j ] for i in pos for j in pos ] )
	# print(pos)

	print(os.path.join(path,worm,'*.tif'))
	flist = glob.glob(os.path.join(path,worm,'*.tif'))
	flist.sort()
	print(flist)

	# imgs = np.array( [ load_stack(fName) for fName in flist ] )
	imgs = load_stack(flist[0])

	mean = []
	std = []
	for idx, p in enumerate( pos ):

		print(p)

		# movie = []
		data = []
		for img in imgs:
		
			imgCrop = img[p[1]-size:p[1]+size,p[0]-size:p[0]+size]
			

			# plt.figure()
			# plt.imshow(global_otsu,interpolation='nearest')
			# # plt.show()
			# plt.figure()
			# plt.imshow(imgCrop*global_otsu,interpolation='nearest')
			# plt.show()

			data.append( np.mean(imgCrop) )

		ax1.plot(data, '-b', lw=.5)

		# mean.append(np.mean(data))
		# std.append(np.std(data))#(np.max(data)-np.min(data)))

		# imsave( os.path.join( path, worm + '_pos'+str(idx)+'.tif' ), np.array(movie) )
	# img = ax2.scatter( pos[:,0], pos[:,1], c=np.array(std)/np.array(mean),cmap=plt.cm.seismic,s=100 )

	std = np.std(imgs,0)
	mean = np.mean(imgs,0)
	print(std.shape)
	img = ax2.imshow( std/mean , interpolation = 'nearest', cmap= 'gray')#, vmin=0.05, vmax = 0.1)
	
	imsave(os.path.join(path,worm+'_std.tif'), std.astype(np.uint16))
	imsave(os.path.join(path,worm+'_mean.tif'), mean.astype(np.uint16))
	imsave(os.path.join(path,worm+'_stdovermean.tif'), ((std/mean)*(2**16-1)).astype(np.uint16))

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
	# ax1.set_ylim((0,2000))

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

	path = 'X:\\Nicola\\160614_noiseTest'
	worm = 'beads_561_noExcFilter'

	path = 'X:\\Nicola\\160620_noiseTest'
	worm = 'Home'

	plotFluorescence( path, worm, ax1, ax2 )
		
	ax2.invert_yaxis()
	plt.show()

