# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:43:34 2016

@author: gritti
"""

import nd2reader as nd2
import matplotlib.pyplot as plt
import numpy as np
from skimage import filters

fname = 'X:\\Simone\\160428_HLH2+mCherry_confocal\\nd025.nd2'

image = nd2.Nd2(fname)

gfp = []
for img in image.select(channels='EGFP'):
    gfp.append(img)
gfp = np.array(gfp)

mcherry = []
for img in image.select(channels='mCherry'):
    mcherry.append(img)
mcherry = np.array(mcherry)

thr = []
for i in mcherry:
    thr.append( filters.threshold_otsu(i) )
print(thr)

idx = 17
plt.imshow(mcherry[idx]>thr[idx])
plt.show()

  
