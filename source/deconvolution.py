import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fftpack
from skimage import restoration
from numpy import linalg

def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(psf_fft*star_fft)))

def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    psf_ifft = linalg.inv(psf_fft)
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(psf_fft*star_fft)))

def makeGaussian(size=100, fwhm = 10, center=None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    g = np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    return g/np.sum(g)

def makeDisk(size = 100, A = 10, radius = 10, center = None):
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    
    return A * ( ( (x-x0)**2 + (y-y0)**2 ) < radius**2 )

def RL_deconv(img, psf, iterations):
	fn = img
	otf = fftpack.fftshift(fftpack.fft2(psf))
	otfhat = fftpack.fftshift(fftpack.fft2(psf[::-1,::-1]))
	# print(np.min(otf),np.max(otf))
	# plt.figure()
	# plt.imshow(np.real(otf))
	for i in np.arange(iterations):
		ffn = fftpack.fftshift(fftpack.fft2(fn))
		Hfn = otf * ffn
		iHfn = fftpack.ifft2(fftpack.ifftshift(Hfn))
		ratio = img // iHfn
		iratio = fftpack.fftshift(fftpack.fft2(ratio))
		res = otfhat * iratio
		ires = fftpack.ifft2(fftpack.ifftshift(res))
		fn = ires * fn
	return np.abs(fn)

psf = makeGaussian(1000, fwhm = 50)
print('psf',np.min(psf),np.max(psf))
plt.figure()
plt.title('psf')
plt.imshow(psf)

obj = 10+makeDisk(1000,A=1000,radius = 100,center = [400,600])
print('obj',np.min(obj),np.max(obj))
plt.figure()
plt.title('obj')
plt.imshow(obj)

blur = signal.fftconvolve(obj,psf,'same')
blur += .00001 * blur.std() * np.random.standard_normal(blur.shape)
# blur -= np.min(blur)
#blur = (blur/np.max(blur)) * np.max(obj)# + 2 * np.random.randn(1000, 1000)
print('blur',np.min(blur),np.max(blur))
plt.figure()
plt.title('blur')
plt.imshow(blur)

deconv = RL_deconv(blur, psf, 50)
# deconv = restoration.unsupervised_wiener(blur,psf)
# #blur = (blur/np.max(blur)) * np.max(obj)# + 2 * np.random.randn(1000, 1000)
# print('deconv',np.min(deconv),np.max(deconv))
plt.figure()
plt.title('deconv')
plt.imshow(deconv)

plt.show()
