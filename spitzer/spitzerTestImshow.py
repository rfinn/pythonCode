#!/usr/bin/env python

from pylab import *
import pyfits

fits=pyfits.open('/Users/rfinn/clusters/spitzer/GroupLirgs/test24.fits')
im=fits[0].data.copy()
#imshow(im,interpolation='nearest')
imshow(im,interpolation='nearest',vmin=0,vmax=0.09)#,cmap='gray')
savefig('test.eps')
