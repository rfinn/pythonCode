#! /usr/bin/env python

from pylab import *
from pyraf import iraf
import glob


files=glob.glob('*.fits')

for file in files:
    iraf.imgets(image=file,param='FILENAME')#get RA of image
    originalName=iraf.imgets.value
    #print file, originalName
    iraf.imrename(file,originalName)


print "All Done!"
