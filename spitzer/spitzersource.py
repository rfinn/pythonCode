#!/usr/bin/env python
import sys, os
import Numeric as N
#import scipy
from math import *
import mystuff as my
#import ppgplot
import random
import sets
import pylab
from pyraf import iraf


class Spitzer24:
    def __init__(self):#individual galaxy properties
        print "dude - 24 micron data!"

    def readfile(self,file):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.id = N.zeros(ngal,'f')
        self.imagex = N.zeros(ngal,'f')
        self.imagey  = N.zeros(ngal,'f')
        self.ra = N.zeros(ngal,'f')
        self.dec = N.zeros(ngal,'f')
        self.f = N.zeros(ngal,'f')
        self.errf = N.zeros(ngal,'f')

        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
	    #print i, line
#	    if mode == 1:
	    (self.id[i],self.imagex[i],self.imagey[i],self.ra[i],self.dec[i],self.f[i],self.errf[i])=(float(t[0]),float(t[8]),float(t[10]),float(t[3]),float(t[5]),float(t[13]),float(t[14]))
            i=i+1
        input.close()
        outfile=open('xy24.dat','w')
        for i in range(len(self.imagex)):
            #print self.imagex[i],self.imagey[i],self.f[i],self.errf[i]
            #string="%8.2f %8.2f %8.1f %8.2f \n" % (self.imagex[i],self.imagey[i],self.f[i],self.errf[i])
	    string="%8.2f %8.2f %8.1f \n" % (self.imagex[i],self.imagey[i],self.f[i])
            outfile.write(string)
        outfile.close()


    def readfile2(self,file):
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            ngal=ngal+1
        input.close()


        self.id = N.zeros(ngal,'f')
        self.imagex = N.zeros(ngal,'f')
        self.imagey  = N.zeros(ngal,'f')
        self.ra = N.zeros(ngal,'f')
        self.dec = N.zeros(ngal,'f')
        self.f = N.zeros(ngal,'f')
        self.errf = N.zeros(ngal,'f')

        input=open(file,'r')
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
	    #print i, line
            (self.id[i],self.imagex[i],self.imagey[i])=(float(t[0]),float(t[1]),float(t[2]))
            i=i+1
        input.close()
        outfile=open('xy24.dat','w')
        for i in range(len(self.imagex)):
            #print self.imagex[i],self.imagey[i],self.f[i],self.errf[i]
            string="%8.2f %8.2f %8.1f %8.2f \n" % (self.imagex[i],self.imagey[i],self.f[i],self.errf[i])
            outfile.write(string)
        outfile.close()

def displaysub(file1,file2,image):

    g24 = Spitzer24()
    os.system('rm xy24.dat')
    g24.readfile(file1)#extracted
    try:
	iraf.images.tv.display(image,2,contrast=0.01)#CL1040
	iraf.images.tv.display(image,frame,contrast=0.01)#CL1040
	iraf.images.tv.tvmark(frame,'xy24.dat',mark='circle', radii=1,color=205)#mark positions of raw file in blue
    except:
	print "Error trying to display.  Make sure ds9 is open."
    g24.readfile2(file2)#detected
    try:
	#iraf.images.tv.display(image,2,contrast=0.01)#CL1040
	#iraf.images.tv.display(image,frame,contrast=0.01)#CL1040
	iraf.images.tv.tvmark(frame,'xy24.dat',mark='circle', radii=3,color=204)#mark positions of raw file in blue
    except:
	print "Error trying to display.  Make sure ds9 is open."


            
#mode=int(sys.argv[1])#equal 1 to plot sources when you are not doing airy ring subtraction, 2 when doing airy ring subtraction
frame=1
mode=2
if (int(mode) == 1):
    file1='mosaic_extract_final.tbl'
    image='mosaic_minus_median_extract.fits'
    file2='mosaic_extract_raw.tbl'
    displaysub(file1,file2,image)

if (int(mode) == 2):
    #file1='residual_mosaic_detect.tbl'
    #image='residual_mosaic_minus_median_detect.fits'
    #file2='residual_mosaic_detect_raw.tbl'
    file1='mosaic_extract.tbl'
    image='mosaic_minus_median_detect.fits'
    file2='mosaic_detect_raw.tbl'
    displaysub(file1,file2,image)


