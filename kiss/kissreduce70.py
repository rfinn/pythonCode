#!/usr/bin/env python
"""
re-written for KISS reduction

useage:

reduce24.py prefix

prefix is the galaxy prefix (e.g. K0225), where the catalog is called for example 
CL1216mosaic_extract.tbl

reduce24.py K0225

"""

import sys, glob, os
#import numarray as n
import numpy as N
import pylab
#from math import *
from pyraf import iraf
#import pyfits
import random
import mystuff
#bcdpath=r13810688/ch1/bcd

SqDegS=0.

#iraf.images()
#iraf.images.imutil()
#iraf.artdata()
#import spitzergetnoise

def findnearest(x1,y1,x2,y2,delta):#x2 is array
	dmin = 100000000000000000000000000000000000.
	matchflag=1
	nmatch=0
	for i in range(len(x2)):
		d = N.sqrt((x1-x2[i])**2 + (y1-y2[i])**2)
		if d < delta:
			nmatch=nmatch+1
		if d < dmin:
			dmin = d
			imatch = i

	
	if dmin > delta:
		imatch = 0
		matchflag = 0
	return imatch, matchflag,nmatch


def initialize():
    os.chdir(bcdpath)
    os.system('mkdir cdf')
    os.system('cp /Users/rfinn/research/KISS/data/cdf/mosaic_70um.nl cdf/.')
    os.system('ln -s /Applications/mopex/cal cal')
    os.system('ls SPITZER*_bcd.fits  > InputImageList.txt')
    os.system('ls SPITZER*_fbcd.fits  > fInputImageList.txt')
    os.system('ls SPITZER*bunc.fits > SigmaList.txt')
    os.system('ls SPITZER*bmask.fits > DmaskList.txt')

def create24nl():
    s="sed -e 's/K0225/"+prefix+"/g' /Applications/mopex/cdf/kiss24K0225.nl > /Applications/mopex/cdf/kiss24"+prefix+".nl"
    os.system(s)

def create70nl():
    s="sed -e 's/K1759/"+prefix+"/g' /Applications/mopex/cdf/kissmips70K1759.nl > /Applications/mopex/cdf/kissmips70"+prefix+".nl"
    os.system(s)
#    print "Date of Observation = "
#    os.system("gethead -f DATE_OBS ../pbcd/SPIT*maic.fits")
#    print "Find pmask file with the first date AFTER the date of observation"
#    print " "
#    print "Loading pbcd mosaic in ds9 for your viewing pleasure"
#    os.system("ds9 ../pbcd/SPIT*maic.fits -zscale -regions /Users/rfinn/research/KISS/RegionsFiles/kiss.reg&")



def extractphot():
    combinepath=bcdpath+'pbcd/Combine/'

    if (SqDegS > 0.1):
	    combinepath=bcdpath
    os.chdir(combinepath)
    os.system('mkdir cal')
    os.system('mkdir cdf')
    os.system('cp /Users/rfinn/research/clusters/spitzer/MIPS24_PRF_HD_center.fits cal/.')
    os.system('cp /Applications/mopex/cal/mips24_prf_mosaic_2.45_4x.fits cal/.')
    #os.system('cp /Users/rfinn/clusters/spitzer/PRF_estimate_cl1216.fits prf_estimate_8_11.nl cal/MIPS24_PRF_HD_center.fits')
    os.system('cp /Users/rfinn/research/clusters/spitzer/*.nl cdf/.')
    #os.system('cp /Users/rfinn/clusters/spitzer/apex_1frame_MIPS24_step2.nl cdf/.')
    #os.system('cp /Users/rfinn/clusters/spitzer/prf_estimate_8_11.nl cdf/.')

    os.system('apex_1frame.pl -n apex_1frame_MIPS24_step1.nl')
    input=open('output_apex_step1/mosaic_extract.tbl','r')
    output=open('output_apex_step1/mosaic_extract_clean.tbl','w')
    i=0
    for line in input:
	if line.find('\\') > -1:
	    output.write(line)
	    continue
	if line.find('|') > -1:
	    output.write(line)
	    continue
	f=line.split()#check to make sure source is far from edge
	if (float(f[1]) < 15.) :
	    continue
	if (float(f[2]) < 15.) :
	    continue
	if (float(f[1]) > 128.) :
	    continue
	if (float(f[2]) > 128.) :
	    continue
	if (i > 3):
	    print "Got 3 sources in mosaic_extract.tbl\n"
	    break
	output.write(line)
	i=i+1
    input.close()
    output.close()
    os.system('prf_estimate.pl -n prf_estimate_8_11.nl')
    os.system('cp prf_estimate/PRF.fits cal/PRF_estimate.fits')
    os.system('apex_1frame.pl -n apex_1frame_MIPS24_step2.nl')
    #os.system('apex_qa.pl -n apex_qa_myredux.nl')#produce point source subtracted image to gauge how well detection went
    #image='apex_qa/Mosaic/residual_mosaic.fits'
    #iraf.images.tv.display(image,3,contrast=0.01)#CL1040


if SqDegS < 0.1:#run for KISS galaxies
	prefix=sys.argv[1]
	bcdpath='/Users/rfinn/research/KISS/data/'+prefix+'/mips70/ch2/bcd/'
	catalogpath='/Users/rfinn/research/KISS/MIPSPhot/' #catalogpath is the directory where the final catalog will be copied to
#spitzersourcpath=bcdpath+'pbcd/Combine/output_apex_step2/'#where final catalog is kept & where spitzersource shoud be run
	spitzersourcpath=bcdpath+'pbcd/Combine/apex_1frame_step2/'#where final catalog is kept & where spitzersource shoud be run
else:
	prefix=sys.argv[1]
	bcdpath='/Users/rfinn/400sqd/completeness/'+prefix+'/'
	catalogpath=bcdpath
	spitzersourcepath=bcdpath

print "Running reduce24.py for ",prefix

#os.chdir(bcdpath)
#run from bcd directory

print 'bcdpath  = ',bcdpath
initialize()
create70nl()
#os.system('flatfield.pl -n flatfield_24_kiss.nl')
#os.system('mosaic.pl -n mosaic_24_kiss.nl')
#create24nl()#creates 24um apex_1frame nl for each galaxy


# didn't use these for KISS - ran apex through gui
#extractphot()#ran this on Ken's images
#extractphotnoring()
#writephotcat()#had to run this on Ken's images
#completeness()#then running this to estimate completeness on Ken's images

#os.chdir(spitzersourcepath)
#os.system('spitzersource.py 1')
os.system('echo All Finished!')


