#!/usr/bin/env python
"""
useage:

reduce24.py bcdpathrvalue prefix

where bcdpathrvalue is the absolute path of the bcd data

prefix is the cluster prefix, where the catalog is called for example 
CL1216mosaic_extract.tbl

reduce24.py r13810944 cl1354

"""

import sys, glob, os
#import numarray as n
import Numeric as N
#import pylab
#from math import *
from pyraf import iraf
#import pyfits
#bcdpath=r13810688/ch1/bcd
iraf.images()
iraf.images.imutil()

def initialize():
    os.system('mkdir cdf')
    os.system('cp /Users/rfinn/clusters/spitzer/flatfield_24_ediscs.nl cdf/flatfield_24_ediscs.nl')
#os.system('cp /Users/rfinn/clusters/spitzer/flatfield_24_ediscs.2.nl cdf/.')
    os.system('cp /Users/rfinn/clusters/spitzer/mosaic_24_ediscs.nl cdf/mosaic_24_ediscs.nl')
    os.system('ln -s ~/mopex_030106/cal cal')
    os.system('mkdir firstframes')
    os.system('mv SPITZER*_0000_1_*.fits firstframes/')
    os.system('ls SPITZER*bcd.fits  > InputImageList.txt')
    os.system('ls SPITZER*bunc.fits > SigmaList.txt')
    os.system('ls SPITZER*bbmsk.fits > DmaskList.txt')


def irafflatten():
    os.system('cp /Users/rfinn/clusters/spitzer/flatsexfiles/* .')
    infile=open('InputImageList.txt','r')
    outfile=open('FlatImageList.txt','w')
    sky=[]
    for line in infile:
	im=line[0:(len(line)-1)]
	mask='mask'+im
	skyim='s'+im
	outline='f'+line
	iraf.imgets(im,'DRIBKGND')
	t=iraf.imgets.value
	sky.append(float(t))
	outfile.write(outline)
    #get object positions using sextractor
	iraf.imarith(im,'-',t,skyim)#subtract sky before running sextractor (otherwise it doesn't detect any objects - don't know why...)
	s='sex '+skyim
	os.system(s)
	x=[]
	y=[]
	catfile=open('test.cat','r')
	for line in catfile:
	    if line.find('#') > -1:
		continue
	    t=line.split()
	    x.append(float(t[10]))
	    y.append(float(t[11]))
	catfile.close()
	x=N.array(x,'f')
	y=N.array(y,'f')
        
	try:#create image of ones same size as image
	    iraf.imarith(im,'/',im,'mask')
	except:
	    iraf.imdelete('mask')
	    iraf.imarith(im,'/',im,'mask')
	print "masking objects"
	for j in range(len(x)): #mask objects and radius around each position using imreplace, radius=10 pix
	    for k in range(11):
		y1=int(y[j]+5-k)
		if y1 < 1:
		    continue
		if y1 > 128:
		    continue
		xmin=int(x[j]-5)
		if xmin < 1:
		    xmin=1
		xmax=int(x[j]+5)
		if xmax > 128:
		    xmax=128

		s='mask['+str(xmin)+':'+str(xmax)+","+str(y1)+":"+str(y1)+"]"
		iraf.imreplace(s,0.)
	iraf.imrename('mask',mask)
	print "updating BPM field in header"
	iraf.hedit(im,fields='BPM',value=mask,add='yes',verify='no',update='yes')
    outfile.close()
    infile.close()
    avesky=N.average(N.array(sky,'f'))
    lthresh=avesky-1.
    hthresh=avesky+.6
    iraf.imcombine('@InputImageList.txt','flat',combine='average',reject='ccdclip',scale='none',zero='mode',lthreshold=lthresh,hthreshold=hthresh,lsigma=2.,hsigma=2.,rdnoise=5.,gain=5.,blank=1.,grow=12.,masktype='badvalue',maskvalue='0')
    t=iraf.imstat('flat',fields='mean',format='no',Stdout=1)
    ave=float(t[0])
    iraf.imarith('flat','/',ave,'nflat')
    iraf.imarith('@InputImageList.txt','/','nflat','@FlatImageList.txt')

def extractphot():
#    combinepath=bcdpath+'pbcd/Combine/'
    combinepath='/Users/rfinn/clusters/spitzer/coma/mips/24/r4740096/ch1/pbcd'
    os.chdir(combinepath)
    os.system('mkdir cal')
    os.system('mkdir cdf')
    os.system('cp /Users/rfinn/clusters/spitzer/MIPS24_PRF_HD_center.fits cal/.')
    os.system('cp /Users/rfinn/clusters/spitzer/apex_1frame_MIPS24_step1.nl cdf/.')
    os.system('cp /Users/rfinn/clusters/spitzer/apex_1frame_MIPS24_step2.nl cdf/.')
    os.system('cp /Users/rfinn/clusters/spitzer/prf_estimate_8_11.nl cdf/.')
    os.system('apex_1frame.pl -n apex_1frame_MIPS24_step1.nl')
    input=open('output_apex_step1/mosaic_extract.tbl','r')
    output=open('output_apex_step1/mosaic_extract_clean.tbl','w')
    i=0
    xmin=100.
    xmax=600.
    ymin=100.
    ymax=2800.
    for line in input:
	if line.find('\\') > -1:
	    output.write(line)
	    continue
	if line.find('|') > -1:
	    output.write(line)
	    continue
	f=line.split()#check to make sure source is far from edge

	if (float(f[1]) < xmin) :
	    continue
	if (float(f[2]) < ymin) :
	    continue
	if (float(f[1]) > xmax) :
	    continue
	if (float(f[2]) > ymax) :
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
    string='cp output_apex_step2/mosaic_extract.tbl '+catalogpath+prefix+'mosaic_extract.tbl'
    os.system(string)


bcdpathrvalue=sys.argv[1]
prefix=sys.argv[2]
bcdpath='/Users/rfinn/clusters/spitzer/'+prefix+'/mips/24/'+bcdpathrvalue+'/ch1/bcd/'
catalogpath='/Users/rfinn/clusters/spitzer/MasterTables/' #catalogpath is the directory where the final catalog will be copied to

#os.chdir(bcdpath)

#initialize()
#irafflatten()
##os.system('flatfield.pl -n flatfield_24_ediscs.nl')
##os.system('flatfield.pl -n flatfield_24_ediscs.2.nl')
#os.system('mosaic.pl -n mosaic_24_ediscs.nl')
extractphot()
os.system('echo All Finished!')


