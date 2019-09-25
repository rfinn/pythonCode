#!/usr/bin/env python
"""
useage !!!
 spitzergetnoise.py 1 
 1= measure sky noise using iraf,
 uses mosaic_minus_median_extract.fits to measure noise
 uses mosaic_extract.tbl to find position of sources
"""
import sys, os
import numarray as N
import scipy
from math import *
from mystuff import *
import ppgplot
import random
from pyraf import iraf
from pyraftest import *
import random
print sys.argv[0],sys.argv[1]

#npoints=sys.argv[3]
npoints=500
#print "npoints = ",npoints
nap=10

xmin=1
xmax=1


#runiraf()


def getpositions(ximage,yimage,xmax,ymax):
    dmin = 5. #minimum distance to object + sqrt(isoarea)
    xpos=[]
    ypos=[]
    d = N.zeros(len(ximage),'f')
    #print len(d)
    i=0
    while i < npoints:    
        xtemp=int(round(random.uniform(5.,125.)))#random.uniform(5.,125.)#*xmax
        ytemp=int(round(random.uniform(5.,140.)))#random.uniform(5.,140.)#*ymax
	if (xtemp > 10.) & (ytemp < (ymax-10)):
	    d = N.sqrt((ximage-xtemp)**2+(yimage-ytemp)**2)
            if (min(d) > dmin):    
		if i > 2:
		    xap=N.array(xpos,'f')
		    yap=N.array(ypos,'f')
		    d2=N.sqrt((xtemp-xap)**2 + (ytemp-yap)**2)
		    if (min(d2) > dmin):
			xpos.append(xtemp)
			ypos.append(ytemp)
			i=i+1

	    else:
		xpos.append(xtemp)
		ypos.append(ytemp)
		i=i+1
    xpos=N.array(xpos,'f')
    ypos=N.array(ypos,'f')
    return xpos,ypos


def calcavesky():
    input=open("noise.dat",'r')
    aperture=[]
    counts=[]
    area=[]
    j=0
    for line in input:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        if line.find('mosaic_minus') > -1: #skip lines with '#' in them
            j=0
            continue
        j=j+1
        if (j > 3):
	    #print j, line
            t = line.split()
            aperture.append(float(t[0]))
            counts.append(float(t[1]))
            area.append(float(t[2]))
    input.close()
    aperture=N.array(aperture,'f')
    counts=N.array(counts,'f')
    area=N.array(area,'f')
    ap=N.zeros(npoints,'f')

    aparea=N.zeros(nap,'f')
    aveap=N.zeros(nap,'f')
    aveaperr=N.zeros(nap,'f')
    avearea=N.zeros(nap,'f')
    aveareaerr=N.zeros(nap,'f')
    #for i in range(len(ap)):
    for i in range(nap):
        #print i, len(ap),aperture[i],aperture[i+1]
        if ( i < (nap-1)):
            ap=N.compress((aperture >= aperture[i]) & (aperture < aperture[i+1]),counts)
            aparea=N.compress((aperture >= aperture[i]) & (aperture < aperture[i+1]),area)
        else:
            ap=N.compress((aperture >= aperture[i]) & (aperture < 20.),counts)
            aparea=N.compress((aperture >= aperture[i]) & (aperture < 20.),area)
        
        #print ap
        #aparea=N.compress((aperture >= aperture[i]) & (aperture < aperture[i+1]),area)
        aveap[i]=N.average(ap)
        aveaperr[i]=scipy.stats.std(ap)
        avearea[i]=N.average(aparea)
        aveareaerr[i]=scipy.stats.std(aparea)
        print "ave sky = %8.4f +/- %8.4f" % (N.average(ap),scipy.stats.std(ap))
        print "ave area = %8.4f +/- %8.4f" % (N.average(aparea),scipy.stats.std(aparea))
    return aveap,aveaperr,avearea,aveareaerr

def solveforab(aveaperr,avearea):#(a,b)=solveforab(aveap,avearea)
    amin=.0
    amax=1
    bmin=.0
    bmax=1
    step=.005
    a=N.arange(amin,amax,step,'f')
    b=N.arange(bmin,bmax,step,'f')
    y=N.zeros(len(aveaperr),'f')
    diff=N.zeros(len(aveaperr),'f')
    mini=1000000000.
    for ai in a:
        for bi in b:
            y=N.zeros(len(avearea),'f')
            y=N.sqrt(avearea)*ai*(1.+bi*N.sqrt(avearea))
            #diff = N.sqrt((y - aveap)**2)
            diff = (y - aveaperr)
            sumdiff = N.sum(abs(diff))
            #print "%5.3f %5.3f %8.2f %8.2f" %(ai,bi,sumdiff,mini)
            if sumdiff < mini:
                afinal=ai
                bfinal=bi
                mini=sumdiff
    return afinal,bfinal

def readfile():
    input=open('mosaic_extract.tbl','r')
    number=[]
    fluxiso = []
    fluxerriso = []
    ximage = []
    yimage = []
    for line in input:
	if line.find('NAXIS1') > -1: #skip lines with '#' in them
	    f=line.split()
	    xmax=float(f[3])
	    continue
	if line.find('NAXIS2') > -1: #skip lines with '#' in them
	    f=line.split()
	    ymax=float(f[3])
	    continue
	if line.find('#') > -1: #skip lines with '#' in them
	    continue
	if line.find('\\') > -1: #skip lines with '#' in them
	    continue
	if line.find('|') > -1: #skip lines with '#' in them
	    continue
	t = line.split()
	number.append(float(t[0]))
	fluxiso.append(float(t[5]))
	fluxerriso.append(float(t[6]))
	ximage.append(float(t[1]))
	yimage.append(float(t[2]))
    number=N.array(number,'f')
    fluxiso=N.array(fluxiso,'f')
    fluxerriso=N.array(fluxerriso,'f')
    ximage=N.array(ximage,'f')
    yimage=N.array(yimage,'f')
    input.close()
    return number,ximage,yimage,fluxiso,fluxerriso,xmax,ymax

def doiraf(ximage,yimage,xmax,ymax):
    (xpos,ypos) = getpositions(ximage,yimage,xmax,ymax)
    coords=open("noisecoords.dat",'w')
    for i in range(len(xpos)):
        coords.write("%8.1f %8.1f \n" % (xpos[i],ypos[i]))
    coords.close()
    sky = open("sky",'w')
    for i in range(npoints):
        sky.write("0.0 \n")
    sky.close()
    aps = open("apertures",'w')
    aps.write("1,1.5,2,2.6,3,3.5,4,4.5,5,5.5")
    aps.close()

    #runiraf()
    iraf.digiphot()
    iraf.daophot()
    ##iraf.apphot()
    image = 'mosaic_minus_median_extract.fits'
    #print image
    os.system("rm noise.dat")
    iraf.digiphot.daophot.phot(image,coords="noisecoords.dat",output="noise.dat",calgorithm='none',skyfile="sky",salgori="file",aperture="apertures",interactive="no",verify='no',verbose='no')
    ##iraf.digiphot.apphot.phot(image,coords="noisecoords.dat",output="noise.dat",calgorithm='none',skyfile="sky",salgori="file",aperture="apertures",interactive="no",verbose='yes')



#contsub=[]
#contsuberr=[]
#contsubisoarea=[]
#realnoise = open("contsub-noise-area.dat",'r')
#for line in realnoise:
#    t = line.split()
#    contsub.append(float(t[0]))
#    contsuberr.append(float(t[1]))
#    contsubisoarea.append(float(t[2]))
#realnoise.close()
#contsub=N.array(contsub,'f')
#contsuberr=N.array(contsuberr,'f')
#contsubisoarea=N.array(contsubisoarea,'f')

def makeplot():
    psplotinit("noise.ps")
    
    DATAMIN = 0.
    DATAMAX = 15.
    
    ppgplot.pgbox("",0.0,0,"L",0.0,0)
#print "making graph, ncl = ",ncl
    path=os.getcwd()
    f=path.split('/')
#print path
#print f
    prefix=f[4]
    title=prefix
    ymin=-.05
    ymax=max(aveaperr)+.1
#ymax=10.
    ppgplot.pgenv(DATAMIN,DATAMAX,ymin,ymax,0)
    ppgplot.pglab("linear size N of aperture (pixel)","rms in Sky (ADU/s)",title)
    ppgplot.pgsci(2)#red
    ppgplot.pgslw(4)  #line width
    x=N.sqrt(avearea)
    y=aveaperr
    ppgplot.pgpt(x,y,7)
#errory(x,y,erry)
    ppgplot.pgsci(1)#black
#ppgplot.pgpt(isoarea,fluxerriso,3)
#x1=N.sqrt(contsubisoarea)
#y1=contsuberr
    
#x1=N.sqrt(isoarea)
#y1=fluxerriso
#y=n*y1
    
#ppgplot.pgpt(x1,y1,1)
#ppgplot.pgsci(4)#blue
#ppgplot.pgpt(x1,y,1)
#ppgplot.pgsci(1)#black
    x=N.arange(0,50,1)
    y=x*(a+b*a*x)
#y=N.sqrt(x)*.02
    ppgplot.pgline(x,y)
#errory(x,y,erry)
    
    ppgplot.pgend()

def runit():
    (number,ximage,yimage,fluxiso,fluxerriso,xmax,ymax)=readfile()
    doiraf(ximage,yimage,xmax,ymax)

    (aveap,aveaperr,avearea,aveareaerr) = calcavesky()
    
    (a,b)=solveforab(aveaperr,avearea)
    print "a,b",a,b
    return aveap,aveaperr,avearea,aveareaerr,a,b

#runit()

