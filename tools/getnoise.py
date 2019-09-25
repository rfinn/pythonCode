#!/usr/bin/env python

'''
GOAL:
  measure noise in sky

INPUTS:
  sky-subtracted image and sextractor catalog
  make sure that the readcatalog() function has the right column information
  set npoints and dmin to reasonable parameters
    when it runs ok, set npoints to 1000 and be prepared to wait
  

USEAGE:
  run getnoise(image,catalog)

NOTES
  should update code to use sextractor segmentation image instead of dmin parameter

  3/30/13 - updating to run on newfirm data

'''



import scipy
import mystuff as my
from pyraf import iraf
import numpy as N
import pylab
import pyfits
import random,os


mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

iraf.digiphot()
#iraf.apphot()
iraf.daophot()

npoints=500
nap=15
ncl=1
dmin = 15. #minimum distance to object + sqrt(isoarea)
aplist="1,2,3,4,5,6,7,8,9,10"
nap=10

def runsextractornewfirm(image):
    defaultsexfile='default.sex.newfirm'
    #cp sextractor files to this director
    os.system('cp '+homedir+'research/LocalClusters/sextractor/default.param .')
    os.system('cp '+homedir+'research/LocalClusters/sextractor/default.conv .')
    os.system('cp '+homedir+'research/LocalClusters/sextractor/default.nnw .')
    #s='sex '+self.sex_image+' -c '+homedir+'research/LocalClusters/sextractor/'+defaultsexfile+' -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+self.unc_image
    s='sex '+image+' -c '+homedir+'research/LocalClusters/sextractor/'+defaultsexfile
    os.system(s)
    
	#save test.cat
    #s='cp test.cat '
    #testcat24=self.prefix+'-test.cat'
    #os.rename('test.cat',testcat24)

def runsextractormips(image):
    #cp sextractor files to this director
    os.system('cp '+homedir+'research/LocalClusters/sextractor/default.param .')
    os.system('cp '+homedir+'research/LocalClusters/sextractor/default.nnw .')
    #s='sex '+self.sex_image+' -c '+homedir+'research/LocalClusters/sextractor/'+defaultsexfile+' -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+self.unc_image
    s='sex '+self.sex_image+' -c '+homedir+'research/LocalClusters/sextractor/'+defaultsexfile+' -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+self.unc_image
    os.system(s)
    
	#save test.cat
    s='cp test.cat '
    testcat24=self.prefix+'-test.cat'
    os.rename('test.cat',testcat24)

def calcavesky():
    input=open("noise.dat",'r')
    aperture=[]
    counts=[]
    area=[]
    j=0
    for line in input:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        #if line.find('.fits') > -1: #skip lines with '#' in them
        if (line.find('mosaic') > -1) | (line.find('.fits') > -1): #skip lines with '#' in them
            j=0
            continue
        j=j+1
        if (j > 3):
            t = line.split()
            try:
                aperture.append(float(t[0]))
                counts.append(float(t[1]))
                area.append(float(t[2]))
            except ValueError:
                print t
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
        aveaperr[i]=pylab.std(ap)
        avearea[i]=N.average(aparea)
        aveareaerr[i]=pylab.std(aparea)
        print "ave sky = %8.4f +/- %8.4f" % (N.average(ap),pylab.std(ap))
        print "ave area = %8.4f +/- %8.4f" % (N.average(aparea),pylab.std(aparea))
    return aveap,aveaperr,avearea,aveareaerr

def getpositions(ximage,yimage,isoarea,im):
    #dmin = 15. #minimum distance to object + sqrt(isoarea)

    iraf.imgets(image=im,param='naxis1')#get RA of image
    t=float(iraf.imgets.value)
    xmax=t
    xcenter=t/2.
    iraf.imgets(image=im,param='naxis2')#get RA of image
    t=float(iraf.imgets.value)
    ymax=t
    ycenter=t/2.

    xpos=[]
    ypos=[]
    da = N.sqrt(isoarea)
    d = N.zeros(len(ximage),'f')
    #print len(d)
    i=0
    while i < npoints:
        xtemp=random.uniform(0,1)*xmax
        ytemp=random.uniform(0,1)*ymax
	dcenter=N.sqrt((xtemp-xcenter)**2+(ytemp-ycenter)**2)
        if (dcenter < 500.):
            j=0
	    d = N.sqrt((ximage-xtemp)**2+(yimage-ytemp)**2)
            #for j in range(len(ximage)):
                #print j,len(ximage),len(yimage),len(da),len(d)
                #d[j] = N.sqrt((ximage[j]-xtemp)**2+(yimage[j]-ytemp)**2)-da[j]
            if (min(d) > dmin):       
                xpos.append(xtemp)
                ypos.append(ytemp)
                i=i+1
                #print 'found a good place to measure sky!',i,npoints
    xpos=N.array(xpos,'f')
    ypos=N.array(ypos,'f')
    return xpos,ypos


def getpositionsLCS(ximage,yimage,isoarea,im,coverage_map):
    iraf.imgets(image=im,param='naxis1')#get RA of image
    t=float(iraf.imgets.value)
    xmax=t
    xcenter=t/2.
    iraf.imgets(image=im,param='naxis2')#get RA of image
    t=float(iraf.imgets.value)
    ymax=t
    ycenter=t/2.
    xpos=[]
    ypos=[]
    da = N.sqrt(isoarea)
    d = N.zeros(len(ximage),'f')
    #print len(d)
    i=0
    covimage=pyfits.open(coverage_map)
    cov_data=(covimage[0].data)
    covimage.close()

    while i < npoints:
        xtemp=random.uniform(0,1)*xmax
        ytemp=random.uniform(0,1)*ymax
        # check coverage map to make sure point is on science area
        # coverage map values > 6 are ok
        if cov_data[ytemp,xtemp] > 6:
	    d = N.sqrt((ximage-xtemp)**2+(yimage-ytemp)**2)
            if (min(d) > dmin):       
                xpos.append(xtemp)
                ypos.append(ytemp)
                i=i+1
                #print 'found a good place to measure sky!',i,npoints
    xpos=N.array(xpos,'f')
    ypos=N.array(ypos,'f')
    return xpos,ypos

def measurephot(xpos,ypos,inimage):
    coords=open("noisecoords.dat",'w')
    for i in range(len(xpos)):
        coords.write("%8.1f %8.1f \n" % (xpos[i],ypos[i]))
    sky = open("sky",'w')
    for i in range(npoints):
        sky.write("0.0 \n")
    sky.close()
    aps = open("apertures",'w')
    aps.write(aplist)
    aps.close()
    #os.system("rm noise.dat")
    #print inimage
    #cdir=os.getcwd()
    #print cdir
    #iraf.phot(image=inimage,coords="noisecoords.dat",output="noise.dat",calgorith='none',skyfile="sky",salgori="file",aperture="apertures",interactive="no",verify='no',verbose='yes')
    print 'phot image=',inimage,' coords=noisecoords.dat output=noise.dat calgorith=none skyfile="sky" salgori=file aperture="apertures" interactive="no" verify="no"'
    #iraf.digiphot.daophot.phot(image=self.mosaic24,coords=output_coords,output=datfile,calgorithm='none',skyfile=skyfile,salgori="file",aperture="apertures",interactive="no",verify='no',verbose='no')

def measurephottest(inimage):
    if os.path.exists('noise.dat'):
        os.remove('noise.dat')
    iraf.phot(image=inimage,coords="noisecoords.dat",output="noise.dat",calgorith='none',skyfile="sky",salgori="file",aperture="apertures",interactive="no",verify='no',verbose='yes')
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


def readcatalogold(catalog):
    input=open(catalog,'r')
    number=[]
    magiso=[]
    magerriso=[]
    fluxiso = []
    fluxerriso = []
    ximage = []
    yimage = []
    isoarea = []
    for line in input:
	if line.find('#') > -1: #skip lines with '#' in them
	    continue
	t = line.split()
	number.append(float(t[0]))
	magiso.append(float(t[2]))
	magerriso.append(float(t[3]))
	isoarea.append(float(t[7]))
	fluxiso.append(float(t[8]))
	fluxerriso.append(float(t[9]))
	ximage.append(float(t[10]))
	yimage.append(float(t[11]))
    number=N.array(number,'f')
    magiso=N.array(magiso,'f')
    magerriso=N.array(magerriso,'f')
    fluxiso=N.array(fluxiso,'f')
    fluxerriso=N.array(fluxerriso,'f')
    ximage=N.array(ximage,'f')
    yimage=N.array(yimage,'f')
    isoarea=N.array(isoarea,'f')
    return ximage,yimage,isoarea,fluxerriso

def readcatalogLCS(catalog):
    input=open(catalog,'r')
    number=[]
    magiso=[]
    magerriso=[]
    fluxiso = []
    fluxerriso = []
    ximage = []
    yimage = []
    isoarea = []
    for line in input:
	if line.find('#') > -1: #skip lines with '#' in them
	    continue
	t = line.split()
        try:
            number.append(float(t[0]))
        except ValueError:
            print "couldn't convert to float:", t[0]
	magiso.append(float(t[11]))
	magerriso.append(float(t[12]))
	isoarea.append(float(t[41]))
	fluxiso.append(float(t[9]))
	fluxerriso.append(float(t[10]))
	ximage.append(float(t[1]))
	yimage.append(float(t[2]))
    number=N.array(number,'f')
    magiso=N.array(magiso,'f')
    magerriso=N.array(magerriso,'f')
    fluxiso=N.array(fluxiso,'f')
    fluxerriso=N.array(fluxerriso,'f')
    ximage=N.array(ximage,'f')
    yimage=N.array(yimage,'f')
    isoarea=N.array(isoarea,'f')
    return ximage,yimage,isoarea,fluxerriso

def readcatalognewfirm(catalog):
    input=open(catalog,'r')
    number=[]
    magiso=[]
    magerriso=[]
    fluxiso = []
    fluxerriso = []
    ximage = []
    yimage = []
    isoarea = []
    for line in input:
	if line.find('#') > -1: #skip lines with '#' in them
	    continue
	t = line.split()
        try:
            number.append(float(t[0]))
        except ValueError:
            print "couldn't convert to float:", t[0]
	magiso.append(float(t[11]))
	magerriso.append(float(t[12]))
	isoarea.append(float(t[36]))
	fluxiso.append(float(t[9]))
	fluxerriso.append(float(t[22]))
	ximage.append(float(t[1]))
	yimage.append(float(t[2]))
    number=N.array(number,'f')
    magiso=N.array(magiso,'f')
    magerriso=N.array(magerriso,'f')
    fluxiso=N.array(fluxiso,'f')
    fluxerriso=N.array(fluxerriso,'f')
    ximage=N.array(ximage,'f')
    yimage=N.array(yimage,'f')
    isoarea=N.array(isoarea,'f')
    return ximage,yimage,isoarea,fluxerriso



def getnoiseLCS(image,catalog,cov_image):#sky-subtracted image
    print 'reading catalog'
    (ximage, yimage,isoarea,fluxerriso)=readcatalogLCS(catalog)
    print 'getting positions where I can measure sky'
    (xpos,ypos)=getpositionsLCS(ximage,yimage,isoarea,image,cov_image)
    print 'measuring photometry'
    measurephot(xpos,ypos,image)
    measurephottest(image)
    print 'calculating average sky in each aperture size'
    (aveap,aveaperr,avearea,aveareaerr) = calcavesky()
    print 'solving for a and b'
    (a,b)=solveforab(aveaperr,avearea)
    print "a,b",a,b
    s=str(image)
    (file,post)=s.split('.')
    print 'plotting results and saving figure as ',file
    plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso)
    return a,b
def getnoise(image,catalog):#sky-subtracted image
    print 'reading catalog'
    (ximage, yimage,isoarea,fluxerriso)=readcatalognewfirm(catalog)
    print 'getting positions where I can measure sky'
    (xpos,ypos)=getpositions(ximage,yimage,isoarea,image)
    print 'measuring photometry'
    measurephot(xpos,ypos,image)
    measurephottest(image)
    print 'calculating average sky in each aperture size'
    (aveap,aveaperr,avearea,aveareaerr) = calcavesky()
    print 'solving for a and b'
    (a,b)=solveforab(aveaperr,avearea)
    print "a,b",a,b
    s=str(image)
    (file,post)=s.split('.')
    print 'plotting results and saving figure as ',file
    plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso)
    return a,b

def plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso):
    outfile=str(file)+'noise.ps'
    pylab.figure()
    #print 'printing avearea',avearea
    #print aveaperr
    x=N.sqrt(avearea)#/pylab.pi
    y=aveaperr
    x=N.array(x,'f')
    y=N.array(y,'f')
    pylab.plot(x,y,'bo',label='Aper Sim')
    x1=N.sqrt(isoarea)
    y1=fluxerriso
    pylab.xlabel("linear size N of aperture (pixel)")
    pylab.ylabel("rms in Sky (ADU/s)")
    
    pylab.plot(x1,y1,'r^',label='IsoArea: Err vs sqrt(area)')
    n=2.
    y=n*y1
    #pylab.plot(x1,y,'k^')
    x=N.arange(0,20,1)

    y=x*(a+b*a*x)
    s='Noise Model (a=%5.3f,b=%5.3f)'%(a,b)
    pylab.plot(x,y,color='c',label=s)
    pylab.legend(loc='upper left',numpoints=1)
    pylab.axis([0,25,0,30])
    pylab.savefig(outfile)
    
#(a,b)=getnoise('check.fits','test.cat')

#inimage='MKW11-WCS-mosaic_minus_median_extract.fits'
#incat='MKW11-test.cat'
#covimage='MKW11-WCS-mosaic_cov.fits'
inimage='tu1634749.fits'
covimage='tu1634751.fits'
incat='test.cat'
getnoise(inimage,incat)
