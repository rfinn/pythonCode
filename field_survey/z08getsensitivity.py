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
from pylab import *
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
    
    t=image.split('.')
    testcat=t[0]+'-test.cat'
    os.rename('test.cat',testcat)

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
        aveaperr[i]=std(ap)
        avearea[i]=N.average(aparea)
        aveareaerr[i]=std(aparea)
        print "ave sky = %8.4f +/- %8.4f" % (N.average(ap),std(ap))
        print "ave area = %8.4f +/- %8.4f" % (N.average(aparea),std(aparea))
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
    outfile=str(file)+'_noise.png'
    figure()
    #print 'printing avearea',avearea
    #print aveaperr
    x=N.sqrt(avearea)#/pi
    y=aveaperr
    x=N.array(x,'f')
    y=N.array(y,'f')
    plot(x,y,'bo',label='Aper Sim')
    x1=N.sqrt(isoarea)
    y1=fluxerriso
    xlabel("linear size N of aperture (pixel)")
    ylabel("rms in Sky (ADU/s)")
    
    plot(x1,y1,'r^',label='IsoArea: Err vs sqrt(area)')
    n=2.
    y=n*y1
    #plot(x1,y,'k^')
    x=N.arange(0,20,1)

    y=x*(a+b*a*x)
    s='Noise Model (a=%5.3f,b=%5.3f)'%(a,b)
    plot(x,y,color='c',label=s)
    legend(loc='upper left',numpoints=1)
    axis([0,25,0,30])
    title(file)
    savefig(outfile)

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
    t=image.split('.')
    ptitle=t[0]
    plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso)
    return a,b


class testimage:
    def __init__(self,inimage,sexcat,a,b):
        self.image=inimage
        self.scat=sexcat
        t=self.image.split('.')
        self.prefix=t[0]
        self.readsexcat(sexcat)
        self.a=a
        self.b=b
    def readsexcat(self, cat):
        infile=open(cat,'r')
        ngal=0
        for line in infile:
            if line.startswith('#'):
                continue
            ngal += 1
        infile.close()

        infile=open(cat,'r')
        xgal=zeros(ngal,'f')
        ygal=zeros(ngal,'f')
        magbest=zeros(ngal,'f')
        magbesterr=zeros(ngal,'f')
        fluxbest=zeros(ngal,'f')
        fluxbesterr=zeros(ngal,'f')
        magiso=zeros(ngal,'f')
        magisoerr=zeros(ngal,'f')
        fluxiso=zeros(ngal,'f')
        fluxisoerr=zeros(ngal,'f')
        fluxradius=zeros(ngal,'f')
        isoarea=zeros(ngal,'f')
        i=0
        for line in infile:
            if line.startswith('#'):
                continue
            t=line.split()
            xgal[i]=float(t[1])
            ygal[i]=float(t[2])

            fluxiso[i]=float(t[9])
            fluxisoerr[i]=float(t[10])
            magiso[i]=float(t[11])
            magisoerr[i]=float(t[12])

            fluxbest[i]=float(t[21])
            fluxbesterr[i]=float(t[22])
            magbest[i]=float(t[23])
            magbesterr[i]=float(t[24])

            isoarea[i]=float(t[36])
            i += 1
        infile.close()
        self.xgal=xgal
        self.ygal=ygal
        self.fluxbest=fluxbest
        self.fluxbesterr=fluxbesterr
        self.magbest=magbest
        self.magbesterr=magbesterr
        self.isoarea=isoarea

    def plotsnrvsmag(self,ptitle):
        figure()
        snr=abs(self.fluxbest/self.fluxbesterr)
        plot(self.magbest,snr,'k.',label='SE Err',markersize=2)

        a=self.a
        b=self.b
        x=N.sqrt(self.isoarea)
        model_err=x*(a+b*a*x)
        snr=abs(self.fluxbest/model_err)
        plot(self.magbest,snr,'c.',label='Ap Sim Err',markersize=2)
        axis([18,24,1,20])
        axhline(y=10,ls='--',label='Target Depth')
        axvline(x=21.6,ls='--',label='_nolegend_')
        gca().set_yscale('log')
        legend(numpoints=1,loc='lower left')
        xlabel('J MAG_BEST')
        ylabel('Signal-to-Noise Ratio')
        title(ptitle)
        outfile=ptitle+'_snr.png'
        savefig(outfile)
#(a,b)=getnoise('check.fits','test.cat')

#inimage='MKW11-WCS-mosaic_minus_median_extract.fits'
#incat='MKW11-test.cat'
#covimage='MKW11-WCS-mosaic_cov.fits'
inimage='tu1634749.fits'
covimage='tu1634751.fits'
incat='test.cat'
#getnoise(inimage,incat)

images=['tu1634749.fits','tu1635392.fits','tu1636446.fits']
images=['rcl1317-F1187N-image-2013-04-01T05h.fits','rcl1317-JX-image-2013-02-02T10h.fits']
for im in images:
    runsextractornewfirm(im)
    t=im.split('.')
    testcat=t[0]+'-test.cat'
    (a,b)=getnoise(im,testcat)
    cl=testimage(im,testcat,a,b)
    cl.plotsnrvsmag(t[0])
#cl1317=testimage(inimage,incat)
#cl1317.plotsnrvsmag()
