#!/scisoft/i386/bin/python
import scipy
import mystuff as my
from pyraf import iraf
import Numeric as N
import pylab
import random,os
npoints=1000
nap=15
ncl=1
def calcavesky():
    input=open("noise.dat",'r')
    aperture=[]
    counts=[]
    area=[]
    #j=0
    for line in input:
        if line.find('#') > -1: #skip lines with '#' in them
            continue
        if line.find('.fits') > -1: #skip lines with '#' in them
            j=0
            continue
        j=j+1
        if (j > 3):
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
        aveaperr[i]=pylab.std(ap)
        avearea[i]=N.average(aparea)
        aveareaerr[i]=pylab.std(aparea)
        print "ave sky = %8.4f +/- %8.4f" % (N.average(ap),pylab.std(ap))
        print "ave area = %8.4f +/- %8.4f" % (N.average(aparea),pylab.std(aparea))
    return aveap,aveaperr,avearea,aveareaerr

def getpositions(ximage,yimage,isoarea,im):
    dmin = 30. #minimum distance to object + sqrt(isoarea)

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
	    d = N.sqrt((ximage-xtemp)**2+(yimage-ytemp)**2)-da
            #for j in range(len(ximage)):
                #print j,len(ximage),len(yimage),len(da),len(d)
                #d[j] = N.sqrt((ximage[j]-xtemp)**2+(yimage[j]-ytemp)**2)-da[j]
            if (min(d) > dmin):       
                xpos.append(xtemp)
                ypos.append(ytemp)
                i=i+1
    xpos=N.array(xpos,'f')
    ypos=N.array(ypos,'f')
    return xpos,ypos

def measurephot(xpos,ypos,image):
    coords=open("noisecoords.dat",'w')
    for i in range(len(xpos)):
        coords.write("%8.1f %8.1f \n" % (xpos[i],ypos[i]))
    sky = open("sky",'w')
    for i in range(npoints):
        sky.write("0.0 \n")
    sky.close()
    aps = open("apertures",'w')
    aps.write("1,2,3,4,5,6,7,8,9,10,11,12,13,14,15")
    aps.close()
    os.system("rm noise.dat")
    iraf.digiphot()
    iraf.daophot()
    print image
    image="j-sky.fits"
    image="n-sky.fits"
    iraf.digiphot.daophot.phot(image,coords="noisecoords.dat",output="noise.dat",skyfile="sky",salgori="file",aperture="apertures",interactive="no",wcsin='logical')

    iraf.apphot()
    iraf.phot('check.fits',coords='noisecoords.dat',output='noise.dat',skyfile='sky',salgori='file',aperture='apertures',interactive='no')

    test=raw_input("enter anything")
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



def getnoise(image,catalog):#sky-subtracted image
    (ximage, yimage,isoarea,fluxerriso)=readcatalog(catalog)
    (xpos,ypos)=getpositions(ximage,yimage,isoarea,image)
    measurephot(xpos,ypos,image)
    iraf.digiphot()
    iraf.daophot()
    iraf.phot('check.fits',coords='noisecoords.dat',output='noise.dat',skyfile='sky',salgori='file',aperture='apertures',interactive='no',calgorith='none',verify='no')

    (aveap,aveaperr,avearea,aveareaerr) = calcavesky()
    (a,b)=solveforab(aveaperr,avearea)
    print "a,b",a,b
    s=str(image)
    (file,post)=s.split('.')
    #plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso)
    return a,b

def plotnoisepylab(aveap,aveaperr,avearea,aveareaerr,a,b,file,isoarea,fluxerriso):
    outfile=str(file)+'noise.ps'
    #print 'printing avearea',avearea
    #print aveaperr
    x=N.sqrt(avearea)
    y=aveaperr
    x=N.array(x,'f')
    y=N.array(y,'f')
    pylab.plot(x,y,'bo')
    x1=N.sqrt(isoarea)
    y1=fluxerriso

    pylab.plot(x1,y1,'r^')
    n=2.
    y=n*y1
    pylab.plot(x1,y,'k^')
    x=N.arange(0,50,1)
    y=x*(a+b*a*x)
    pylab.plot(x,y)
    pylab.savefig(outfile)
    
#(a,b)=getnoise('check.fits','test.cat')

