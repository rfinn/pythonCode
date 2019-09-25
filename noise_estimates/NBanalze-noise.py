#!/usr/bin/env python
#useage !!!
# analyze-noise.py 1 j-sky.fits testj.cat
# run              ,
# 1= measure sky noise,
#input image is background subtracted image where background comes from sextractor,
# sextractor catalog
import sys, os
import Numeric as N
import scipy
from math import *
from mystuff import *
import ppgplot
import random
from pyraf import iraf
import random
print sys.argv[0],sys.argv[1]

#npoints=sys.argv[3]
npoints=1000
print "npoints = ",npoints
nap=15
filt = open("filter",'r')
for line in filt:
    filte=int(line)
filt.close()
if (filte == 8):
    ncl = 0
if (filte == 6):
    ncl = 1
if (filte == 4):
    ncl = 2
if (filte == 1):
    ncl = 3
if (filte == 2):
    ncl = 4
if (filte == 1201):
    ncl = 5

if (filte == 1184):
    ncl = 6
print "ncl = ",ncl,filte
def checkposition(x,y):
    if (ncl == 3):#CLJ0023
        pos = checkfilter1(x,y)
    if (ncl == 2):#CL1216
        pos = checkfilter4(x,y)
    if (ncl == 1):#CL1054-12
        pos =checkfilter6(x,y)
    if (ncl == 0):#CL1040
        pos = checkfilter8(x,y)
    if (ncl == 4):#RXJ0152
        pos = checkfilter2(x,y)
    if (ncl == 5):#RXJ0152
        pos = checkfilterms1054(x,y)
    if (ncl == 6):#RXJ0152
        pos = checkfilterHDFN1(x,y)
    if (x < 15):
        pos = 0
    if (y < 15):
        pos = 0
    if (x > (c.xmax[ncl]-15.)):
        pos = 0
    if (y > (c.ymax[ncl]-15.)):
        pos = 0
    return pos

def checkfilter4(x,y):#CL1216
    pos = 1
    #check for bottom right corner
    if ((x > 439) & (y < 431)):
        d=sqrt((x-434)**2 + (y-426)**2)
        if(d > 470) :
            pos = 0
    #check for top left corner
    if ((x < 200) & (y > 590)):
        d=sqrt((x-348)**2 + (y-354)**2)
        if (d > 465):
            pos = 0
    #check for bottom left corner
    if ((x < 100) & (y < 100)):
        d=sqrt((x-405)**2 + (y-343)**2)
        if (d > 470):
            pos = 0
    xystars = [(489.,757.),(75.53,454.1),(492.,386.),(170.,454.1),(181.,336.9),(347.6,615.9),(746.4,496.6),(520.5,499.),(703.1,573.8),(580.1,592.7),(347.6,615.9),(93.,760.),(444.14,139.14)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        if d < 4:
            pos = 0
    return pos
def checkfilter6(x,y):
    pos = 1
    #check for top right and left corners
    if (y > 500):
	d=sqrt((x-430)**2 + (y-423)**2)
        if(d > 495):
            pos = 0
    #check for bottom left corner
    if ((y < 250)):
	d=sqrt((x-465)**2 + (y-518)**2)
        if (d > 500):
            pos = 0
    #check for bottom right corner
    if ((x > 514) & (y < 384)):
	d=sqrt((x-514)**2 + (y-384)**2)
        if(d > 400):
            pos = 0
    #bottom edge
    if (y < 10):
        pos = 0
    #reject stars
    xystars = [(498.5,43.3),(100.7,91.5),(716.,125.),(186.9,134.5),(299.,141.),(69.,147.),(389.,219.6),(520.9,235.),(745.7,311.6),(704.,342.),(248.2,422.5),(441.8,484.),(292.,574.),(879.,598.),(203.6,618.7),(349.,620.),(469.6,676.),(676.4,659.7),(766.5,704.),(510.5,828.7),(323.19,634.62),(297.,84.),(827.,653.),(704.,286.),(715.,68.)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 0
    return pos

def checkfilter8(x,y):
    #check for right side
    pos = 1
    if ((x > 334)):
	d=sqrt((x-334)**2 + (y-464)**2)
        if(d > 440):
            pos = 0
    #check for larger inscribed circle
    d=sqrt((x-400)**2 + (y-400)**2)
    if (d > 430.):
        pos = 0
    #reject foreground galaxies that have ap problems & stars
    xystars = [(535.,794.),(505.,312.),(362.4,714.8),(75.53,454.1),(170.,454.1),(170.,666.3),(639.9,417.2),(421.8,189.6),(278.8,797.2),(443.2,715.6),(362.,631.)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        if d < 4:
            pos = 0
    return pos

def checkfilter1(x,y):
    pos = 1
    d=sqrt((x-485)**2 + (y-460)**2)
    if (d > 480.):
	pos=0
    if (y > 810):
        pos = 0
    if (x > 800):
        pos = 0
    #reject stars
    xystars = [(389.27,444.69),(634.71,281.38),(415.,261.)]
    for (xs,ys) in xystars:
        d = sqrt((xs-x)**2 + (ys-y)**2)
        if d < 4:
            pos = 0
    return pos

def checkfilter2(x,y):
    pos = 1
    if (x > 450):
        d=sqrt((x-391)**2 + (y-470)**2)
        if ( d > 510):
            pos = 0

    if (x < 450) & (y < 450):
        d=sqrt((x-485)**2 + (y-528)**2)
        if ( d > 510):
            pos = 0

    if (y > 516):
        d=sqrt((x-481)**2 + (y-418)**2)
        if ( d > 510):
            pos = 0

    if (y < 20):
        pos = 0
    xystars = [(666.,308.),(814.,458.),(420.,151.),(341.,756.)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 0
    return pos


def checkfilterms1054(x,y):
    pos = 1
    d = sqrt((437.-x)**2 + (475.-y)**2)
    if (d > 400):
        pos = 0
    #xystars = [(666.,308.),(814.,458.),(420.,151.),(341.,756.)]
    #for (xs,ys) in xystars:
    #  d = sqrt((xs-x)**2 + (ys-y)**2)
    #  if d < 4:
    #      pos = 0
    return pos

def checkfilterHDFN1(x,y):
    pos = 1
    if (y > 880.):
        pos = 0
    d = sqrt((457.-x)**2 + (484.-y)**2)
    if (d > 440):
        pos = 0
    xystars = [(502.,306.),(458.,306.),(202.,382.),(240.,546.),(194.,346.),(286.,650.)]
    for (xs,ys) in xystars:
      d = sqrt((xs-x)**2 + (ys-y)**2)
      if d < 4:
          pos = 0
    return pos


def getpositions(ximage,yimage,isoarea):
    dmin = 30. #minimum distance to object + sqrt(isoarea)
    xmax=c.xmax[ncl] #max x dimension
    ymax = c.ymax[ncl] #max y dimension
    xpos=[]
    ypos=[]
    da = N.sqrt(isoarea)
    d = N.zeros(len(ximage),'f')
    #print len(d)
    i=0
    while i < npoints:
        xtemp=random.uniform(0,1)*xmax
        ytemp=random.uniform(0,1)*ymax
        if (checkposition(xtemp,ytemp) > 0):
            j=0
            for j in range(len(ximage)):
                #print j,len(ximage),len(yimage),len(da),len(d)
                d[j] = N.sqrt((ximage[j]-xtemp)**2+(yimage[j]-ytemp)**2)-da[j]
            if (min(d) > dmin):       
                xpos.append(xtemp)
                ypos.append(ytemp)
                i=i+1
    xpos=N.array(xpos,'f')
    ypos=N.array(ypos,'f')
    return xpos,ypos

class Cluster:
    def __init__(self):
        #print "dude - a cluster!"
        self.id = [1040,1054,1216,0023,0152]


        self.filt=[1,8,6,4,2,1201,1184]#pisces filter positions
        self.z=N.array([0.704,0.748,0.794,0.845,0.833,.832])
        self.xc=N.array([394.,524.,389.,555.,428.,454.])
        self.yc=N.array([348.,519.,341.,417.,422.,417.])
        self.sig=N.array([418.,504.,1018.,415.,1000.,1170.])
        self.errsigp=N.array([41.,113.,46.,0.1*415.,0.1*self.sig[4],150.])
        self.errsigm=N.array([41.,65.,46.,0.1*415,0.1*self.sig[4],150.])
        self.sfrmin=N.array([0.5,0.5,0.5,0.24,0.24,.5])
        self.dL=N.array([3005.,3236.,3483.,3760.,3694.,3688.]) #Mpc/h
        self.dA=N.array([5.02,5.14,5.25,5.35,5.33,5.33]) #kpc/arcsec/h
        self.flux0=N.array([8.55,8.46,9.54,9.21,12.8,12.8])#flux ZP times 1d-17 in ergs/s/cm^2
        self.fluxmin=N.array([0.17,0.17,0.081,0.11,0.11,.2])#3sig narrow-band flux in adu/s
        #self.xmax = N.array([791.,909.,822.,861.,858.,934.],'f')
        #self.ymax = N.array([821.,872.,769.,861.,908.,910.],'f')
        self.xmax = N.array([791.,909.,822.,861.,973.,934.,938.],'f')
        self.ymax = N.array([821.,872.,769.,861.,973.,910.,901.],'f')

def calcavesky():
    input=open("noise.dat",'r')
    aperture=[]
    counts=[]
    area=[]
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
#    ap1=N.compress((aperture > 0) & (aperture < 2),counts)
#    area1=N.compress((aperture > 0) & (aperture < 2),area)
#    ap2=N.compress((aperture > 1) & (aperture < 3),counts)
#    area2=N.compress((aperture > 1) & (aperture < 3),counts)
#    ap3=N.compress((aperture > 2) & (aperture < 4),counts)
#    area3=N.compress((aperture > 2) & (aperture < 4),counts)
#    ap4=N.compress((aperture > 3) & (aperture < 5),counts)
#    area4=N.compress((aperture > 3) & (aperture < 5),counts)
#    ap5=N.compress((aperture > 4) & (aperture < 6),counts)
#    area5=N.compress((aperture > 4) & (aperture < 6),counts)
#    ap6=N.compress((aperture > 5) & (aperture < 7),counts)
#    area6=N.compress((aperture > 5) & (aperture < 7),counts)
#    ap7=N.compress((aperture > 6) & (aperture < 8),counts)
#    area7=N.compress((aperture > 6) & (aperture < 8),counts)
#    ap8=N.compress((aperture > 7) & (aperture < 9),counts)
#    area8=N.compress((aperture > 7) & (aperture < 9),counts)
#    ap9=N.compress((aperture > 8) & (aperture < 10),counts)
#    area9=N.compress((aperture > 8) & (aperture < 10),counts)
#    ap10=N.compress((aperture > 9) & (aperture < 11),counts)
#    area10=N.compress((aperture > 9) & (aperture < 11),counts)

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
c = Cluster()

input=open(sys.argv[3],'r')
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
#ncl=0
#iraf=0
doiraf=int(sys.argv[1])
if (doiraf > 0):
    (xpos,ypos) = getpositions(ximage,yimage,isoarea)
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
    image = sys.argv[2]
    print image
    #image="j-sky.fits"
    #image="n-sky.fits"
    iraf.digiphot.daophot.phot(image,coords="noisecoords.dat",output="noise.dat",skyfile="sky",salgori="file",aperture="apertures",interactive="no")

(aveap,aveaperr,avearea,aveareaerr) = calcavesky()

(a,b)=solveforab(aveaperr,avearea)
print "a,b",a,b
contsub=[]
contsuberr=[]
contsubisoarea=[]
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
psplotinit("noise.ps")

DATAMIN = 0.
DATAMAX = 30.

ppgplot.pgbox("",0.0,0,"L",0.0,0)
print "making graph, ncl = ",ncl
if ncl == 0:
    title="CL1040"
    #ymin=-.05
    #ymax=.8
    #a=.25
    #b=.03*a
    #n=3.

    #a=.15
    #b=.1*a
    n=2.
if ncl == 1:
    title="CL1054"
    #a=.015
    #b=.15*a
    n=1.7
if ncl == 2:
    title="CL1216"
    #a=.008
    #b=.165*a
    n=1.5

if ncl == 3:
    title="CLJ0023"
    #a=.008
    #b=.165*a
    n=2.5

if ncl == 4:
    title="RXJ0152"
    #a=.008
    #b=.165*a
    n=1.3

if ncl == 5:
    title="MS1054"
    #a=.008
    #b=.165*a
    n=1.45
if ncl == 6:
    title="HDFN1"
    #a=.008
    #b=.165*a
    n=1.45
ymin=-.05
ymax=max(aveaperr)
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
x1=N.sqrt(isoarea)
y1=fluxerriso
y=n*y1

ppgplot.pgpt(x1,y1,1)
ppgplot.pgsci(4)#blue
ppgplot.pgpt(x1,y,1)
ppgplot.pgsci(1)#black
x=N.arange(0,50,1)
y=x*(a+b*a*x)
#y=N.sqrt(x)*.02
ppgplot.pgline(x,y)
#errory(x,y,erry)

ppgplot.pgend()
