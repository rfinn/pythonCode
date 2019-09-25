#!/usr/bin/env python

#writing some scripts to aid in the geometric distortion correction 

import Numeric as N
from pyraf import iraf
import glob
import os
iraf.images()
iraf.images.imutil()
iraf.immatch()


RAcenter=132.82483
DECcenter=11.8000
#dimensions of pisces image
imxmax=1024
imymax=1024

def plotsdsspos(im):
    infile=open('opencl.astrom.star21','r')
    dra=[]
    ddec=[]
    iraf.imgets(image=im,param='RA')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    RAcenter=15.*(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    iraf.imgets(image=im,param='DEC')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    DECcenter=(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    for line in infile:
	t=line.split()
	ra=(float(t[0])-RAcenter)*N.cos(N.pi/180.*float(t[1]))#correct delta ra for cos declination
	dec=(DECcenter-float(t[1]))#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	dra.append(ra)#save offsets from center position
	ddec.append(dec)#save offsets from center
    dra=N.array(dra,'f')
    ddec=N.array(ddec,'f')
    infile.close()

    xcenter=512
    ycenter=512
    plate=0.5#arcsec/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/plate
    y=ddec/plate
    #center at (512,512) in pixel coordinates
    x=x+xcenter
    y=y+ycenter
#reflect about y axis and then swap x & y to match the orientation of pisces
    x=imxmax-x#reflect about x
    yt=x#swap x & t
    xt=y
    xt=imxmax-xt
    
    dy=-26.
    dx=20.

    for i in range(1000):
	yt=yt+dy
	xt=xt+dx
	outfile=open('sdsscoords.dat','w')
	for i in range(len(x)):
	    s="%8.2f %8.2f \n"%(xt[i],yt[i])
	    outfile.write(s)
	outfile.close()

	iraf.display(im,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=7)
	iraf.tvmark(2,'sdsscoords.dat',color=204,radii=17)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sdss coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sdss coords\n")
	dy=float(dy)

def plot2masspos(im):
    infile=open('fp_2mass.fp_psc1374.tbl','r')
    dra=[]
    ddec=[]
    iraf.imgets(image=im,param='RA')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    RAcenter=15.*(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    iraf.imgets(image=im,param='DEC')#get RA of image
    t=iraf.imgets.value
    t1=t.split(':')
    for i in range(len(t1)):#convert to floats
	t1[i]=float(t1[i])
    DECcenter=(t1[0]+t1[1]/60.+t1[2]/3600.)#convert to degrees
    for line in infile:
	if line.find('\\') > -1:
	    continue
	if line.find('=') > -1:
	    continue
	if line.find('|') > -1:
	    continue
	t=line.split()
	if (float(t[6]) >15.):
	    continue
	ra=(float(t[1])-RAcenter)*N.cos(N.pi/180.*float(t[2]))#correct delta ra for cos declination
	dec=(DECcenter-float(t[2]))#makes east to the left
	ra=ra*3600.#convert to arcsec
	dec=dec*3600.
	dra.append(ra)#save offsets from center position
	ddec.append(dec)#save offsets from center
    dra=N.array(dra,'f')
    ddec=N.array(ddec,'f')
    infile.close()

    xcenter=512
    ycenter=512
    plate=0.5#arcsec/pixel
    #convert offsets to pixels on pisces at 90"
    x=dra/plate
    y=ddec/plate
    #center at (512,512) in pixel coordinates
    x=x+xcenter
    y=y+ycenter
#reflect about y axis and then swap x & y to match the orientation of pisces
    x=imxmax-x#reflect about x
    yt=x#swap x & t
    xt=y
    xt=imxmax-xt
    
    dy=-26.
    dx=20.

    for i in range(1000):
	yt=yt+dy
	xt=xt+dx
	outfile=open('2masscoords.dat','w')
	for i in range(len(x)):
	    s="%8.2f %8.2f \n"%(xt[i],yt[i])
	    outfile.write(s)
	outfile.close()

	iraf.display(im,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=5)
	iraf.tvmark(2,'2masscoords.dat',color=206,radii=7.)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust positions \n")
	    flag=float(flag)
	if flag > .1:
	    break
	dx=raw_input("enter new x offset for sdss coords\n")
	dx=float(dx)
	dy=raw_input("enter new y offset for sdss coords\n")
	dy=float(dy)

def runsextractor(image):
    s="sex "+str(image)
    for i in range(1000):
	os.system(s)
	os.system('getxycat.pl')
	iraf.display(image,2)
	iraf.tvmark(2,'testxy.cat',color=207,radii=7)
	try:
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=float(flag)
	except ValueError:
	    print "Sorry, didn't get that.  Let's try one more time."
	    flag= raw_input("how does it look?  1=keep, 0=adjust default.sex \n")
	    flag=float(flag)
	if flag > .1:
	    break
	flag= raw_input("Edit default.sex.  Hit any key when ready \n")	
def runall(prefixes,cat):
    for prefix in prefixes:
	images=glob.glob(prefix)
	iraf.imgets(images[0],'FILTER')
	filter=str(iraf.imgets.value)

	for im in images:
	    runsextractor(im)

	    t=im.split('.')
	    if cat < 0.1:
		out='sdssmatch-'+str(t[0])
		refcat='sdsscoords.dat'
		s='cat '+out+' >> sdssmatch'+filter
		plotsdsspos(im)
	    if cat > .1:
		out='2massmatch-'+str(t[0])
		refcat='2masscoords.dat'
		s='cat '+out+' >> 2massmatch'+filter
		plot2masspos(im)
	    #iraf.xyxymatch(input='testxy.cat',reference='sdsscoords.dat',output=out,tolerance=20.,refpoints='refpoints',interactive='yes')
	    iraf.xyxymatch(input='testxy.cat',reference=refcat,output=out,tolerance=10.,refpoints='refpoints',interactive='no')
	    os.system(s)

	infile='sdssmatch'+filter
	infile='2massmatch'+filter
	iraf.geomap(input=infile,database='PiscesBok',transform=filter,xmin=1.,xmax=1024,ymin=1.,ymax=1024.)

def testone():
    im='gmqopen0004.fits'
    runsextractor(im)
    plot2masspos(im)

def detectfield():
    images=glob.glob('EGSs-*.fits')
    catalogs=[]
    for im in images:
	runsextractor(im)
	(name,junk)=im.split('.')
	(pre,filt)=name.split('-')
	s='cp default.sex default.sex.'+str(filt)
	os.system(s)
	s='cp test.cat test.'+str(filt)'.cat'
	os.system(s)
	s='test.'+str(filt)+'.cat'
	catalogs.append(s)

    return catalogs


class galaxy:
    def __init__(self):#individual galaxy properties
        #print "dude - a galaxy!"
        self.magj = []
        self.errmagj = []
        self.isoarea = []
        self.fj = []
        self.errfj = []
        self.x = []
        self.y = []
        self.fwhm = []
        self.flag = []
        self.seclass  = []
        self.f1j = []
        self.f2j = []
        self.f3j = []
        self.f4j = []
        self.f5j = []
        self.f6j = []
        self.f7j = []
        self.errf1j = []
        self.errf2j = []
        self.errf3j = []
        self.errf4j = []
        self.errf5j = []
        self.errf6j = []
        self.errf7j = []
        self.mag6j = []
        self.errmag6j = []

    def readcat(self):
        name="testj.cat"
        j=0
        for line in open(name):
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            t=line.split()
            #print t[0]
            #    number 1 magiso 2 magisoerr 3 magauto 6 magautoerr 7 fluxiso 9 fluxisoerr 10 x 11 y 12 magbest 13 magbesterr 14 elong 15 ellip 16 fwhm 17 flag 18 class 19 fluxap1 20 fluxap2 21 fluxap3 22 fluxap1err 23 fluxap2err 24 fluxap3err 25
            self.mag.append(float(t[1]))#magiso
            self.errmag.append(float(t[2]))
            self.isoarea.append(float(t[7]))
            self.f.append(float(t[8]))
            self.errf.append(float(t[9]))
            self.x.append(float(t[10]))
            self.y.append(float(t[11]))
            self.fwhm.append(float(t[16]))
            self.flag.append(float(t[17]))
            self.seclass.append(float(t[18]))
            self.f1.append(float(t[19]))
            self.f2.append(float(t[20]))
            self.f3.append(float(t[21]))
            self.f4.append(float(t[22]))
            self.f5.append(float(t[23]))
            self.f6.append(float(t[24]))
            self.f7.append(float(t[25]))
            self.errf1.append(float(t[26]))
            self.errf2.append(float(t[27]))
            self.errf3.append(float(t[28]))
            self.errf4.append(float(t[39]))
            self.errf5.append(float(t[30]))
            self.errf6.append(float(t[31]))
            self.errf7.append(float(t[32]))
            self.mag6.append(float(t[38])) #mag in 2" aperture
            self.errmag6.append(float(t[45]))

    def conv2array(self):
        self.mag = N.array(self.mag,'f')
        self.errmag = N.array(self.errmag,'f')
        self.isoarea = N.array(self.isoarea,'f')
        self.f = N.array(self.f,'f')
        self.errf = N.array(self.errf,'f')
        #self.errf = N.sqrt(self.isoarea)*c.noisea[ncl]*(1.+c.noiseb[ncl]*N.sqrt(self.isoarea))
        try:
            self.errmag=2.5*N.log10(1. + self.errf/self.f)
        except:
            self.errmag=99.
        self.x = N.array(self.x,'f')
        self.y = N.array(self.y,'f')
        self.fwhm = N.array(self.fwhm,'f')
        self.flag = N.array(self.flag,'f')
        self.seclass  = N.array(self.seclass,'f')
        self.f1 = N.array(self.f1,'f')
        self.f2 = N.array(self.f2,'f')
        self.f3 = N.array(self.f3,'f')
        self.f4 = N.array(self.f4,'f')
        self.f5 = N.array(self.f5,'f')
        self.f6 = N.array(self.f6,'f')
        self.f7 = N.array(self.f7,'f')
        self.errf1 = N.array(self.errf1,'f')
        self.errf2 = N.array(self.errf2,'f')
        self.errf3 = N.array(self.errf3,'f')
        self.errf4 = N.array(self.errf4,'f')
        self.errf5 = N.array(self.errf5,'f')
        self.errf6 = N.array(self.errf6,'f')
        #area=3.1415*(2./.18)**2
        #self.errf6 = N.sqrt(area)*c.noisea[ncl]*(1.+c.noiseb[ncl]*N.sqrt(area))
        self.errf7 = N.array(self.errf7,'f')
        self.mag6 = N.array(self.mag6,'f')
        self.errmag6 = N.array(self.errmag6,'f')
        #self.errmag6=2.5*N.log10(1. + abs(self.errf6/self.f6))


cat=1#0=sdss, 1=2mass
#prefixes=['mqopen3*']#,'mqopen1*','mqopen2*','mqopen3*','mqopen4*','mqopen6*']
#runall(prefixes,cat)

#testone()

catalogs=detectfield()
