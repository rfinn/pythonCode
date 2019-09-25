#!/usr/bin/env python
"""
useage 
spitzer1.py mode

mode = 0 : calculate quantities, match ediscs, write mastertable
mode = 1 : read master table

"""
import sys, os
import numpy as N
import numpy as num
#import numarray as n
import glob
import scipy
from scipy import stats
from math import *
import pylab
import mystuff as my
#import ppgplot
#import random
import sets
#import steve
import matplotlib

from pyraf import iraf
import pyfits
from scipy.interpolate import splrep,splev,splint
import poisson
import matplotlib.font_manager
import numrecipes as nr


blendall=[]
blendall24=[]
iraf.imcoords()
#frame=int(sys.argv[1])
frame=1
catalogpath='/Users/rfinn/research/clusters/ediscs/catalogs/'
catalogpath='/home/rfinn/research/clusters/ediscs/catalogs/'
delta=3.0 #max allowed offset (in arcseconds) between matched sources
delta=2.0 #max allowed offset (in arcseconds) between matched sources
print "matching tolerance is %3.1f arcseconds"%(delta)
delta=delta/3600.
h=.7#H100
Lsol=3.826e33#normalize by solar luminosity
bellconv=9.8e-11#converts Lir (in L_sun) to SFR/yr
bellconv=4.5e-44#Kenn 98 conversion from erg/s to SFR/yr
#lirmin=N.log10(10**10.7/3.*5)#log10 lir min for comparing all clusters
lirmin=10.91#log10 lir min for comparing all clusters
#lirminloz=10.75#log10 lir min for comparing all clusters
lirminloz=10.75#log10 lir min for comparing all clusters
lirminhiz=10.95#log10 lir min for comparing all clusters
#Mvcut=-19.6
Mvcut=-20.1
evolveMv=1
Mvcuthz=-20.5
minmass=6.1e10
minmassloz=4.24e10
minmasshizblue=3.3e10
minmasslozblue=2.3e10
#minmasshizblue=minmass
#minmasslozblue=minmassloz

#mode=float(sys.argv[1])
mode=1.
matplotlib.rc('legend',numpoints=3,fontsize=20,markerscale=1)

f24IRz=[]
f24IRconv=[]
errf24IRconv=[]
#infile2=open('/home/rfinn/research/clusters/spitzer/MasterTables/fluxconv24toIR.v2.dat','r')
infile2=open('/home/rfinn/research/clusters/spitzer/MasterTables/fluxconv24toIR.dat','r')
for line in infile2:
	if line.find('#') > -1:
		continue
	t=line.split()
	for j in range(len(t)):
		t[j]=float(t[j])
	f24IRz.append(t[0])
	f24IRconv.append(t[6])
	errf24IRconv.append(t[7])
	#f24IRconv.append(t[1])
	#errf24IRconv.append(t[2])
f24IRz=N.array(f24IRz,'f')
f24IRconv=N.array(f24IRconv,'f')
errf24IRconv=N.array(errf24IRconv,'f')
infile2.close()
pylab.cla()
pylab.clf()
pylab.plot(f24IRz,f24IRconv,'k-')
y=f24IRconv+errf24IRconv
pylab.plot(f24IRz,y,'k--')
y=f24IRconv-errf24IRconv
pylab.plot(f24IRz,y,'k--')
pylab.axis([0.2,1.,3.,17.])
pylab.xlabel(r'$\rm z$',fontsize=32)
pylab.ylabel(r'$\rm F(3-1100 \mu m)/\nu F_\nu (24\mu m)$',fontsize=32)
pylab.savefig('Ave24toIRconvvsz.eps')



def calcMvcut(z):
	if evolveMv >0.1:
		mvcut = -20.5 + (0.8-z)#makes Mvcut range from -20.5 at z=0.8 to -20.1 at z=0.4
	else:
		mvcut=-20.1
	return mvcut

def getsfrfromlir(lir):
	sfr=lir*bellconv*Lsol
	return sfr

def yaxissfr(ax1):#when plotting lir, but on log scale
	y1, y2=ax1.get_ylim()
	#print y1,y2
	ax2=pylab.twinx()
	ax2.set_ylim(getsfrfromlir(y1),getsfrfromlir(y2))
	#print 'y limit in SFR =',getsfrfromlir(y1),getsfrfromlir(y2)
	ax2.set_yscale('log')



def yaxissfrlog(ax1):#when plotting log(lir)
	y1, y2=ax1.get_ylim()
	ax2=pylab.twinx()
	ax2.set_ylim(getsfrfromlir(10.**y1),getsfrfromlir(10.**y2))
	#ax2.set_yscale('log')

#def xaxissfrlog(ax1):#when plotting log(lir)
#	y1, y2=ax1.get_xlim()
#	ax2=twiny()
#	ax2.set_xlim(getsfrfromlir(10.**y1),getsfrfromlir(10.**y1))
#	ax2.set_xscale('log')


def xaxissfrlog(ax1):#when plotting log(lir)
	locs, labels = pylab.xticks()
	ax1=pylab.gca()
	ticks=ax1.xaxis.get_major_ticks()
	a=ticks[0]
	xfontsize=a.label1.get_fontsize()
	y1, y2=ax1.get_ylim()
	xmin,xmax,ymin,ymax=ax1.axis()
	ylab=y2+0.03*(y2-y1)
	for x in locs:
	    if x < xmin:
		continue
	    if x > xmax:
		continue
	    x2=getsfrfromlir(10.**x)
	    s='$%3.0f$'%(x2)
	    pylab.text(x,ylab,s,fontsize=xfontsize,horizontalalignment='center')#,transform=ax1.transAxes)


tck24=splrep(f24IRz,f24IRconv)
morphout=open('ToBeClassified.dat','w')
def readLaiMorph():
	infile=open('/home/rfinn/research/clusters/spitzer/lai-morphologies.dat','r')
	Lid=[]#lai id
	Lmorph=[]
	Lvistype=[]
	for line in infile:
		if line.find('#') > -1:
			continue
		t=line.split()#EdiscsID              G.R   R.F  A.L   V.D  Final  Hstvisnumtype
		Lid.append(t[0])
		Lvistype.append(float(t[6]))
		morph=[]
		cols=[1,2,4]#columns corresponding to GR, RF, and VD
		for j in cols:#save all 3 classifications from GR, RF, and VD
			
			s= t[j]#0=E/S0, 1=sp, 2=irr, 3=int, 4=pec, 5=compact
			if s.find('E/S0') > -1:
				morph.append(0)
			elif s.find('sp') > -1:
				morph.append(1)
			elif s.find('irr') > -1:
				morph.append(2)
			elif s.find('int') > -1:
				morph.append(3)
			elif s.find('pec') > -1:
				morph.append(3)
			elif s.find('comp') > -1:
				morph.append(2)
		#Lmorph.append((N.array(morph),'f'))
		Lmorph.append(morph)
	infile.close()

	infile=open('/home/rfinn/research/clusters/spitzer/combine_additionalclassifications.txt','r')
	for line in infile:
		if line.find('#') > -1:
			continue
		t=line.split()#EdiscsID              G.R   R.F  A.L   V.D  Final  Hstvisnumtype
		Lid.append(t[1])
		Lvistype.append(float(t[2]))
		morph=[]
		for i in range(3):#save all 3 classifications from GR, RF, and VD
			j=i+5
			s= t[j]#0=E/S0, 1=sp, 2=irr, 3=int, 4=pec, 5=compact
			if s.find('E/S0') > -1:
				morph.append(0)
			elif s.find('sp') > -1:
				morph.append(1)
			elif s.find('irr') > -1:
				morph.append(2)
			elif s.find('int') > -1:
				morph.append(3)
			elif s.find('pec') > -1:
				morph.append(3)
			elif s.find('comp') > -1:
				morph.append(2)
		#Lmorph.append((N.array(morph,'f')))
		Lmorph.append(morph)
	infile.close()
	#Lmorph=N.array(Lmorph,'f')
	Lvistype=N.array(Lvistype,'f')
	print 'Lai morphologies',len(Lid),len(Lmorph),len(Lvistype)
	#print 'Lai morphologies, shape',N.shape(Lid),N.shape(Lmorph),N.shape(Lvistype)
	#print 'Lai morphologies, Lmorph = ',Lmorph
	return Lid,Lmorph,Lvistype

def findnearestold(x1,y1,x2,y2,delta):
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

def findnearest(x1,y1,x2,y2,delta):#use where command
	matchflag=1
	nmatch=0
	d=N.sqrt((x1-x2)**2 + (y1-y2)**2)#x2 and y2 are arrays
	t=pylab.where(d<delta)
	matches=t[0]
	if len(matches) > 0:
		nmatch=len(matches)
		if nmatch > 1:
			try:
				z=[]
				for im in matches:
					z.append(d[int(im)])
				z=pylab.array(z,'f')
				imatch=matches[pylab.where(z == z.min())]
			except IndexError:
				print len(d),matches,N.array(matches,'i'),d[matches[0]],d[matches[1]],d[matches],d[pylab.array(matches,'i')]
		else:
			imatch=matches[0]
			
	else:
		imatch = 0
		matchflag = 0


	return imatch, matchflag,nmatch


class GemsGalaxies:
    def __init__(self):#use pyfits to read in GEMS data from Eric's file
        print "dude - reading gems data!"

	"""
	Columns of fits table
	['NUMBER',
	'X',
	'Y',
	'RA',
	'DEC',
	'APMR',  Apparent R mag
	'MR',
	'Z',
	'DZ',
	'B',  Absolute B
	'RS',  seems to be all zeros
	'U',   Absolute U
	'V',   Absolute V
	'US',  all zeros
	'GS',
	'U280',
	'COMP',
	'MASS',  stellar mass?
	'FNO',
	'XCTS',
	'XLUM',
	'GEMSID',
	'GEMSRA',
	'GEMSDEC',
	'GEMSMAG',
	'GEMSDMAG',
	'GEMSRE',
	'GEMSDRE',
	'GEMSN',
	'GEMSDN',
	'GEMSAR',
	'WHMIPS', 0 if not on 24um, 1 is 24um detection; 2=non detection
	'SFR',
	'FLUX24',
	'TIR',
	'TUV',
	'BRATE',
	'LOCALDENS',
	'VTYPE']
	
	"""

	infile='/home/rfinn/research/clusters/spitzer/GEMS/for_rose.fits'
	hdulist=pyfits.open(infile)
	tbdata=hdulist[1].data
	z=tbdata.field('Z')#table data
	dz=tbdata.field('DZ')
	flux24=tbdata.field('FLUX24')
	whmips=tbdata.field('WHMIPS')#0 if object not on 24um image; 1 = 24um detection; 2 = detection
	gsfr=tbdata.field('SFR')#24 micron flux in micro Jy
	gtir=tbdata.field('TIR')
	gvtype=tbdata.field('VTYPE')#0=E,1=S0,2=Sa,3=Sbc,4=Sdm,5=Irr,6=Pec/Int,7=comp,X=no image
	#print gvtype
	gmass=tbdata.field('MASS')
	gemsid=tbdata.field('GEMSID')
	gMr=tbdata.field('APMR')#apparent R magnitude
	gMv=tbdata.field('V')#absolute V mag
	gMb=tbdata.field('B')#absolute B mag
	gMu=tbdata.field('U')#absolute B mag
	hdulist.close()

	z2=1.*z
	dz2=1.*dz
	flux242=1.*flux24
	whmips2=1*whmips
	gsfr2=1.*gsfr
	gtir2=1.*gtir
	gvtype2=1.*gvtype
	#print gvtype
	gmass2=1.*gmass
	gemsid2=gemsid
	gMr2=1.*gMr
	gMv2=1.*gMv
	gMb2=1.*gMb
	gMu2=1.*gMu


	print "Got gems galaxies",len(z2),len(whmips2),len(flux242),len(gtir2),len(gvtype2)
	#keep only galaxies w/24 micron detection
	#self.z=N.compress(self.whmips > 0.1,self.z)
	self.z=pylab.compress(whmips > 0.1,z2)
	self.dz=pylab.compress(whmips > 0.1,dz2)
	self.flux24=pylab.compress(whmips > 0.1,flux242)
	self.gsfr=pylab.compress(whmips > 0.1,gsfr2)
	print "len gtir2, whmips2, gsfr,flux24 = ",len(gtir2),len(whmips2),len(self.gsfr),len(self.flux24)
	#gtir=[]
	#for i in range(len(whmips2)):
	#	gt=gtir2[i]
	#	if whmips2[i] > 0.1:
	#		gtir.append(gt[0])
	#self.gtir=N.array(gtir,'d')
	self.gtir=pylab.compress(whmips > 0.1,gtir2)
	#temp=[]
	#for i in range(len(whmips)):
	#	if (whmips[i] > 0.1):
	#		temp.append(gvtype[i])
	#self.gvtype=temp
	self.gvtype=gvtype2[pylab.where(whmips > 0.1)]
	#print self.gvtype
	self.gmass=gmass2[pylab.where(whmips > 0.1)]
	#self.gemsid=N.compress(whmips > 0.1,gemsid2)
	self.gMr=gMr2[pylab.where(whmips > 0.1)]#+5.*pylab.log10(h)
	self.gMv=pylab.compress(whmips > 0.1,gMv2)#+5.*pylab.log10(h)
	self.gMb=pylab.compress(whmips > 0.1,gMb2)#+5.*pylab.log10(h)
	self.gMu=pylab.compress(whmips > 0.1,gMu2)#+5.*pylab.log10(h)
	self.whmips=N.compress(whmips > 0.1,whmips2)


	zmin=0.42
	zmax=0.8
	cutz=1
	if cutz:#keep only galaxies in certain z slice 
		self.flux24=N.compress((self.z > zmin) & (self.z < zmax),self.flux24)
		self.gsfr=N.compress((self.z > zmin) & (self.z < zmax),self.gsfr)
		self.gtir=pylab.compress((self.z > zmin) & (self.z < zmax),self.gtir)
		self.gvtype=self.gvtype[pylab.where(((self.z > zmin) & (self.z < zmax)))]
		self.gmass=N.compress((self.z > zmin) & (self.z < zmax),self.gmass)
		#self.gemsid=N.compress((self.z > zmin) & (self.z < zmax),self.gemsid)
		self.gMr=N.compress((self.z > zmin) & (self.z < zmax),self.gMr)
		self.gMv=N.compress((self.z > zmin) & (self.z < zmax),self.gMv)
		self.gMb=N.compress((self.z > zmin) & (self.z < zmax),self.gMb)
		self.gMu=N.compress((self.z > zmin) & (self.z < zmax),self.gMu)
		self.whmips=N.compress((self.z > zmin) & (self.z < zmax),self.whmips)
		self.dz=N.compress((self.z > zmin) & (self.z < zmax),self.dz)
		self.z=N.compress((self.z > zmin) & (self.z < zmax),self.z)
		self.Mvcut=calcMvcut(self.z)

	#calc stellar mass
	self.stellmass=10.**((5.48-self.gMb)/2.5)*10.**(1.737*(self.gMb-self.gMv)-0.942)#from B band
	#self.stellmass=10.**((4.82-self.gMv)/2.5)*10.**(1.305*(self.gMb-self.gMv)-0.628)#from V band

	#compress visual classifications

	#self.redflag= (self.gMb-self.gMv)>(.022*(self.gMv+20)+.65)#true for red galaxies
	self.redflag= (self.gMu-self.gMv)>(1.15-0.31*self.z-0.08*(self.gMv+20))#true for red galaxies
	#calculate SFR in same way we do for ediscs galaxies
	self.Lir=N.zeros(len(self.z),'d')

	#self.gMrabs=N.zeros(len(self.z),'d')
	for i in range(len(self.z)):
		f24conv=splev(self.z[i],tck24)
		conversion=1.e-6*1.e-23*(4.*N.pi*my.dLcm(self.z[i],h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s
		self.Lir[i]=self.flux24[i]*conversion*f24conv/Lsol
		#self.gMrabs[i]=self.gMr[i] - 5.*pylab.log10(my.dL(self.z[i],h))-25.
		#if i < 20:
		#	print "Gems galaxies Mr",i,self.gMr[i],self.gMrabs[i]
	self.SFRir=self.Lir*bellconv*Lsol

	pylab.clf()
	pylab.hist(self.z,50)
	pylab.savefig('gemszdist.eps')

	zmin=0.42
	zmax=0.59
	if cutz:#keep only galaxies in lowz slice
		self.lzflux24=N.compress((self.z > zmin) & (self.z < zmax),self.flux24)
		self.lzgsfr=N.compress((self.z > zmin) & (self.z < zmax),self.gsfr)
		self.lzgtir=pylab.compress((self.z > zmin) & (self.z < zmax),self.gtir)

		self.lzgvtype=self.gvtype[pylab.where((self.z > zmin) & (self.z < zmax))]
		self.lzgmass=N.compress((self.z > zmin) & (self.z < zmax),self.gmass)
		#self.gemsid=N.compress((self.z > zmin) & (self.z < zmax),self.gemsid)
		self.lzgMr=N.compress((self.z > zmin) & (self.z < zmax),self.gMr)
		self.lzgMv=N.compress((self.z > zmin) & (self.z < zmax),self.gMv)
		self.lzgMb=N.compress((self.z > zmin) & (self.z < zmax),self.gMb)
		self.lzgMu=N.compress((self.z > zmin) & (self.z < zmax),self.gMu)
		self.lzwhmips=N.compress((self.z > zmin) & (self.z < zmax),self.whmips)
		self.lzdz=N.compress((self.z > zmin) & (self.z < zmax),self.dz)
		self.lzz=N.compress((self.z > zmin) & (self.z < zmax),self.z)

		self.lzLir=pylab.compress((self.z > zmin) & (self.z < zmax),self.Lir)
		self.lzSFRir=pylab.compress((self.z > zmin) & (self.z < zmax),self.SFRir)
		self.lzstellmass=pylab.compress((self.z > zmin) & (self.z < zmax),self.stellmass)
		self.lzredflag=pylab.compress((self.z > zmin) & (self.z < zmax),self.redflag)
		self.lzMvcut=calcMvcut(self.lzz)

	zmin=0.63
	zmax=0.8

	if cutz:#keep only galaxies in certain z slice 
		self.hzflux24=N.compress((self.z > zmin) & (self.z < zmax),self.flux24)
		self.hzgsfr=N.compress((self.z > zmin) & (self.z < zmax),self.gsfr)
		self.hzgtir=pylab.compress((self.z > zmin) & (self.z < zmax),self.gtir)
		self.hzgvtype=self.gvtype[pylab.where((self.z > zmin) & (self.z < zmax))]
		self.hzgmass=N.compress((self.z > zmin) & (self.z < zmax),self.gmass)
		#self.gemsid=N.compress((self.z > zmin) & (self.z < zmax),self.gemsid)
		self.hzgMr=N.compress((self.z > zmin) & (self.z < zmax),self.gMr)
		self.hzgMv=N.compress((self.z > zmin) & (self.z < zmax),self.gMv)
		self.hzgMb=N.compress((self.z > zmin) & (self.z < zmax),self.gMb)
		self.hzgMu=N.compress((self.z > zmin) & (self.z < zmax),self.gMu)
		self.hzwhmips=N.compress((self.z > zmin) & (self.z < zmax),self.whmips)
		self.hzdz=N.compress((self.z > zmin) & (self.z < zmax),self.dz)
		self.hzz=N.compress((self.z > zmin) & (self.z < zmax),self.z)
		self.hzLir=pylab.compress((self.z > zmin) & (self.z < zmax),self.Lir)
		self.hzSFRir=pylab.compress((self.z > zmin) & (self.z < zmax),self.SFRir)
		self.hzstellmass=pylab.compress((self.z > zmin) & (self.z < zmax),self.stellmass)
		self.hzredflag=pylab.compress((self.z > zmin) & (self.z < zmax),self.redflag)
		self.hzMvcut=calcMvcut(self.hzz)
		

	temp=N.compress(self.Lir > 10.**lirmin,self.Lir)
	self.avelir = pylab.average(temp)
	self.avelirerr = pylab.std(temp)/pylab.sqrt(1.*len(temp))


	lcomp=10.95
	x=len(N.compress((self.lzgMv < Mvcut) & (self.lzLir > 10.**lcomp),self.lzLir))
	y=len(N.compress((self.lzgMv < Mvcut),self.lzLir))
	(a,b,c)=my.ratioerror(x,y)
	print x,y,a,b,c
	print "Fraction of 0.42 < z < 0.6 GEMS galaxies with LIR > %5.2f, Mv < %5.1f = %5.2f / %5.2f = %5.2f + %5.2f - %5.2f"%(lcomp,Mvcut,x,y,a,b,c)

	x=len(N.compress((self.hzgMv < Mvcut) & (self.hzLir > 10.**lcomp),self.hzLir))
	y=len(N.compress((self.hzgMv < Mvcut),self.hzLir))
	(a,b,c)=my.ratioerror(x,y)
	print "Fraction of 0.6 < z < 0.8 GEMS galaxies with LIR > %5.2f, Mv < -20.1 + (0.8-z) = %5.2f / %5.2f = %5.2f + %5.2f - %5.2f"%(lcomp,x,y,a,b,c)
#	lcomp=10.95
#	x=len(N.compress((self.lzgMv < Mvcut) & (self.lzLir > 10.**lcomp),self.lzLir))
#	y=len(N.compress((self.lzgMv < Mvcut),self.lzLir))
#	(a,b,c)=my.ratioerror(x,y)
#	print x,y,a,b,c
#	print "Fraction of 0.42 < z < 0.6 GEMS galaxies with LIR > %5.2f, Mv < -19.6 = %5.2f / %5.2f = %5.2f + %5.2f - %5.2f"%(lcomp,x,y,a,b,c)

	x=len(N.compress((self.gMv < Mvcut) & (self.Lir > 10.**lcomp),self.Lir))
	y=len(N.compress((self.gMv < Mvcut),self.Lir))
	(a,b,c)=my.ratioerror(x,y)
	print "Fraction of 0.42 < z < 0.8 GEMS galaxies with LIR > %5.2f, Mv < -19.6 = %5.2f / %5.2f = %5.2f + %5.2f - %5.2f"%(lcomp,x,y,a,b,c)

	bv=self.gMb-self.gMv
	x=(N.compress((self.gMv < self.Mvcut) & (self.Lir > 10.**lcomp),bv))
	print "Average B-V of 0.42 < z < 0.8 GEMS galaxies with LIR > %5.2f, Mv < -19.6 = %5.2f + %5.2f"%(lcomp,pylab.average(x),pylab.std(x))
	self.bvave=pylab.average(x)
	self.bvstd=pylab.std(x)/pylab.sqrt(1.*len(x))

	ub=self.gMu-self.gMb
	x=(N.compress((self.gMv < self.Mvcut) & (self.Lir > 10.**lcomp),ub))
	print "Average U-B of 0.42 < z < 0.8 GEMS galaxies with LIR > %5.2f, Mv < -19.6 = %5.2f + %5.2f"%(lcomp,pylab.average(x),pylab.std(x))
	self.ubave=pylab.average(x)
	self.ubstd=pylab.std(x)/pylab.sqrt(1.*len(x))

	print "Number of galaxies in loz GEMS sample = ",len(self.lzLir)
	print "Number of galaxies in hiz GEMS sample = ",len(self.hzLir)

	#print len(self.Lir),len(self.gtir)
	#self.gtir=N.array(self.gtir,'d')
	gtir1=[]
	gtir2=[]
	for i in range(len(self.Lir)):
		c=self.gtir[i]
		#print "testing in gems ",c,self.gtir[i]
		#to run on mac pro
		try:
			gtir1.append(c[0])
			gtir2.append(c[1])
		except IndexError:
			gtir1.append(c)
			gtir2.append(c)
	gtir1=N.array(gtir1,'d')
	gtir2=N.array(gtir2,'d')

# on mac pro, use this
#	for i in range(len(self.Lir)):
#		c=self.gtir[i]
#		#print "testing in gems ",c
#		gtir1.append(c)
#		#gtir2.append(c[1])
	gtir1=N.array(gtir1,'d')
	gtir2=N.array(gtir2,'d')
	pylab.cla()
	pylab.clf()

	pylab.plot(pylab.log10(self.Lir),pylab.log10(gtir1),'bo')
	#pylab.plot(pylab.log10(self.Lir),pylab.log10(gtir2),'ro')
	xl=N.arange(10.,12.5,.01)
	yl=xl
	pylab.plot(xl,yl,'k-')
	r=(gtir1/self.Lir)
	scale=pylab.average(r)
	
	print "scale = ",scale,' gtir1/Lir +/-',pylab.std(r)
	pylab.plot(xl,(yl-0.15),'k--')
	pylab.plot(xl,(yl+pylab.log10(scale)),'g--')
	pylab.xlabel('Our calculation of Lir')
	pylab.ylabel('Gems Lir')
	pylab.savefig('lircomparison.eps')


	pylab.cla()
	pylab.clf()
	pylab.plot(self.SFRir,self.gsfr,'bo')
	xl=N.arange(min(self.SFRir),max(self.SFRir),.1)
	yl=xl
	pylab.plot(xl,yl,'k-',label='y=x')
	pylab.plot(xl,0.45*yl,'k--',label='y=0.45*x')
	pylab.legend(loc='best')
	pylab.xlabel('Our calculation of SFR')
	pylab.ylabel('Gems SFR')
	pylab.savefig('GEMSSFRcomparison.eps')

class ediscs24:
    def __init__(self):#individual galaxy properties
        print "dude - 24 micron data!"
	self.acsflag=acsflag
	self.prefix=prefix
	self.fullprefix=fullprefix
	self.idra=idra
	self.iddec=iddec
	#self.mipsimage='/home/rfinn/research/clusters/spitzer/'+prefix+'/mips/24/'+rvalue+'/ch1/bcd/pbcd/Combine/mosaic.fits'
	self.mipsimage='/home/rfinn/research/clusters/spitzer/final-images/'+prefix+'_final24.fits'
	self.zcl=z
	self.Mvcut = calcMvcut(self.zcl)#-20.5 + (0.8-self.zcl)#makes Mvcut range from -20.5 at z=0.8 to -20.1 at z=0.4
	self.sigma=sigma
	self.apcor4=apcor4
	self.apcor4=1.67
	self.apcor3=apcor3
	self.apcor2=apcor2
	self.f80=f80
	self.errsigmap=errsigmap
	self.errsigmam=errsigmam
	self.rac=float(rac)
	self.decc=float(decc)
	self.DL = my.dL(z,h)
	#self.r200=1.41*self.sigma/1000.*1./N.sqrt(.7+.3*(1+self.zcl)**3)/h
	self.r200=1.73*self.sigma/1000.*1./N.sqrt(.7+.3*(1+self.zcl)**3)/h
	dA=my.DA(self.zcl,h)#kpc/arcsec
	self.r200arcmin=self.r200*1000./dA/60.
	print "R200 (armin) = ",self.r200arcmin,self.sigma
	dzmin=1000.
	for i in range(len(f24IRz)):
		dz=abs(self.zcl-f24IRz[i])
		if dz < dzmin:
			dzmin=dz
			conv=f24IRconv[i]			
			errconv=errf24IRconv[i]
	self.f24IRconv=conv
	self.errf24IRconv=errconv
	print self.prefix," at z= %6.4f: flux 24 to IR conv = %4.1f +/- %4.1f (%4.1f )"%(self.zcl,self.f24IRconv,self.errf24IRconv,self.errf24IRconv/self.f24IRconv*100)

    def readmaster(self):	
	mastertable=self.fullprefix+'mastertable.dat'
	self.readmastertable(mastertable)#comment this line and uncomment line above

    def readmastertable(self,file):
        print file
	nmemb=0.
	nspecmemb=0.
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
	    t=line.split()
	    #if float(t[42]) < 0.5:#keep members only
		#    continue

	
            ngal=ngal+1
        input.close()


	self.specmembflag = N.zeros(ngal,'f')
	self.matchflagspec = N.zeros(ngal,'f')
	self.fflat24 = N.zeros(ngal,'f')
	self.errfflat24 = N.zeros(ngal,'f')
	self.matchflagediscsflat24 = N.zeros(ngal,'f')
	

        self.ediscsID = []
	self.ra = N.zeros(ngal,'f')
	self.dec = N.zeros(ngal,'f')
	self.xcorr = N.zeros(ngal,'f')
	self.ycorr = N.zeros(ngal,'f')
	self.starflag = N.zeros(ngal,'f')
	self.EW = N.zeros(ngal,'f')#Halpha EW
	self.EWerr = N.zeros(ngal,'f')
	self.SFR = N.zeros(ngal,'f')#Halpha SFR
	self.SFRerr = N.zeros(ngal,'f')
	self.matchflagha = N.zeros(ngal,'f')#match flag ha
	self.SFflag = N.zeros(ngal,'f')#if SF 
	self.gimtype = N.zeros(ngal,'f')
	self.matchflagmorphgimtype = N.zeros(ngal,'f')	
	self.vistype = N.zeros(ngal,'f')
	self.matchflagvistype = N.zeros(ngal,'f')
	self.matchflag24 = N.zeros(ngal,'f')
	self.flux24 = N.zeros(ngal,'f')
	self.flux24err = N.zeros(ngal,'f')
	self.nmatchediscs24 = N.zeros(ngal,'f')

	self.misoV = N.zeros(ngal,'f')
	self.misoeVapsim = N.zeros(ngal,'f')
	self.misoR = N.zeros(ngal,'f')
	self.misoeRapsim = N.zeros(ngal,'f')
	self.misoI = N.zeros(ngal,'f')
	self.misoeIapsim = N.zeros(ngal,'f')
	self.misoJ = N.zeros(ngal,'f')
	self.misoeJapsim = N.zeros(ngal,'f')	
	self.misoK = N.zeros(ngal,'f')
	self.misoeKapsim = N.zeros(ngal,'f')
	self.magV = N.zeros(ngal,'f')
	self.mageVapsim = N.zeros(ngal,'f')
	self.magR = N.zeros(ngal,'f')
	self.mageRapsim = N.zeros(ngal,'f')
	self.magI = N.zeros(ngal,'f')
	self.mageIapsim = N.zeros(ngal,'f')
	self.magJ = N.zeros(ngal,'f')
	self.mageJapsim = N.zeros(ngal,'f')
	self.magK = N.zeros(ngal,'f')
	self.mageKapsim = N.zeros(ngal,'f')
	self.membflag = N.zeros(ngal,'f')
	self.newspecmatchflag = N.zeros(ngal,'f')	
	self.defmembflag = N.zeros(ngal,'f')
	self.specz = N.zeros(ngal,'f')
	self.spectype = N.zeros(ngal,'f')
	self.specEWOII = N.zeros(ngal,'f')
	self.matchflagspecediscs = N.zeros(ngal,'f')
	self.specEWOIIflag = N.zeros(ngal,'f')
	self.bestz = N.zeros(ngal,'f')
	self.lowz = N.zeros(ngal,'f')
	self.highz = N.zeros(ngal,'f')
	self.wmin = N.zeros(ngal,'f')
	self.Pclust = N.zeros(ngal,'f')
	self.LUlowzclust = N.zeros(ngal,'f')
	self.LUbestzclust = N.zeros(ngal,'f')
	self.LUhighzclust = N.zeros(ngal,'f')	
	self.LBlowzclust = N.zeros(ngal,'f')
	self.LBbestzclust = N.zeros(ngal,'f')
	self.LBhighzclust = N.zeros(ngal,'f')	
	self.LVlowzclust = N.zeros(ngal,'f')
	self.LVbestzclust = N.zeros(ngal,'f')
	self.LVhighzclust = N.zeros(ngal,'f')	
	self.LRlowzclust = N.zeros(ngal,'f')
	self.LRbestzclust = N.zeros(ngal,'f')
	self.LRhighzclust = N.zeros(ngal,'f')	
	self.LIlowzclust = N.zeros(ngal,'f')
	self.LIbestzclust = N.zeros(ngal,'f')
	self.LIhighzclust = N.zeros(ngal,'f')	
	self.LJlowzclust = N.zeros(ngal,'f')
	self.LJbestzclust = N.zeros(ngal,'f')
	self.LJhighzclust = N.zeros(ngal,'f')	
	self.LKlowzclust = N.zeros(ngal,'f')
	self.LKbestzclust = N.zeros(ngal,'f')
	self.LKhighzclust = N.zeros(ngal,'f')
	self.fluxK = N.zeros(ngal,'f')
	self.UBlowzclust = N.zeros(ngal,'f')
	self.UBbestzclust = N.zeros(ngal,'f')
	self.UBhighzclust = N.zeros(ngal,'f')
	self.BVlowzclust = N.zeros(ngal,'f')
	self.BVbestzclust = N.zeros(ngal,'f')
	self.BVhighzclust = N.zeros(ngal,'f')
	self.UVlowzclust = N.zeros(ngal,'f')
	self.UVbestzclust = N.zeros(ngal,'f')
	self.UVhighzclust = N.zeros(ngal,'f')
	self.matchflagediscsirac = N.zeros(ngal,'f')
	self.iracf1 = N.zeros(ngal,'f')
	self.iracf2 = N.zeros(ngal,'f')
	self.iracf3 = N.zeros(ngal,'f')
	self.iracf4 = N.zeros(ngal,'f')
	self.erriracf1 = N.zeros(ngal,'f')
	self.erriracf2 = N.zeros(ngal,'f')
	self.erriracf3 = N.zeros(ngal,'f')
	self.erriracf4 = N.zeros(ngal,'f')
	self.iracsexflag0 = N.zeros(ngal,'f')
	self.iracsexflag1 = N.zeros(ngal,'f')
	self.iracwch1 = N.zeros(ngal,'f')
	self.iracwch2 = N.zeros(ngal,'f')
	self.iracwch3 = N.zeros(ngal,'f')
	self.iracwch4 = N.zeros(ngal,'f')
	self.iracwmin = N.zeros(ngal,'f')
	self.nmatchediscsirac = N.zeros(ngal,'f')
	self.L24 = N.zeros(ngal,'d')
	self.errL24 = N.zeros(ngal,'d')
	self.LHa = N.zeros(ngal,'d')
	self.errLHa = N.zeros(ngal,'d')
	self.snr24 = N.zeros(ngal,'d')
	self.imagex24 = N.zeros(ngal,'d')
	self.imagey24 = N.zeros(ngal,'d')
	self.fap1 = N.zeros(ngal,'d')#24 um aperture fluxes
	self.fap2 = N.zeros(ngal,'d')
	self.fap3 = N.zeros(ngal,'d')
	self.fap4 = N.zeros(ngal,'d')#2.6 pixel radius
	self.fap5 = N.zeros(ngal,'d')
	self.fap6 = N.zeros(ngal,'d')
	self.fap7 = N.zeros(ngal,'d')
	self.fap8 = N.zeros(ngal,'d')
	self.fap9 = N.zeros(ngal,'d')
	self.fap10 = N.zeros(ngal,'d')
	self.errfap1 = N.zeros(ngal,'d')
	self.errfap2 = N.zeros(ngal,'d')
	self.errfap3 = N.zeros(ngal,'d')
	self.errfap4 = N.zeros(ngal,'d')
	self.errfap5 = N.zeros(ngal,'d')
	self.errfap6 = N.zeros(ngal,'d')
	self.errfap7 = N.zeros(ngal,'d')
	self.errfap8 = N.zeros(ngal,'d')
	self.errfap9 = N.zeros(ngal,'d')
	self.errfap10 = N.zeros(ngal,'d')


	self.LVlow = N.zeros(ngal,'f')
	self.LVbest = N.zeros(ngal,'f')
	self.LVhigh = N.zeros(ngal,'f')	
	self.LRlow = N.zeros(ngal,'f')
	self.LRbest = N.zeros(ngal,'f')
	self.LRhigh = N.zeros(ngal,'f')	
	self.UBlow = N.zeros(ngal,'f')
	self.UBbest = N.zeros(ngal,'f')
	self.UBhigh = N.zeros(ngal,'f')
	self.BVlow = N.zeros(ngal,'f')
	self.BVbest = N.zeros(ngal,'f')
	self.BVhigh = N.zeros(ngal,'f')
	self.UVlow = N.zeros(ngal,'f')
	self.UVbest = N.zeros(ngal,'f')
	self.UVhigh = N.zeros(ngal,'f')

	self.ra24 = N.zeros(ngal,'f')
	self.dec24 = N.zeros(ngal,'f')
	self.flux80flag = N.zeros(ngal,'f')


        input=open(file,'r')
        i=0
	k=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
            t=line.split()
	    #if float(t[42]) < 0.5:#keep members only
		#    continue


	    #print "number of fields in a line = ",len(t)
	    #print i,t[len(t)-1]
	    self.ediscsID.append(t[0])

	    for j in range(1,len(t)):
		t[j]=float(t[j])



	    (self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EW[i],self.EWerr[i],self.SFR[i],self.SFRerr[i],self.matchflagha[i],self.SFflag[i],self.gimtype[i],self.matchflagmorphgimtype[i],self.vistype[i],self.matchflagvistype[i],self.matchflag24[i],self.flux24[i],self.flux24err[i],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.specz[i],self.spectype[i],self.specEWOII[i],self.matchflagspecediscs[i],self.specEWOIIflag[i],self.bestz[i],self.lowz[i],self.highz[i],self.wmin[i],self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i],self.nmatchediscsirac[i],self.L24[i],self.errL24[i],self.LHa[i],self.errLHa[i],self.snr24[i],self.imagex24[i],self.imagey24[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i],self.LVlow[i],  self.LVbest[i],  self.LVhigh[i],  self.LRlow[i],  self.LRbest[i],  self.LRhigh[i], self.UBlow[i],  self.UBbest[i],  self.UBhigh[i],  self.BVlow[i],  self.BVbest[i],  self.BVhigh[i],  self.UVlow[i],  self.UVbest[i],  self.UVhigh[i],self.ra24[i],self.dec24[i],self.flux80flag[i])=t[1:]
	    #print i,self.defmembflag[i],self.matchflagvistype[i],self.vistype[i],self.matchflag24[i]
	    if (self.defmembflag[i] > 0.) & (self.matchflagvistype[i] > 0.):
		k=k+1
		#print k,self.vistype[i],self.matchflag24[i]
	    i=i+1


	input.close()

	self.photmembflag=N.zeros(len(self.membflag),'f')
	for i in range(len(self.membflag)):
	    if (self.membflag[i] > 0.) & (self.wmin[i] > 0.4):
		self.photmembflag[i]=1.
	self.supermembflag=N.zeros(len(self.membflag),'f')
	for i in range(len(self.membflag)):
		if self.newspecmatchflag[i] > 0.:
			self.supermembflag[i]=1.
		elif (self.defmembflag[i] > 0.) & (self.wmin[i] > 0.4):
			self.supermembflag[i]=1.
	    

	pre=self.prefix
	#if pre.find('1216') > -1:#scale GTO data
	#    self.flux24=self.flux24*146.
	#    self.flux24err=self.flux24err*146.
	#    self.L24=self.L24*146.
	#    self.errL24=self.errL24*146.


	self.MR = 4.28 - 2.5*N.log10(self.LRbestzclust)-25.+5.*log10(h)
	self.MU = 5.66 - 2.5*N.log10(self.LUbestzclust)-25.+5.*log10(h)
	#self.errMRlo = self.MR - (4.28 - 2.5*N.log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*N.log10(self.Rlumhi) -25.) - self.MR 
	self.MV = 4.82 - 2.5*N.log10(self.LVbestzclust) -25.+5.*log10(h)
	#self.errMvlo = self.MV - (4.82 - 2.5*N.log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*N.log10(self.vlumhi) - 25. ) - self.MV
	self.MB = 5.48 - 2.5*N.log10(self.LBbestzclust) -25.+5.*log10(h)#check zeropoint here...
	#B-band mass-to-light ratio (Bell et al 2003)
	self.stellmass=10.**(1.737*(self.MB-self.MV)-0.942)*self.LBbestzclust*1.e10/h**2#luminosities are in 10^10 h^-2 Lsol
	#self.stellmass=10.**(1.305*(self.MB-self.MV)-0.628)*self.LVbestzclust*1.e10/h**2#M* from MV (Bell et al 2003)
	self.MRbestz = 4.28 - 2.5*N.log10(self.LRbest)-25.+5.*log10(h)
	#self.errMRlo = self.MR - (4.28 - 2.5*N.log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*N.log10(self.Rlumhi) -25.) - self.MR 
	self.MVbestz = 4.82 - 2.5*N.log10(self.LVbest) -25.+5.*log10(h)
	self.MUbestz = 5.66 - 2.5*N.log10(self.LVbest+self.UVbest)-25.+5.*log10(h)
	#self.errMvlo = self.MV - (4.82 - 2.5*N.log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*N.log10(self.vlumhi) - 25. ) - self.MV

	#self.redflag=(self.MB-self.MV)>(-.022*(self.MV+20)+.65)#true if galaxy is redder than B-V cut
	
	self.redflag= (self.MU-self.MV)>(1.15-0.31*self.zcl-0.08*(self.MV+20))#true for red galaxies

	self.L24=self.L24/Lsol
	self.errL24=self.errL24/Lsol
	#self.f24c=self.fap4*self.apcor4#corrected to infinite aperture
	#self.errf24c=self.errfap4*self.apcor4
	self.f24c=self.flux24
	self.errf24c=self.flux24err
	self.fmin=min(self.f24c)#[pylab.where(abs(self.f24c/self.errf24c)>2.5)])
	self.fmax=max(self.f24c)#[pylab.where(abs(self.f24c/self.errf24c)>2.5)])
	conversion=1.e-6*1.e-23*(4.*N.pi*my.dLcm(self.zcl,h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s
	Lir=self.f24c*conversion*self.f24IRconv#factor of ten converts get scaling to go from from nu fnu(24) to L(IR)
	self.Lir=N.array(Lir,'d')/Lsol#convert to solar luminosities
	errLir=self.errf24c*conversion*self.f24IRconv
	self.errLir=N.array(errLir,'d')/Lsol

	self.Lir80=self.f80*conversion*self.f24IRconv/Lsol#luminosity corresponding to f at 80% comple
	print self.prefix," log10(Lir80) = %6.3f" %(pylab.log10(self.Lir80))


	#calculate Lir for bestz - best photoz

	f24conv=splev(self.bestz,tck24)
	conversion=1.e-6*1.e-23*(4.*N.pi*my.dLcm(self.bestz,h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s
	print len(f24conv),len(conversion)
	Lir=self.f24c*conversion*f24conv#factor of ten converts get scaling to go from from nu fnu(24) to L(IR)
	self.Lirbestz=N.array(Lir,'d')/Lsol#convert to solar luminosities
	errLir=self.errf24c*conversion*self.f24IRconv
	self.errLirbestz=N.array(errLir,'d')/Lsol


	for i in range(len(self.ediscsID)):#set Lir of undetected sources to Lir80, 80% completeness limit
		if (abs(self.matchflag24[i]) < 0.1):
			self.Lir[i]=self.Lir80
			self.errLir[i]=self.Lir80

	sfrconv=bellconv*Lsol
	self.SFRir=self.Lir*sfrconv
	self.SFRirerr=self.errLir*sfrconv
	print "Got ",len(self.membflag)," galaxies"

	self.SFR80=self.Lir80*sfrconv
	print self.prefix," SFR80 = %5.2f" %(self.SFR80)
	#for i in range(len(self.ra)):
	#	if self.matchflag24[i] > 0.:
	#		print i,self.Lir[i],self.errLir[i],self.f24c[i],self.errf24c[i],self.matchflag24[i]
	#print self.Lir

	#create field sample
	#zmin=0.43
	#if self.z > 0.65:
	#	self.fieldz=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.8) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.specz)
	#	self.fieldztot=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.8) & (self.specz > zmin) & (abs(self.matchflag24) < 0.1),self.specz)
	#	self.fieldf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.8) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.fap4)*self.apcor4#convert to inf ap
	#	self.fielderrf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.8) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.errfap4)*self.apcor4#convert to inf ap
	#elif self.z <= 0.65:
	#	self.fieldz=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.6) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.specz)
	#	self.fieldztot=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.6) & (self.specz > zmin) & (abs(self.matchflag24) < 0.1),self.specz)
	#	self.fieldf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.6) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.fap4)*self.apcor4#convert to inf ap
	#	self.fielderrf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.6) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.errfap4)*self.apcor4#convert to inf ap
	#	#self.fieldMr=N.compress((abs(self.newspecmatchflag) < 0.1) & (self.specz < 0.6) & (self.specz > zmin) & (self.matchflag24 > 0.1),self.Mr)

	#self.fieldz=N.compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > 0.1),self.specz)
	self.fieldmatchflag24=N.compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.matchflag24)
	self.fieldz=N.compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1)  & (self.matchflag24 > -0.1),self.specz)
#	self.fieldf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.z) < 0.1) & (self.matchflag24 > 0.1),self.fap4)*self.apcor4#convert to inf ap
#	self.fielderrf24c=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.z) < 0.1) & (self.matchflag24 > 0.1),self.errfap4)*self.apcor4#convert to inf ap
	self.fieldf24c=N.compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.flux24)
	self.fielderrf24c=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.flux24err)
	self.fieldMr=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.MRbestz)
	self.fieldMv=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.MVbestz)
	self.fieldBV=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.BVbest)
	self.fieldUB=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.UBbest)
	self.fieldra=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.ra)
	self.fielddec=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.dec)


	index=N.arange(len(self.ediscsID))
	self.fieldindex=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),index)

	ids=[]
	for i in index:
		ids.append(self.ediscsID[i])
	self.fieldediscsID=ids
	self.fieldLir=N.zeros(len(self.fieldz),'d')
	self.fielderrLir=N.zeros(len(self.fieldz),'d')
	for i in range(len(self.fieldz)):
		f24conv=splev(self.fieldz[i],tck24)
		conversion=1.e-6*1.e-23*(4.*N.pi*my.dLcm(self.fieldz[i],h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s

		self.fieldLir[i]=self.fieldf24c[i]*conversion*f24conv/Lsol
		self.fielderrLir[i]=self.fielderrf24c[i]*conversion*f24conv/Lsol
		dL=my.dL(self.fieldz[i],.7)

	self.fieldSFRir=self.fieldLir*bellconv*Lsol
	self.fielderrSFRir=self.fielderrLir*bellconv*Lsol
	print "Number of spec field galaxies = ",len(self.fieldz)
	#print "w/redshifts = ",self.fieldz


	nmultimatch=N.compress((self.matchflag24 > 0.1) & (self.nmatchediscs24 > 1.1),self.nmatchediscs24)
	n24match=len(N.compress(self.matchflag24 > 0.1,self.matchflag24))
	print "Number of ediscs sources w/more than 1 optical match = ",len(nmultimatch),n24match
	print "number of matches = ",nmultimatch
	self.nmultimatch=len(nmultimatch)
	self.n24match=n24match

	bia='ForBianca/'+self.prefix+'.ediscs.v6.dat'
	outfile=open(bia,'w')
	string='#1-ediscsID  2-photmembflag  3-specmembflag  4-matchflag24 5-log10(Lir erg/s) 6-f24 7-errf24 8-zspec\n'
	outfile.write(string)
	string='#ediscsID: new ediscs ID\n'
	outfile.write(string)
	string='#photmembflag: photmemb and wmin > 0.4 \n'
	outfile.write(string)
	string='#specmembflag: -1=no spec; 0=non-member; 1=member\n'
	outfile.write(string)
	string='#matchflag24: -1=not on 24um image; 0=on image but no detection; 1=24um detection \n'
	outfile.write(string)
	string='#log10(Lir), where Lir is in erg/s and spans 8-1000 um (valid for galaxies with matchflag24 = 1); \n'
	outfile.write(string)
	string='#Lir = 80% completeness limit for non-detections (matchflag24 = 0); \n'
	outfile.write(string)
	string='#Lir = -99 for galaxies not on 24um image (matchflag24 = -1); \n'
	outfile.write(string)
	string='#f24: 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
	outfile.write(string)
	string='#errf24: error in 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
	outfile.write(string)
	string='#f24 = -99 for galaxies not on 24um image (matchflag24 = -1); \n'
	outfile.write(string)
	string='#zspec = spectroscopic redshift; only valid for galaxies with specmembflag > -1\n'
	outfile.write(string)

	for i in range(len(self.Lir)):
		if self.matchflag24[i] > .1:
			string="%s %4.1f %4.1f %4.1f %8.4f %8.1f %8.1f %6.4f \n"%(self.ediscsID[i],self.photmembflag[i],self.newspecmatchflag[i],self.matchflag24[i],N.log10(self.Lir[i]*Lsol),self.flux24[i],self.flux24err[i],self.specz[i])
		elif abs(self.matchflag24[i]) < .1:
			string="%s %4.1f %4.1f %4.1f %8.4f %8.1f %8.1f %6.4f \n"%(self.ediscsID[i],self.photmembflag[i],self.newspecmatchflag[i],self.matchflag24[i],N.log10(self.Lir80*Lsol),self.f80,self.f80,self.specz[i])#put in 80% completeness limit for galaxies not detected on 24um image
		else:
			string="%s %4.1f %4.1f %4.1f %8.4f %8.1f %8.1f %6.4f \n"%(self.ediscsID[i],self.photmembflag[i],self.newspecmatchflag[i],self.matchflag24[i],-99.,-99.,-99.,self.specz[i])#Lir = -99 for galaxies not on 24um image
		outfile.write(string)
	outfile.close()


	file2=str(fullprefix)+'.spec.reg'
	output111=open(file2,'w')
	file3=str(prefix)+'.spec.dat'
	output112=open(file3,'w')
	output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	#output112.write('#ediscsID memb specz spectype RA Dec xcorr ycorr\n')
	output112.write('#ediscsID memb specz spectype RA Dec VisType matchflag24 Lir f24 errf24 SFRIR v errv r errr i erri j errj\n')
	string='#matchflag24: -1=not on 24um image; 0=on image but no detection; 1=24um detection \n'
	output112.write(string)
	string='#log10(Lir), where Lir is in erg/s and spans 8-1000 um (valid for galaxies with matchflag24 = 1); \n'
	output112.write(string)
	string='#Lir = 80% completeness limit for non-detections (matchflag24 = 0); \n'
	output112.write(string)
	string='#Lir = -99 for galaxies not on 24um image (matchflag24 = -1); \n'
	output112.write(string)
	string='#f24: 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
	output112.write(string)
	string='#errf24: error in 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
	output112.write(string)
	string='#f24 = -99 for galaxies not on 24um image (matchflag24 = -1); \n'
	output112.write(string)
	for i in range(len(self.newspecmatchflag)):
		if self.newspecmatchflag[i] > -.5:
			if self.newspecmatchflag[i] < .1:
				string1 = "circle(%12.8f, %12.8f, 2\") \n"%(self.ra[i],self.dec[i])
				output111.write(string1)
				string1 = "%s   0  %9.4f %5.1f %12.8f %12.8f %8.2f %8.2f %8.4e %8.4e %8.4e %8.4e %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i],self.f24c[i],self.errf24c[i],self.SFRir[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i], self.magI[i],self.mageIapsim[i], self.magJ[i],self.mageJapsim[i])
				#string1 = "%s   0  %6.4f %4.1f %12.8f %12.8f\n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i])
				output112.write(string1)
				
			else:
				string1 = "circle(%12.8f, %12.8f, 2\") # color=red \n"%(self.ra[i],self.dec[i])
				output111.write(string1)

				string1 = "%s   1  %9.4f %5.1f %12.8f %12.8f %8.2f %8.2f %8.4e %8.4e %8.4e %8.4e %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i],self.f24c[i],self.errf24c[i],self.SFRir[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i], self.magI[i],self.mageIapsim[i], self.magJ[i],self.mageJapsim[i])
				output112.write(string1)
				string1 = "%s   1  %6.4f %4.1f %12.8f %12.8f %6.1f %6.1f\n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i])


	output111.close()
	output112.close()



	#for patti to investigate mergers.  Including all galaxies that are on 24um image, are on HST image, and have spec
	# will create a ds9 file for each cluster showing members in green, nonmembers in red
	# one data file per clusters that includes various flags and Lir

	if self.acsflag > 0.1:
		file2=str(fullprefix)+'.HST.24.spec.reg'
		output111=open(file2,'w')
		file3=str(prefix)+'.HST.24.spec.dat'
		output112=open(file3,'w')
		file4=str(prefix)+'.HST.24em.spec.reg'
		output113=open(file4,'w')
		file5=str(prefix)+'.HST.24em.spec.dat'
		output114=open(file5,'w')
		output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
		output113.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
		output112.write('#ediscsID memb specz spectype RA Dec VisType matchflag24 Lir f24 errf24 SFRIR v errv r errr i erri j errj\n')
		string='#matchflag24: -1=not on 24um image; 0=on image but no detection; 1=24um detection \n'
		output112.write(string)
		string='#log10(Lir), where Lir is in erg/s and spans 8-1000 um (valid for galaxies with matchflag24 = 1); \n'
		output112.write(string)
		string='#Lir = 80% completeness limit for non-detections (matchflag24 = 0); \n'
		output112.write(string)
		string='#Lir = -99 for galaxies not on 24um image (matchflag24 = -1); \n'
		output112.write(string)
		string='#f24: 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
		output112.write(string)
		string='#errf24: error in 24um flux in microJy (set to 80% completeness limit for non-detections); \n'
		output112.write(string)
		string='#f24 = -99 for galaxies not on 24um image (matchflag24 = -1); \n'


		#create list of galaxies w/24um emission
		for i in range(len(self.matchflag24)):
			if (self.matchflag24[i] > -0.1):
				if (self.supermembflag[i] > .1):
					string1 = "%12.8f %12.8f %s %2i %2i %9.4f %8.2f %8.2f %2i %6.2f\n"%(self.ra[i],self.dec[i], self.ediscsID[i],int(self.matchflag24[i]),int(self.newspecmatchflag[i]),self.specz[i],self.f24c[i],self.SFRir[i],int(self.flux80flag[i]),self.magI[i])
					output114.write(string1)
					if (self.matchflag24[i]>.1):
						string1 = "circle(%12.8f, %12.8f, 2\") # color=red %s\n"%(self.ra[i],self.dec[i], self.ediscsID[i])
						output113.write(string1)
		output113.close()
		output114.close()
		for i in range(len(self.newspecmatchflag)):
			if self.newspecmatchflag[i] > -.5:
				#print "found spec match"
				if (self.matchflag24[i] > -0.1):
					#print "gal on 24 um image, HST flag =", self.matchflagvistype[i]
					if (self.matchflagvistype[i] > 0.1):
						#print "gal on hst image"
						if self.newspecmatchflag[i] < .1:
							string1 = "circle(%12.8f, %12.8f, 2\") \n"%(self.ra[i],self.dec[i])
							output111.write(string1)
							string1 = "%s   0  %9.4f %5.1f %12.8f %12.8f %8.2f %8.2f %8.4e %8.4e %8.4e %8.4e %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i],self.f24c[i],self.errf24c[i],self.SFRir[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i], self.magI[i],self.mageIapsim[i], self.magJ[i],self.mageJapsim[i])
							output112.write(string1)
							
						else:
							string1 = "circle(%12.8f, %12.8f, 2\") # color=red \n"%(self.ra[i],self.dec[i])
							output111.write(string1)
							string1 = "%s   1  %9.4f %5.1f %12.8f %12.8f %8.2f %8.2f %8.4e %8.4e %8.4e %8.4e %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i],self.f24c[i],self.errf24c[i],self.SFRir[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i], self.magI[i],self.mageIapsim[i], self.magJ[i],self.mageJapsim[i])
							#string1 = "%s   1  %6.4f %4.1f %12.8f %12.8f %8.2f %8.2f %8.4e\n"%(self.ediscsID[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i])
							output112.write(string1)

		output111.close()
		output112.close()



	greg='forgreg/'+self.prefix+'.gr.v1.dat'
	outfile=open(greg,'w')
	string='#ediscsID  supermembflag photmembflag  specmembflag  matchflag24 F(24)  errF(24)  log10(Lir erg/s) SFR errSFR\n'
	outfile.write(string)
	string='#ediscsID: new ediscs ID\n'
	outfile.write(string)
	string='#supermembflag: either photmemb or specmemb and wmin > 0.4 \n'
	outfile.write(string)
	string='#photmembflag: photmemb and wmin > 0.4 \n'
	outfile.write(string)
	string='#specmembflag: -1=no spec; 0=non-member; 1=member\n'
	outfile.write(string)
	string='#matchflag24: -1=not on 24um image; 0=on image but no detection; 1=24um detection \n'
	outfile.write(string)
	string='#Total 24um Flux (uJy) \n'
	outfile.write(string)
	string='#error in Total 24um Flux (uJy) \n'
	outfile.write(string)
	string='#log10(Lir), where Lir is in erg/s and spans 8-1000 um \n'
	outfile.write(string)
	string='#SFR Msun/yr  \n'
	outfile.write(string)
	string='#err SFR \n'
	outfile.write(string)
	for i in range(len(self.Lir)):
		string="%s %4.1f %4.1f %4.1f %4.1f %8.4f %8.4f %8.4f %8.4f %8.4f \n"%(self.ediscsID[i],self.supermembflag[i],self.photmembflag[i],self.newspecmatchflag[i],self.matchflag24[i],self.f24c[i],self.errf24c[i],N.log10(self.Lir[i]*Lsol),self.SFRir[i],self.SFRirerr[i])
		outfile.write(string)
	outfile.close()

	self.drdeg=N.sqrt(((self.ra-self.rac)*N.cos(self.dec*N.pi/180.))**2+(self.dec-self.decc)**2)#dr, projected distance from cluster center, in degrees
	self.drarcmin=self.drdeg*60.#dr in arcmin
	self.dr=self.drarcmin/self.r200arcmin#dr/ R200

	self.fielddr=N.compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.dr)

	#calculate fraction of OII members w/no 24 micron
	n=len(N.compress((self.newspecmatchflag > 0.1)& ((self.spectype-3) < 1.1)& (abs(self.matchflag24) < .1),self.matchflag24))
	nt=len(N.compress((self.newspecmatchflag > 0.1)& ((self.spectype-3.)<1.1) & (self.matchflag24 > -.5),self.matchflag24))
	print n,nt, "in calculating OII fraction"
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb OII w/no 24 = %i, Numb OII = %i, frac = %6.3f - %6.3f + %6.3f"%(n,nt,a,b,c)
	self.nOIIno24=n
	self.nOII=nt

	#calculate fraction of 24 micron - detected members w/no OII
	n=len(N.compress((self.newspecmatchflag > 0.1)& ((self.spectype-3.)<1.1) & (self.matchflag24 > .1),self.matchflag24))
	nt=len(N.compress((self.newspecmatchflag > 0.1)& (self.matchflag24 > 0.1),self.matchflag24))
	print "printing n,nt ",n,nt
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb 24 spec memb w/emission lines = %i, Numb 24 = %i, frac = %6.3f - %6.3f +%6.3f"%(n,nt,a,b,c)
	self.n24specnoOII=n
	self.n24spec=nt
	self.nspec=len(N.compress((self.newspecmatchflag > 0.1),self.matchflag24))
	     

	#names=N.compress((self.newspecmatchflag > 0.1)& (self.specEWOIIflag < 0.1)& (self.matchflag24 > .1),self.ediscsID)
	for k in range(len(self.matchflag24)):
		if ((self.newspecmatchflag[k] > 0.1)& (self.specEWOIIflag[k] < 0.1)& (self.matchflag24[k] > .1)):
			s='%s \n'%(self.ediscsID[k])
			vandanaout.write(s)
	#for i in range(50):
	#	print self.prefix, i,"dr = ",self.dr[i],self.rac,self.decc,self.ra[i],self.drdeg[i],self.drarcmin[i],self.r200arcmin
    def readoptIRphot(self,file):#reading in ediscs photometry file
	    input=open(file,'r')
            #get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    ngal=ngal+1
	    input.close()
	    print "Number of galaxies in optIR file = ",ngal
		    
	    self.ediscsID = []
	    xcorr = N.zeros(ngal,'f')#pixel x-position of distortion-corrected (fss_avg) I image
	    ycorr = N.zeros(ngal,'f')#pixel x-position of distortion-corrected (fss_avg) I image
	    self.ra = N.zeros(ngal,'f')
	    self.dec = N.zeros(ngal,'f')
	    self.misoB = N.zeros(ngal,'f')
	    self.misoeBapsim = N.zeros(ngal,'f')
	    self.misoV = N.zeros(ngal,'f')
	    self.misoeVapsim = N.zeros(ngal,'f')
	    self.misoR = N.zeros(ngal,'f')
	    self.misoeRapsim = N.zeros(ngal,'f')
	    self.misoI = N.zeros(ngal,'f')
	    self.misoeIapsim = N.zeros(ngal,'f')
	    self.misoJ = N.zeros(ngal,'f')
	    self.misoeJapsim = N.zeros(ngal,'f')
	    self.misoK = N.zeros(ngal,'f')
	    self.misoeKapsim = N.zeros(ngal,'f')
	    self.EWha =-99.*N.ones(ngal,'f')
	    self.EWhaerr =-99.*N.ones(ngal,'f')
	    self.SFR =-99.*N.ones(ngal,'f')
	    self.SFRerr =-99.*N.ones(ngal,'f')
	    self.lha =-99.*N.ones(ngal,'f')
	    self.lhaerr =-99.*N.ones(ngal,'f')
	    self.matchhaediscs = N.zeros(ngal,'i')
	    self.matchflaghaediscs = N.zeros(ngal,'f')
	    self.mautoB= N.zeros(ngal,'f')
	    self.mautoeBapsim= N.zeros(ngal,'f')
	    self.mautoV= N.zeros(ngal,'f')
	    self.mautoeVapsim= N.zeros(ngal,'f')
	    self.mautoR= N.zeros(ngal,'f')
	    self.mautoeRapsim= N.zeros(ngal,'f')
	    self.mautooI= N.zeros(ngal,'f')
	    self.mautoeIapsim= N.zeros(ngal,'f')
	    self.mautoJ= N.zeros(ngal,'f')
	    self.mautoeJapsim= N.zeros(ngal,'f')
	    self.mautoK= N.zeros(ngal,'f')
	    self.mautoeKapsim= N.zeros(ngal,'f')

	    self.magB = N.zeros(ngal,'f')
            self.mageBapsim = N.zeros(ngal,'f')

	    self.magV = N.zeros(ngal,'f')
            self.mageVapsim = N.zeros(ngal,'f')
            self.magR = N.zeros(ngal,'f')
            self.mageRapsim = N.zeros(ngal,'f')
            self.magI = N.zeros(ngal,'f')
            self.mageIapsim = N.zeros(ngal,'f')
            self.magJ = N.zeros(ngal,'f')
            self.mageJapsim = N.zeros(ngal,'f')
            self.magK = N.zeros(ngal,'f')
            self.mageKapsim = N.zeros(ngal,'f')
	    self.SF_flag = N.zeros(ngal,'f')
	    #initialize other flags
	    self.hstvisnumtype = N.zeros(ngal,'f')
	    self.matchmorphvisediscs = N.zeros(ngal, 'i') 
	    self.matchflagmorphvisediscs = N.zeros(ngal, 'i')
	    self.hstgimtype = N.zeros(ngal,'f')
	    self.matchmorphgimediscs = N.zeros(ngal,'i')
	    self.matchflagmorphgimediscs = N.zeros(ngal,'i')
	    self.fflat24 = N.zeros(ngal,'f')
            self.errfflat24 = N.zeros(ngal,'f')
            self.matchflagediscsflat24 = N.zeros(ngal, 'i')
	    self.matchediscsflat24 = N.zeros(ngal,'i')
	    self.fluxK = N.zeros(ngal,'d')

	    if file.find('vrijk') > -1:

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
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])
			    xcorr[i]=float(t[2])
			    ycorr[i]=float(t[3])

			    self.misoV[i]=float(t[128])
			    self.misoeVapsim[i]=float(t[134])
			    self.misoR[i]=float(t[136])
			    self.misoeRapsim[i]=float(t[142])
			    self.misoI[i]=float(t[144])
			    self.misoeIapsim[i]=float(t[150])
			    self.misoJ[i]=float(t[152])
			    self.misoeJapsim[i]=float(t[158])
			    self.misoK[i]=float(t[160])
			    self.misoeKapsim[i]=float(t[166])
			    self.mautoV[i]=float(t[168])
			    self.mautoeVapsim[i]=float(t[174])
			    self.mautoR[i]=float(t[176])
			    self.mautoeRapsim[i]=float(t[182])
			    self.mautooI[i]=float(t[184])
			    self.mautoeIapsim[i]=float(t[190])
			    self.mautoJ[i]=float(t[192])
			    self.mautoeJapsim[i]=float(t[198])
			    self.mautoK[i]=float(t[200])
			    self.mautoeKapsim[i]=float(t[206])
		    
			    self.magV[i] = float(t[8])
			    self.mageVapsim[i] = float(t[14])
			    self.magR[i] = float(t[16])
			    self.mageRapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magJ[i] = float(t[32])
			    self.mageJapsim[i] = float(t[38])
			    self.magK[i] = float(t[40])
			    self.mageKapsim[i] = float(t[46])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal

	    if file.find('bvik') > -1:

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
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])
			    xcorr[i]=float(t[2])
			    ycorr[i]=float(t[3])

			    self.misoB[i]=float(t[104])
			    self.misoeBapsim[i]=float(t[110])

			    self.misoV[i]=float(t[112])
			    self.misoeVapsim[i]=float(t[118])
			    self.misoI[i]=float(t[120])
			    self.misoeIapsim[i]=float(t[126])
			    self.misoK[i]=float(t[128])
			    self.misoeKapsim[i]=float(t[134])

			    self.mautoB[i]=float(t[136])
			    self.mautoeBapsim[i]=float(t[142])
			    self.mautoV[i]=float(t[144])
			    self.mautoeVapsim[i]=float(t[150])
			    self.mautooI[i]=float(t[152])
			    self.mautoeIapsim[i]=float(t[158])
			    self.mautoK[i]=float(t[160])
			    self.mautoeKapsim[i]=float(t[166])

			    self.magB[i] = float(t[8])
			    self.mageBapsim[i] = float(t[14])
			    self.magV[i] = float(t[16])
			    self.mageVapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magK[i] = float(t[32])
			    self.mageKapsim[i] = float(t[38])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal

	    if file.find('bvijk') > -1:

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
			    self.ediscsID.append(t[0])
			    self.ra[i]=float(t[6])*15.
			    self.dec[i]=float(t[7])
			    xcorr[i]=float(t[2])
			    ycorr[i]=float(t[3])

			    self.misoB[i]=float(t[128])
			    self.misoeBapsim[i]=float(t[134])
			    self.misoV[i]=float(t[136])
			    self.misoeVapsim[i]=float(t[142])
			    self.misoI[i]=float(t[144])
			    self.misoeIapsim[i]=float(t[150])
			    self.misoJ[i]=float(t[152])
			    self.misoeJapsim[i]=float(t[158])
			    self.misoK[i]=float(t[160])
			    self.misoeKapsim[i]=float(t[166])

			    self.mautoB[i]=float(t[168])
			    self.mautoeBapsim[i]=float(t[174])
			    self.mautoV[i]=float(t[176])
			    self.mautoeVapsim[i]=float(t[182])
			    self.mautooI[i]=float(t[184])
			    self.mautoeIapsim[i]=float(t[190])
			    self.mautoJ[i]=float(t[192])
			    self.mautoeJapsim[i]=float(t[198])
			    self.mautoK[i]=float(t[200])
			    self.mautoeKapsim[i]=float(t[206])
		    
			    self.magB[i] = float(t[8])
			    self.mageBapsim[i] = float(t[14])
			    self.magV[i] = float(t[16])
			    self.mageVapsim[i] = float(t[22])
			    self.magI[i] = float(t[24])
			    self.mageIapsim[i] = float(t[30])
			    self.magJ[i] = float(t[32])
			    self.mageJapsim[i] = float(t[38])
			    self.magK[i] = float(t[40])
			    self.mageKapsim[i] = float(t[46])
		    #print self.magK[i]
			    i = i+1
		    input.close()
		    self.nediscs=ngal
	    #convert xcorr and ycorr to ra and dec using wcsctran
	    #final wcs images are in 
	    wcsimage='/home/rfinn/research/clusters/ediscs/images/wcs-final/'+self.fullprefix+'i.fss_avg.fits'
	    outfile=open('wcsin','w')
	    for i in range(len(xcorr)):
		    s="%12.8f %12.8f \n"%(xcorr[i],ycorr[i])
		    outfile.write(s)
	    outfile.close()
	    try:
		    os.system('rm wcsout')
	    except:
		    print "running wcsctran for first time"
	    iraf.wcsctran('wcsin','wcsout',image=wcsimage,inwcs='physical',outwcs='world',verbose='no')
	    wcsra=N.zeros(len(xcorr),'f')
	    wcsdec=N.zeros(len(xcorr),'f')
	    outfile=open('wcsout','r')
	    i=0
	    for line in outfile:
		    f=line.split()
		    wcsra[i]=float(f[0])
		    wcsdec[i]=float(f[1])
		    i=i+1
	    outfile.close()
	    #now set ra and dec to wcsra and wcsdec
	    if self.prefix.startswith('cl1059'):
		    self.ra=wcsra
		    self.dec=wcsdec
	    if self.prefix.startswith('cl1353'):
		    self.ra=wcsra
		    self.dec=wcsdec


	    file=str(self.fullprefix)+'.wcsradec.reg'
	    output111=open(file,'w')
	    output111.write("global color=blue font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	    for i in range(len(wcsra)):
		    string1 = "circle(%12.8f, %12.8f, 3\") \n"%(wcsra[i],wcsdec[i])
		    output111.write(string1)
	    output111.close()



    def measurenearest(self):#measure distances to all other ediscs galaxies
	    #if min < radius of 24um photometry, then set contamflag24 = 1
	    rapphot24=3./3600.#3"radius
	    ngal=len(self.ediscsID)
	    self.contamflag24=N.zeros(ngal,'i')
	    for i in range(ngal):
		    d=N.sqrt((self.ra[i]-self.ra)**2+(self.dec[i]-self.dec)**2)

		    dsort=N.take(d,N.argsort(d))#sort array
		    dmin = dsort[1]#skip first entry (distance of gal to itself)
		    if dmin < rapphot24:
			    self.contamflag24[i]=1
	    contam24=N.compress(self.matchflag24 > .1,self.contamflag24)
	    n24=len(contam24)
	    contam=N.sum(self.contamflag24)*1./(1.*ngal)
	    contam24=N.sum(contam24)*1./(1.*n24)
	    print self.prefix, 'Fraction of galaxies with another ediscs galaxy w/in 24um aperture = %6.4f '%(contam)
	    print self.prefix, 'Fraction of 24um galaxies with another ediscs galaxy w/in 24um aperture = %6.4f '%(contam24)
	    blendall.append(contam)
	    blendall24.append(contam24)
    def read24(self,file,mipsimage):
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


        self.id24 = N.zeros(ngal,'f')
        self.imagex24 = N.zeros(ngal,'f')
        self.imagey24  = N.zeros(ngal,'f')
        self.ra24 = N.zeros(ngal,'f')
        self.dec24 = N.zeros(ngal,'f')
        self.f24 = N.zeros(ngal,'d')#flux
        self.errf24 = N.zeros(ngal,'d')
        self.fap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.fap2 = N.zeros(ngal,'d')#flux
        self.fap3 = N.zeros(ngal,'d')#flux
        self.fap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.fap5 = N.zeros(ngal,'d')#flux
        self.fap6 = N.zeros(ngal,'d')#flux
        self.fap7 = N.zeros(ngal,'d')#flux
        self.fap8 = N.zeros(ngal,'d')#flux
        self.fap9 = N.zeros(ngal,'d')#flux
        self.fap10 = N.zeros(ngal,'d')#flux
        self.errfap1 = N.zeros(ngal,'d')#flux in aperture 1 (1,1.5,2,2.6,3,3.5,4,4.5,5.,5.5) pixels
        self.errfap2 = N.zeros(ngal,'d')#flux
        self.errfap3 = N.zeros(ngal,'d')#flux
        self.errfap4 = N.zeros(ngal,'d')#flux in ap 4 - this is one w/ap cor of 1.67 (Calzetti et al 2007)
        self.errfap5 = N.zeros(ngal,'d')#flux
        self.errfap6 = N.zeros(ngal,'d')#flux
        self.errfap7 = N.zeros(ngal,'d')#flux
        self.errfap8 = N.zeros(ngal,'d')#flux
        self.errfap9 = N.zeros(ngal,'d')#flux
        self.errfap10 = N.zeros(ngal,'d')#flux
        self.snr24 = N.zeros(ngal,'d')#SNR calculated by mopex
        self.ndeblend = N.zeros(ngal,'f')#SNR calculated by mopex


        input=open(file,'r')
        i=0
	print 'reading ',file
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue

            t=line.split()
            #(self.id24[i],self.imagex24[i],self.imagey24[i],self.ra24[i],self.dec24[i],self.f24[i],self.errf24[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
	    #print i,'length of t = ',len(t)
	    #print t
	    #print file,i,line

	    (self.id24[i],self.ndeblend[i],self.ra24[i],self.dec24[i],self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i],self.snr24[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i])=(float(t[0]),float(t[2]),float(t[3]),float(t[5]),float(t[8]),float(t[10]),float(t[13]),float(t[14]),float(t[18]),float(t[25]),float(t[26]),float(t[27]),float(t[28]),float(t[29]),float(t[30]),float(t[31]),float(t[32]),float(t[33]),float(t[34]),float(t[45]),float(t[46]),float(t[47]),float(t[48]),float(t[49]),float(t[50]),float(t[51]),float(t[52]),float(t[53]),float(t[54]))
            i=i+1
        input.close()#44 -> 43
	print 'number of galaxies in 24 micron catalog = ',i,len(self.f24),len(self.fap4)
	#pylab.errorbar(self.f24,self.fap4*1.67,yerr=self.errf24,xerr=self.errfap4*1.67)
	pylab.plot(self.f24,self.fap4*1.67,'bo')
	pylab.plot(self.f24,self.fap4,'ro')
	pylab.plot(self.f24,self.fap3,'go')
	pylab.plot(self.f24,self.fap2,'yo')
	pylab.plot(self.f24,self.fap1,'co')
	x=N.arange(0.,300.,10.)
	y=x
	pylab.plot(x,y,'k-')
	y=2.*x
	pylab.plot(x,y,'k-')
	y=3.*x
	pylab.plot(x,y,'k-')
	y=4.*x
	pylab.plot(x,y,'k-')
	y=5.*x
	pylab.plot(x,y,'k-')
	pylab.xlabel('Mopex F(24)')
	pylab.ylabel('Ap Flux (r=2.6 pix) w/apcor=1.67')
	pylab.axis([0.,50.,0.,50.])
	s=str(self.prefix)+'fluxcomp.eps'
	pylab.savefig(s)
        outfile=open('xy24.dat','w')
        for i in range(len(self.imagex24)):
            #print self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i]
            string="%6.2f %6.2f %8.1f %6.2f \n" % (self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i])
            outfile.write(string)
        outfile.close()

	dL=my.dLcm(z,.7)
	self.l24=self.f24*1.e-6#convert to Jy
	self.l24=self.l24*1.e-23#convert to erg/s/cm^2/Hz
	self.l24=self.l24*4*N.pi*dL**2#convert to erg/s/Hz
	self.l24=self.l24*2.448e12#multiply by bandwidth in Hz to get erg/s
	self.errl24=self.l24*(self.errf24/self.f24)
	#delta=2.5 #max allowed offset (in arcseconds) between matched sources
	#delta=delta/3600.
	self.ra24=self.ra24+mydra24[self.prefix]/3600.
	self.dec24=self.dec24+myddec24[self.prefix]/3600.

	x1=self.ra
	y1=self.dec
	x2=self.ra24
	y2=self.dec24


	
	
	self.matchediscs24 = N.zeros(len(x1), 'i')
	self.matchflag24 = N.zeros(len(x1), 'i')
	self.flux80flag = N.zeros(len(x1), 'i')
	self.nmatchediscs24 = N.zeros(len(x1), 'i')
	self.dediscs24 = N.zeros(len(x1),'f')
	#transform ediscs ra & dec onto mips image, pixel values in order to check if ediscs image lies on 24um image
	outfile=open('wcsin','w')
	for i in range(len(x1)):
		s="%12.8f %12.8f \n"%(x1[i],y1[i])
		outfile.write(s)
	outfile.close()
	try:
		os.system('rm wcsout')
	except:
		print "running wcsctran for first time"
	iraf.wcsctran('wcsin','wcsout',image=mipsimage,inwcs='world',outwcs='physical',verbose='no')
	xpix=N.zeros(len(x1),'f')
	ypix=N.zeros(len(x1),'f')
	outfile=open('wcsout','r')
	i=0
	for line in outfile:
		f=line.split()
		#print i, line, f,len(x1),len(xpix)
		xpix[i]=float(f[0])
		ypix[i]=float(f[1])
		i=i+1
	outfile.close()
	buffer=0.
	xmin=1.
	iraf.imgets(image=mipsimage,param='naxis1')#get RA of image
	xmax=float(iraf.imgets.value)
	ymin=1.
	iraf.imgets(image=mipsimage,param='naxis2')#get RA of image
	ymax=float(iraf.imgets.value)
	nomatchedi24=0
	rashift=0.#alignment shifts to bring 24 and Iband image into better alignment
	decshift=0.
	#if s.find('1059')>-1:
	#	rashift=.0003
	#	decshift=.00028
	
	for i in range(len(x1)):#loop over ediscs sources
		if (xpix[i] < (1.+buffer)):#check to see if transformed coord of ediscs source lies on 24um image
			self.matchflag24[i]=-1.
			continue
		if (xpix[i] > (xmax-buffer)):
			self.matchflag24[i]=-1.
			continue
		if (ypix[i] < (1.+buffer)):
			self.matchflag24[i]=-1.
			continue
		if (ypix[i] > (ymax-buffer)):
			self.matchflag24[i]=-1.
			continue
		s=self.prefix

		#(self.matchediscs24[i],self.matchflag24[i],self.nmatchediscs24[i]) = findnearest((x1[i]+rashift),(y1[i]+rashift),x2,y2,delta)#match ediscs source to 24um source


	#find x and y coords of 24um sources on iband image.  Then check to see if x24 and y24 are greater than min(x-iband), min(y-iband) - this corrects for 24um sources that lie within buffer around usable region of iband image
	outfile=open('wcsin2','w')
	for i in range(len(x2)):#24um ra
		s="%12.8f %12.8f \n"%(x2[i],y2[i])
		outfile.write(s)
	outfile.close()
	try:
		os.system('rm wcsout2')
	except:
		print "running wcsctran for first time"
	s3='/home/rfinn/research/clusters/ediscs/images/wcs-final/'+self.fullprefix+'*i.fss_avg.fits'
	iband=glob.glob(s3)
        #this is not going to work for cl105411 and cl105412
	#12/07/08 - updated to use full prefix, so following two if statements are superfluous
	#if self.prefix.startswith('cl105411'):
	#	s3='/home/rfinn/research/clusters/ediscs/images/wcs-final/cl1054-11*i.fss_avg.fits'
	#	iband=glob.glob(s3)
	#if self.prefix.startswith('cl105412'):
	#	s3='/home/rfinn/research/clusters/ediscs/images/wcs-final/cl1054-12*i.fss_avg.fits'
	#	iband=glob.glob(s3)
	print iband
	#09/18/09 changed outwcs from 'physical' to 'logical'
	iraf.wcsctran('wcsin2','wcsout2',image=iband[0],inwcs='world',outwcs='logical',verbose='no')
	x24pix=N.zeros(len(x2),'f')
	y24pix=N.zeros(len(x2),'f')
	outfile=open('wcsout2','r')
	i=0
	for line in outfile:
		f=line.split()
		#print i, line, f,len(x1),len(xpix)
		x24pix[i]=float(f[0])
		y24pix[i]=float(f[1])
		i=i+1
	outfile.close()

	xminiband=min(self.xcorr)#find max and min x-y coordinate of iband image
	xmaxiband=max(self.xcorr)
	yminiband=min(self.ycorr)
	ymaxiband=max(self.ycorr)
	nomatchedi24=0
	nora=[]
	nodec=[]
	nmulti=0
	multimatches=[]
	totalmatch=0
	for i in range(len(x2)):#loop over 24um positions  x1=ediscs  x2=24um

		if x24pix[i] < xminiband:
			continue
		if x24pix[i] > xmaxiband:
			continue
		if y24pix[i] < yminiband:
			continue
		if y24pix[i] > ymaxiband:
			continue
		(matchid,flag,nmatch) = findnearest(x2[i],y2[i],x1,y1,delta)

			#if self.nmatchediscs24[j] < .1:
		if flag > 0:
			try:
				j=int(matchid)#id of ediscs match
			except TypeError:
				print 'trouble w/matchid = ',matchid
				j=int(matchid[0])
			dnew=N.sqrt((x2[i]-x1[j])**2+(y2[i]-y1[j])**2)
			if self.nmatchediscs24[j] > 0.1:#ediscs galaxy already has a 24um match
				print 'got multiple match to ediscs galaxy!'
				if dnew > self.dediscs24[j]:#new match is not closer to ediscs galaxy
					self.nmatchediscs24[j] += 1
					continue
			self.matchediscs24[j]=i
			self.matchflag24[j]=1
			if (self.fap4[i]*self.apcor4 > self.f80):
				self.flux80flag[j]=1
			self.nmatchediscs24[j] += 1
			self.dediscs24[j]=dnew
			#if ediscs galaxy is already matched to a 24um source, then determine which 24um source is closer
				
		#(self.matchediscs24[i],self.matchflag24[i],self.nmatchediscs24[i]) = findnearest((x1[i]+rashift),(y1[i]+rashift),x2,y2,delta)#match ediscs source to 24um source
		if nmatch > 1:
			nmulti += 1
			multimatches.append(nmatch)
		if flag < .5:
			nomatchedi24 += 1
			nora.append(x2[i])
			nodec.append(y2[i])
		else:
			totalmatch += 1
	nora=N.array(nora,'f')
	nodec=N.array(nodec,'f')
	pylab.cla()
	pylab.clf()
	pylab.plot(nora,nodec,'bo')
	s=prefix+'NoEdiMatch.eps'
	pylab.savefig(s)
	self.nomatchedi24=nomatchedi24
	print self.prefix,": Number of 24 um sources w/ediscs matches = ",totalmatch
	print self.prefix,": Number of 24 um sources w/no ediscs match = ",nomatchedi24
	print self.prefix,": Number of 24 um sources w/multiple ediscs matches = %i (%5.2f per cent)"%(nmulti,float(nmulti)/float(totalmatch)*100)
	print self.prefix,": For multiple matches, number of matches = ",multimatches

	self.nedi24match=totalmatch
	self.nedi24nomatch=nomatchedi24
	self.nedi24multimatch=nmulti

	#print self.prefix,": Number of 24 um sources w/no ediscs match = ",nomatchedi24
	#print self.prefix,": Number of 24 um sources w/multiple ediscs matches = ",nmulti
	#print self.prefix,": For multiple matches, number of matches = ",multimatches

	file2=str(fullprefix)+'.NoEdiMatch.reg'
	output111=open(file2,'w')
	output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	for i in range(len(nora)):
		string1 = "circle(%12.8f, %12.8f, 2\") \n"%(nora[i],nodec[i])
		output111.write(string1)
	output111.close()


    def readirac(self,prefix):
	self.matchediscsirac = N.zeros(len(self.ra), 'i')
	self.matchflagediscsirac = N.zeros(len(self.ra), 'i')
	self.nmatchediscsirac = N.zeros(len(self.ra), 'i')

        print 'Reading IRAC data'
        file='/home/rfinn/research/clusters/spitzer/irac-catalogs/'+str(prefix)+'.irac.cat'
	try:
		input=open(file,'r')
		ngal=0        #get number of galaxies

		for line in input:
			if line.find('#') > -1: #skip lines with '#' in them
				continue
			if line.find('\\') > -1: #skip lines with '#' in them
				continue
			if line.find('|') > -1: #skip lines with '#' in them
				continue
			ngal=ngal+1
		input.close()


		self.iracra = N.zeros(ngal,'f')#irac ra
		self.iracdec = N.zeros(ngal,'f')#irac dec
		self.iracf1 = N.zeros(ngal,'f')#irac flux channel 1
		self.erriracf1 = N.zeros(ngal,'f')
		self.iracf2  = N.zeros(ngal,'f')#irac flux channel 2
		self.erriracf2 = N.zeros(ngal,'f')
		self.iracf3 = N.zeros(ngal,'f')#irac flux channel 3
		self.erriracf3 = N.zeros(ngal,'f')
		self.iracf4 = N.zeros(ngal,'f')#irac flux channel 4
		self.erriracf4 = N.zeros(ngal,'f')
		self.iracsexflag0 = N.zeros(ngal,'f')
		self.iracsexflag1 = N.zeros(ngal,'f')
		self.iracwch1 = N.zeros(ngal,'f')
		self.iracwch2 = N.zeros(ngal,'f')
		self.iracwch3 = N.zeros(ngal,'f')
		self.iracwch4 = N.zeros(ngal,'f')	
		self.iracwmin = N.zeros(ngal,'f')#relative exposure on irac image

		print "IRAC file = ",file
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
			(self.iracra[i],self.iracdec[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i])=(float(t[1]),float(t[2]),float(t[5]),float(t[6]),float(t[7]),float(t[8]),float(t[9]),float(t[10]),float(t[11]),float(t[12]),float(t[13]),float(t[14]),float(t[15]),float(t[16]),float(t[17]),float(t[18]),float(t[19]))
			i=i+1
		input.close()

		x1=self.ra
		y1=self.dec
		x2=self.iracra
		y2=self.iracdec
	

		for i in range(len(x1)):
			(self.matchediscsirac[i],self.matchflagediscsirac[i],self.nmatchediscsirac[i]) = findnearest(x1[i],y1[i],x2,y2,delta)

	except IOError:
		print "No IRAC data yet!"
		ngal=len(self.ra)
		self.iracra = N.zeros(ngal,'f')#irac ra
		self.iracdec = N.zeros(ngal,'f')#irac dec
		self.iracf1 = N.zeros(ngal,'f')#irac flux channel 1
		self.erriracf1 = N.zeros(ngal,'f')
		self.iracf2  = N.zeros(ngal,'f')#irac flux channel 2
		self.erriracf2 = N.zeros(ngal,'f')
		self.iracf3 = N.zeros(ngal,'f')#irac flux channel 3
		self.erriracf3 = N.zeros(ngal,'f')
		self.iracf4 = N.zeros(ngal,'f')#irac flux channel 4
		self.erriracf4 = N.zeros(ngal,'f')
		self.iracsexflag0 = N.zeros(ngal,'f')
		self.iracsexflag1 = N.zeros(ngal,'f')
		self.iracwch1 = N.zeros(ngal,'f')
		self.iracwch2 = N.zeros(ngal,'f')
		self.iracwch3 = N.zeros(ngal,'f')
		self.iracwch4 = N.zeros(ngal,'f')	
		self.iracwmin = N.zeros(ngal,'f')#relative exposure on irac image

	#delta=2. #max allowed offset (in arcseconds) between matched sources


    def readwide24(self,file,file2):
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


        self.idwide24 = N.zeros(ngal,'f')
        self.imagexwide24 = N.zeros(ngal,'f')
        self.imageywide24  = N.zeros(ngal,'f')
        self.rawide24 = N.zeros(ngal,'f')
        self.decwide24 = N.zeros(ngal,'f')
        self.fwide24 = N.zeros(ngal,'f')#flux
        self.errfwide24 = N.zeros(ngal,'f')

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

	    (self.idwide24[i],self.rawide24[i],self.decwide24[i],self.imagexwide24[i],self.imageywide24[i],self.fwide24[i],snr)=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
	    #print "Number of fields in first line = ",len(t)," expecting 5"
	    self.errfwide24[i]=self.fwide24[i]/snr
	    #if i < 10:
	    #	    print line
	#	    for ta in t:
	#		    print ta
	#	    print i, self.fwide24[i]
		
            i=i+1
        input.close()

	delta=1. #max allowed offset (in arcseconds) between matched sources
	delta=delta/3600.
	x1=self.ra24
	y1=self.dec24
	x2=self.rawide24
	y2=self.decwide24
	
	self.matchwide24 = N.zeros(len(x1), 'i')
	self.nmatchwide24 = N.zeros(len(x1), 'i')
	self.matchflagwide24 = N.zeros(len(x1), 'i')
	j=0
	for i in range(len(x1)):
		(self.matchwide24[i],self.matchflagwide24[i],self.nmatchwide24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)
		if self.matchflagwide24[i] > 0:
			if j < 10:
				print j, self.ra24[i],self.dec24[i],self.rawide24[self.matchwide24[i]],self.decwide24[self.matchwide24[i]]
				print j, self.f24[i], self.fwide24[self.matchwide24[i]]
				j=j+1
	i=0
	fw=[]
	fwerr=[]
	f=[]
	ferr=[]
	for i in range(len(self.ra24)):
	    if (self.matchflagwide24[i] > 0):
		f.append(self.f24[i])
		ferr.append(self.errf24[i])
		fw.append(self.fwide24[self.matchwide24[i]])
		fwerr.append(self.fwide24[self.matchwide24[i]])

	f=N.array(f,'d')
	ferr=N.array(ferr,'d')
	fw=N.array(fw,'d')
	fwerr=N.array(fwerr,'d')
	print "average ratio of core/wide flux = ",N.average(f/fw)," +/-",pylab.std(f/fw)
	pylab.plot(f,fw,'y^',label="",markersize=5)
	pylab.xlabel('F(24)-core')
	pylab.ylabel('F(24)-wide')
	#pylab.show()
	pylab.clf()
	pylab.cla()
	pylab.plot(f,fw,'b^',label="Scan vs Phot",markersize=5)

	pylab.xlabel('F(24)-core (uJy)')
	pylab.ylabel('F(24)-wide (uJy)')
	x=N.arange(0.,(max(f)+100.),100.)
	y=x
	pylab.plot(x,y,'k',label='y=x')
	y=0.7*x
	pylab.plot(x,y,'k--',label="y=0.7 x")
	t=pylab.legend(loc='upper left')				
	pylab.axis([100.,300.,100.,300.])
	d=N.compress(f < 200.,f)
	w=N.compress(f < 200.,fw)
	r=(w-d)/d
	print "STD for f < 250 uJy (wide-deep)/deep = ",pylab.std(r)
	#pylab.errorbar(f,fw,yerr=fwerr,xerr=ferr,fmt=None,ecolor='k')
	pylab.savefig(file2)

    def magnitudetofluxconversion(self,prefix,plothaflag):
	    i=0
	    pylab.cla()
	    pylab.clf()
	    log=open('Klog','w')
	    for i in range(len(self.magK)):
		    self.fluxK[i] = 640/10**(self.magK[i]/2.5)
		    print >>log, i, self.magK[i],self.fluxK[i]
		    i=i+1
	    log.close()
	    ra=[]
	    dec=[]
	    xbad = N.compress(self.magK < -90.,self.ra)
	    ybad = N.compress(self.magK < -90.,self.dec)
	    pylab.plot(self.ra,self.dec,'ko',label="full phot",markersize=3)
	    pylab.plot(xbad,ybad,'ro',label="No NIR",markersize=3)
	    pylab.plot(self.ra24,self.dec24,'bs',label="24um",markersize=7)
	    if (plothaflag > 0.1):
		    pylab.plot(self.raha,self.decha,'c.',label="Halpha",markersize=3)

	    i=0
	    for i in range(len(self.ra)):
		    if (self.matchflaghaediscs[i] > 0):
			    if (self.matchflag24[i] > 0):
				    ra.append(self.ra[i])
				    dec.append(self.dec[i])
				    #print i,"got a match"
		    i=i+1

	    ra=N.array(ra,'d')
	    dec=N.array(dec,'d')
	    pylab.plot(ra,dec,'r^',label="24-Ha-Edi",markersize=9)
	    t=pylab.legend(loc='upper left')					       
	    outfile=str(prefix)+'.eps'
	    pylab.savefig(outfile)
		    
    def readsextractor(self,file):#read sextractor output
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

        self.x2 = N.zeros(ngal,'f')
        self.y2  = N.zeros(ngal,'f')
        self.f2 = N.zeros(ngal,'f')#flux
        self.errf2 = N.zeros(ngal,'f')

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
	    self.f2[i]=float(t[8])
	    self.errf2[i]=float(t[9])
	    #self.f2[i]=float(t[35])
	    #self.errf2[i]=float(t[42])
	    self.x2[i]=float(t[10])
	    self.y2[i]=float(t[11])	

            i=i+1
        input.close()


    def readhafile(self,file):#read halpha output and match to ediscs
	input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()
		
	self.haID = []
	self.raha = N.zeros(ngal,'f')
	self.decha = N.zeros(ngal,'f')
	#self.Fluxn = N.zeros(ngal,'f')
	#self.Fluxnerr = N.zeros(ngal,'f')
	#self.FluxJ = N.zeros(ngal,'f')
	#self.FluxJerr = N.zeros(ngal,'f')
	#self.magJ = N.zeros(ngal,'f')
	#self.magJerr = N.zeros(ngal,'f')
	#self.ContSub = N.zeros(ngal,'f')
	#self.ContSuberr = N.zeros(ngal,'f')
	#self.SNR = N.zeros(ngal,'f')
	self.EWha = N.zeros(ngal,'f')
	self.EWhaerr = N.zeros(ngal,'f')
	self.SFR = N.zeros(ngal,'d')
	self.SFRerr = N.zeros(ngal,'d')
	self.lha = N.zeros(ngal,'d')
	self.lhaerr = N.zeros(ngal,'d')
	self.SF_flag = N.zeros(ngal,'f')

	input=open(file,'r')
	#get number of galaxies
	i=0
	for line in input:
		if line.find('#') > -1:
			continue
		t=line.split()
		shortname=t[0]
		self.haID.append(shortname)
		self.raha[i]=((float(shortname[0:2])+float(shortname[2:4])/60.+(float(shortname[4:7])/10.)/3660.)*15.)#convert from hours to degrees
		
		if float(shortname[7:10]) < 0:
			self.decha[i]=float(shortname[7:10])-float(shortname[10:12])/60.-(float(shortname[12:])/10.-0.)/3660#to align w/mips
		else:
			self.decha=float(shortname[7:10])+float(shortname[10:12])/60.+float(shortname[12:])/10./3660
		self.SFR[i]=float(t[14])
		self.SFRerr[i]=float(t[15])
		self.EWha[i]=float(t[12])
		self.EWhaerr[i] = float(t[13])
		self.SF_flag[i] = float(t[16])
		i=i+1
	input.close()
       
	self.matchhaediscs = N.zeros(len(self.ediscsID), 'i')
	self.matchflaghaediscs = -1.*N.ones(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.haID)):
			haname=str(self.haID[j])
			if ediscsname.find(haname) > -1:
				#print "matched Halpha name to ediscs ID \n",(self.raha[j]-self.ra[i])*3600./15.,(self.decha[j]-self.dec[i])*3600.
				self.matchhaediscs[i]=j
				self.matchflaghaediscs[i]=1
				self.raha[j]=self.ra[i]
				self.decha[j]=self.dec[i]
				if (self.EWha[j] < 0.):
					self.matchflaghaediscs[i]=0.
				break
	self.lha=self.SFR/(7.9e-42)
	self.lhaerr=self.SFRerr/7.9e-42


    def readha2file(self,file):#read h alpha output of cl0023
	    input=open(file,'r')
	#get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    ngal=ngal+1
		    input.close()
		    racenter=(0.+23.0/60.+51.8/60./60.)*15.
		    deccenter=4.+22./60.+41./60./60.
		    self.raha = N.zeros(ngal,'f')
		    self.decha = N.zeros(ngal,'f')
	#self.x = N.zeros(ngal,'f')
	#self.y = N.zeros(ngal,'f')
	#self.Fluxn = N.zeros(ngal,'f')
	#self.Fluxnerr = N.zeros(ngal,'f')
	#self.FluxJ = N.zeros(ngal,'f')
	#self.FluxJerr = N.zeros(ngal,'f')
	#self.Contsub = N.zeros(ngal,'f')
	#self.Contsuberr = N.zeros(ngal,'f')
	#self.EW(Ha) = N.zeros(ngal,'f')
	#self.EW(Ha)err = N.zeros(ngal,'f')
	#self.Flux(Ha) =  N.zeros(ngal,'f')
	#self.Flux(Ha)err = N.zeros(ngal,'f')
	#self.L(Ha) = N.zeros(ngal,'f')
	#self.L(Ha)err = N.zeros(ngal,'f')
		    self.SFR = N.zeros(ngal,'f')
		    self.SFRerr = N.zeros(ngal,'f')
		    
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
		    dra=float(t[1])/60./60.
		    dec=float(t[2])/60./60.
		    self.raha[i]=racenter+dra
		    self.decha[i]=deccenter+dec
		    self.SFR[i]=float(t[17])
		    self.SFRerr[i]=float(t[18])
			    
		    i=i+1
		    input.close()


    def readphotoz(self,file):
        input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

        self.idnew =[]
	self.idold =[]
	self.bestz = N.zeros(ngal,'f')
 	self.membflag = N.zeros(ngal,'f')
	self.lowz = N.zeros(ngal,'f')
        self.highz = N.zeros(ngal,'f')
        self.wmin = N.zeros(ngal,'f')
        self.Pclust = N.zeros(ngal,'f')

        self.LVlow = N.zeros(ngal,'f')
        self.LVbest = N.zeros(ngal,'f')
        self.LVhigh = N.zeros(ngal,'f')
        self.LRlow = N.zeros(ngal,'f')
        self.LRbest = N.zeros(ngal,'f')
        self.LRhigh = N.zeros(ngal,'f')
        self.UBlow = N.zeros(ngal,'f')
        self.UBbest = N.zeros(ngal,'f')
        self.UBhigh = N.zeros(ngal,'f')
        self.BVlow = N.zeros(ngal,'f')
        self.BVbest = N.zeros(ngal,'f')
        self.BVhigh = N.zeros(ngal,'f')
        self.UVlow = N.zeros(ngal,'f')
        self.UVbest = N.zeros(ngal,'f')
        self.UVhigh = N.zeros(ngal,'f')



        self.LUlowzclust = N.zeros(ngal,'f')
        self.LUbestzclust = N.zeros(ngal,'f')
        self.LUhighzclust = N.zeros(ngal,'f')
        self.LBlowzclust = N.zeros(ngal,'f')
        self.LBbestzclust = N.zeros(ngal,'f')
        self.LBhighzclust = N.zeros(ngal,'f')
        self.LVlowzclust = N.zeros(ngal,'f')
        self.LVbestzclust = N.zeros(ngal,'f')
        self.LVhighzclust = N.zeros(ngal,'f')
        self.LRlowzclust = N.zeros(ngal,'f')
        self.LRbestzclust = N.zeros(ngal,'f')
        self.LRhighzclust = N.zeros(ngal,'f')
        self.LIlowzclust = N.zeros(ngal,'f')
        self.LIbestzclust = N.zeros(ngal,'f')
        self.LIhighzclust = N.zeros(ngal,'f')
        self.UBlowzclust = N.zeros(ngal,'f')
        self.UBbestzclust = N.zeros(ngal,'f')
        self.UBhighzclust = N.zeros(ngal,'f')
        self.BVlowzclust = N.zeros(ngal,'f')
        self.BVbestzclust = N.zeros(ngal,'f')
        self.BVhighzclust = N.zeros(ngal,'f')
        self.UVlowzclust = N.zeros(ngal,'f')
        self.UVbestzclust = N.zeros(ngal,'f')
        self.UVhighzclust = N.zeros(ngal,'f')
        self.Membflag = N.zeros(ngal,'f')
        self.LJlowzclust = N.zeros(ngal,'f')
        self.LJbestzclust = N.zeros(ngal,'f')
        self.LJhighzclust = N.zeros(ngal,'f')
        self.LKlowzclust = N.zeros(ngal,'f')
        self.LKbestzclust = N.zeros(ngal,'f')
        self.LKhighzclust = N.zeros(ngal,'f')
        self.xcorr = N.zeros(ngal,'f')
        self.ycorr = N.zeros(ngal,'f')
        self.starflag = N.zeros(ngal,'f')

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
	    self.idnew.append(t[0])
	    self.idold.append(t[1])
	    self.bestz[i] = float(t[3])
	    self.membflag[i]= float(t[56])
	 
	    self.lowz[i] = float(t[2])
            self.highz[i] = float(t[4])
            self.wmin[i] = float(t[30])
            self.Pclust[i] = float(t[31])

            self.LVlow[i] = float(t[12])#rest-frame luminosities using best z, good for field
            self.LVbest[i] = float(t[13])
            self.LVhigh[i] = float(t[14])
            self.LRlow[i] = float(t[15])
            self.LRbest[i] = float(t[16])
            self.LRhigh[i] = float(t[17])

            self.UBlow[i] = float(t[21])
            self.UBbest[i] = float(t[22])
            self.UBhigh[i] = float(t[23])
            self.BVlow[i] = float(t[24])
            self.BVbest[i] = float(t[25])
            self.BVhigh[i] = float(t[26])
            self.UVlow[i] = float(t[27])
            self.UVbest[i] = float(t[28])
            self.UVhigh[i] = float(t[29])


            self.LUlowzclust[i] = float(t[32])
            self.LUbestzclust[i] = float(t[33])
            self.LUhighzclust[i] = float(t[34])
            self.LBlowzclust[i] = float(t[35])
            self.LBbestzclust[i] = float(t[36])
            self.LBhighzclust[i] = float(t[37])
            self.LVlowzclust[i] = float(t[38])
            self.LVbestzclust[i] = float(t[39])
            self.LVhighzclust[i] = float(t[40])
            self.LRlowzclust[i] = float(t[41])
            self.LRbestzclust[i] = float(t[42])
            self.LRhighzclust[i] = float(t[43])
            self.LIlowzclust[i] = float(t[44])
            self.LIbestzclust[i] = float(t[45])
            self.LIhighzclust[i] = float(t[46])
            self.UBlowzclust[i] = float(t[47])
            self.UBbestzclust[i] = float(t[48])
            self.UBhighzclust[i] = float(t[49])
            self.BVlowzclust[i] = float(t[50])
            self.BVbestzclust[i] = float(t[51])
            self.BVhighzclust[i] = float(t[52])
            self.UVlowzclust[i] = float(t[53])
            self.UVbestzclust[i] = float(t[54])
            self.UVhighzclust[i] = float(t[55])
            self.LJlowzclust[i] = float(t[57])
            self.LJbestzclust[i] = float(t[58])
            self.LJhighzclust[i] = float(t[59])
            self.LKlowzclust[i] = float(t[60])
            self.LKbestzclust[i] = float(t[61])
            self.LKhighzclust[i] = float(t[62])
            self.xcorr[i] = float(t[102])
            self.ycorr[i] = float(t[103])
            self.starflag[i] = float(t[104])
	    i = i+1
	    
	input.close()


    def match24sextractor(self,file):
        input=open(file,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

    def readgimmorphology(self,file1,file2): #read gim morphology and match to ediscs
        input=open(file1,'r')
	#get number of galaxies
	ngal=0
	for line in input:
	    if line.find('#') > -1:
		continue
	    ngal=ngal+1
	input.close()

	self.hstgimID = [] #id for gim
	self.hstgimstargal = N.zeros(ngal,'f') #star(0) or galaxy (1)
	self.hstgimtype = N.zeros(ngal,'f') #galaxy type

	input=open(file1,'r') #opens hstgim id file
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue   

	    t=line.split()
	    self.hstgimID.append(t[0])

	    i=i+1
	input.close()

	input=open(file2,'r') #opens hstgim morphology data
        i=0
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue   

	    t=line.split()
	    self.hstgimstargal[i]=float(t[5])
	    self.hstgimtype[i]=float(t[12])
	    
	    i=i+1
	input.close()

        self.matchmorphgimediscs = N.zeros(len(self.ediscsID), 'i') #match morphology to ediscs
        self.matchflagmorphgimediscs = N.zeros(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.hstgimID)):
			hstgimname=str(self.hstgimID[j])
			if ediscsname.find(hstgimname) > -1:
				self.matchmorphgimediscs[i]=j
				self.matchflagmorphgimediscs[i]=1
				break



    def readvismorphology(self,file): #read vis morphology and match to ediscs
        input=open(file,'r')
        #get number of galaxies
        ngal=0
        for line in input:
            if line.find('#') > -1:
                continue
            ngal=ngal+1
        input.close()

        self.hstvisID = []
        self.hstvisra = N.zeros(ngal,'f')
        self.hstvisdec = N.zeros(ngal,'f')
        self.hstvisIauto =  N.zeros(ngal,'f')
        self.hstvisnumtype = N.zeros(ngal,'f')
        self.hstvisflagfeatureless = N.zeros(ngal,'f')
        self.hstvisflagedgeon = N.zeros(ngal,'f')
        self.hstvisflagtoosmall = N.zeros(ngal,'f')
        self.hstvisflagLSB = N.zeros(ngal,'f')
        self.hstvisflagcosray = N.zeros(ngal,'f')
        self.hstvisflagdust = N.zeros(ngal,'f')
        self.hstvisflagdistorted = N.zeros(ngal,'f')


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
            self.hstvisID.append(t[0])
            self.hstvisra[i] = float(t[1])
            self.hstvisdec[i] = float(t[2])
            self.hstvisIauto[i] = float(t[3])
            self.hstvisnumtype[i] = float(t[4])
            self.hstvisflagfeatureless[i] = float(t[5])
            self.hstvisflagedgeon[i] = float(t[6])
            self.hstvisflagtoosmall[i] = float(t[7])
            self.hstvisflagLSB[i] = float(t[8])
            self.hstvisflagcosray[i] = float(t[9])
            self.hstvisflagdust[i] = float(t[10])
            #self.hstvisdistorted[i] = float(t[11])

            i = i+1
        input.close()
	
        self.matchmorphvisediscs = N.zeros(len(self.ediscsID), 'i') #match morphology to ediscs
        self.matchflagmorphvisediscs = N.zeros(len(self.ediscsID), 'i')
	i=0
	for i in range(len(self.ediscsID)):
		j=0
		ediscsname=str(self.ediscsID[i])
		    #print "checking ha id \n"
		for j in range(len(self.hstvisID)):
			hstvisname=str(self.hstvisID[j])
			if ediscsname.find(hstvisname) > -1:
				self.matchmorphvisediscs[i]=j
				self.matchflagmorphvisediscs[i]=1
				break


    def readspectroscopy(self,file):#for reading .v1.spe files
        input=open(file,'r')
        #get number of galaxies
	ngal=0
	for line in input:
		if line.find('#') > -1:
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue
		if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			continue
		ngal=ngal+1
       	input.close()
	
	self.specediscsID = []
	self.specmembership = N.zeros(ngal,'i')
	self.specz = N.zeros(ngal,'f')
	self.spectype = N.zeros(ngal,'f')
	self.specEWOII = N.zeros(ngal,'f')
	self.specEWOIIflag = N.zeros(ngal,'f')

		
		
	input=open(file,'r')
	i=0
	for line in input:
		if line.find('#') > -1: #skip lines with '#' in them
			continue
		if line.find('\\') > -1: #skip lines with '#' in them
			continue
		if line.find('|') > -1: #skip lines with '#' in them
			continue
		if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			continue
		t = line.split()
		name=t[4]
		if name.find(':') > -1:
			self.specediscsID.append(name[:(len(name)-1)])
		else:
			self.specediscsID.append(name)
		try:
			self.specz[i]=float(t[6])
		except:
			self.specz[i]=-99.
		#memb=t[3]
		try:
			self.specmembership[i] = float(t[3])
			#if self.specmembership[i] > 0.1:#test to make sure galaxy is w/in 3 sigma
			#	dv=(my.getvfromz(self.specz[i])-my.getvfromz(self.zcl))/self.sigma
			#	if abs(dv) > 3:
			#		self.specmembership[i]=0.
				
		except:
			self.specmembership[i] = 0.
		try:
			self.spectype[i] = float(t[9])
		except:
			self.spectype[i] = -99.
		try:			
			self.specEWOII[i]=float(line[80:85])
			self.specEWOIIflag[i]=1.
		except:
			self.specEWOII[i]=-99.	
		i = i+1
       	input.close()

   

	self.matchspecediscs = N.zeros(len(self.ediscsID),'i')
	self.matchflagspecediscs = N.zeros(len(self.ediscsID),'i')

	i=0
	for i in range(len(self.idold)):
		j=0
		ediscsoldname=str(self.idold[i])
		    #print "checking ha id \n"
		for j in range(len(self.specediscsID)):
			specname=str(self.specediscsID[j])
			if ediscsoldname.find(specname) > -1:
				try:
					self.matchspecediscs[i]=j
					self.matchflagspecediscs[i]=1
					break
				except IndexError:
					print "Index error in match spec"
					print "len(self.matchspecediscs)",len(self.matchspecediscs),len(self.ediscsID)
					print "i, j",i,j
    def readspectroscopyv2(self,file):
	    #update this for reading v24

	    input=open(file,'r')
         #get number of galaxies
	    ngal=0
	    for line in input:
		    if line.find('#') > -1:
			    continue
		    if line.find('\\') > -1: #skip lines with '#' in them
			    continue
		    if line.find('|') > -1: #skip lines with '#' in them
			    continue
		    if line[0:4].find('cl') < 0: #skip lines with out 'cl1040' in them
			    continue
		    ngal=ngal+1
	    input.close()
       
	    self.specediscsID = []
	    self.specmembership = N.zeros(ngal,'i')
	    self.specz = N.zeros(ngal,'f')
	    self.spectype = N.zeros(ngal,'f')
	    self.specEWOII = N.zeros(ngal,'f')
	    self.specEWOIIflag = N.zeros(ngal,'f')

               
               
	    input=open(file,'r')
	    i=0
	    for line in input:
		    if line.find('#') > -1: #skip lines with '#' in them
			    continue
		    if line.find('\\') > -1: #skip lines with '#' in them
			    continue
		    if line.find('|') > -1: #skip lines with '#' in them
			    continue
		    if line[0:4].find('c') < 0: #skip lines with out 'cl1040' in them
			    continue
		    t = line.split()
		    name=t[3]#old ediscs ID
		    if name.find(':') > -1:
			    #self.specediscsID.append(name[:(len(name)-1)])
			    self.specediscsID.append(name)#keep colon so we can treat multiple matches differently
		    else:
			    self.specediscsID.append(name)
		    try:
			    self.specz[i]=float(t[9])
		    except:
			    n=t[9]
			    if n.find(':') > -1:#redshift has : at end, indicating uncertainty
				    b=n.split(':')
				    self.specz[i]=float(b[0])
			    else:
				    self.specz[i]=-99.
		    memb=t[8]
		    try:
			    memb=int(memb)
			    self.specmembership[i]=memb
		    except:
			    self.specmembership[i]=0.
		    spec=t[10]
		    if spec.find('3') > -1:
			    self.specEWOIIflag[i] = 1.
		    if spec.find('4') > -1:
			    self.specEWOIIflag[i] = 1.
		    if spec.find('2') > -1:#weak emission lines
			    self.specEWOIIflag[i] = 1.
		    try:
			    self.spectype[i] = float(t[10])
		    except:
			    self.spectype[i] = -99.
			    self.specEWOII[i]=-99. 
		    i = i+1
	    input.close()


	    self.matchspecediscs = N.zeros(len(self.ediscsID),'i')
	    self.matchflagspecediscs = N.zeros(len(self.ediscsID),'i')

	    i=0
	    for i in range(len(self.idold)):
                j=0
                ediscsoldname=str(self.idold[i])
                    #print "checking ha id \n"
                for j in range(len(self.specediscsID)):
                        specname=str(self.specediscsID[j])
                        if specname.find(ediscsoldname) > -1:
				if specname.find(':') > -1:#need to deal w/multiple matches
					self.matchspecediscs[i]=j
					self.matchflagspecediscs[i]=1
					#look for other matches
					matchindex=[j]
					membflag=[specmembership[j]]
					for k in range(j,j+5):#assumes a max of 5 multiple matches
						t=self.specediscsID[k].split(':')
						name2=t[0]
						if specname.find(name2) > -1:
							matchindex.append(k)
							membflag.append(specmembership[k])
					#keep the member
					for k in (matchindex):
						k=int(k)
						try:
							if float(membflag[k]) > 0.1:
								self.matchspecediscs[i]=k
								self.matchflagspecediscs[i]=1
								break
						except:
							print 'Problem matching crazy spec catalogs, just thought you should know'

							
				else:
					self.matchspecediscs[i]=j
					self.matchflagspecediscs[i]=1
					break
	    print "Number of matches with spectroscopy = ",N.sum(self.matchflagspecediscs)


    def memb(self):
	    nspec=0.
	    nspecphotagree=0
	    nspecphotdisagree=0

	    self.defmembflag=N.zeros(self.nediscs,'i')
	    self.newspecmatchflag=-1.*N.ones(self.nediscs,'i')
	    for i in range(len(self.ediscsID)):
		    if (self.membflag[i] > 0) & (self.matchflagspecediscs[i] < 1.):
			    self.defmembflag[i] = 1
			    continue
		    if (self.matchflagspecediscs[i] > 0.):
			    #print "got here",self.matchflag24[i]
			    if (self.matchflag24[i] > 0.):
				    nspec=nspec + 1
			    j=self.matchspecediscs[i]
			    diff=abs(self.specmembership[j]-self.membflag[i])
			    if (self.matchflag24[i] > 0.):
				    if (diff < .1):
					    nspecphotagree=nspecphotagree + 1.
				    else:
					    nspecphotdisagree=nspecphotdisagree + 1.
			    if self.specmembership[j] > 0.:
				    self.defmembflag[i]=1.
			            continue
		    #print i, "defmembflag = ",self.defmembflag[i]
		    

	    for i in range(len(self.newspecmatchflag)):
		    if (self.matchflagspecediscs[i] > 0.5):
			    if (self.specmembership[self.matchspecediscs[i]] > 0.5):
				    self.newspecmatchflag[i] = 1.
			    elif (self.specmembership[self.matchspecediscs[i]] < 0.5):
				    self.newspecmatchflag[i] = 0.
			    j=self.matchspecediscs[i]
			    #print i,j,self.specmembership[j],self.newspecmatchflag[i]

	    print "24 micron sources only:  Nspec, nphotspecagree, %, nphotspecdisagree, %"
	    try:
		    print nspec, nspecphotagree, nspecphotagree/(1.*nspec), nspecphotdisagree,nspecphotdisagree/(1.*nspec)
	    except ZeroDivisionError:
		    print "Warning:  nspec = 0!!!"
	    nspecmemb=len(N.compress((self.matchflag24>0.) & (self.newspecmatchflag > 0.),self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress((self.matchflag24>0.) & (self.newspecmatchflag > 0.) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b,c)=my.ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of spec members & 24micron that are photoz members = %5.3f - %5.3f +%5.3f"%(a,b,c)
	    nspecmemb=len(N.compress((self.matchflag24>0.) & (abs(self.newspecmatchflag) < .1) ,self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress((self.matchflag24>0.) & (abs(self.newspecmatchflag) < .1) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b,c)=my.ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of spec non-members & 24 micron that are photoz members = %5.3f - %5.3f +%5.3f"%(a,b,c)," Nspecmemb = ",nspecmemb

	    nspecmemb=len(N.compress( (abs(self.newspecmatchflag) < .1) ,self.newspecmatchflag))
	    nspecmembNphotmemb=len(N.compress( (abs(self.newspecmatchflag) < .1) & (self.membflag > 0.),self.newspecmatchflag))
	    (a,b,c)=my.ratioerror(1.*nspecmembNphotmemb,1.*nspecmemb)
	    print "Fraction of ALL spec non-members that are photoz members = %5.3f - %5.3f +%5.3f"%(a,b,c)," Nspec nonmemb = ",nspecmemb
    def readflat24(self,file):
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
                   
                   
            self.idflat24 = N.zeros(ngal,'f')
            self.imagexflat24 = N.zeros(ngal,'f')
            self.imageyflat24  = N.zeros(ngal,'f')
            self.raflat24 = N.zeros(ngal,'f')
            self.decflat24 = N.zeros(ngal,'f')
            self.fflat24 = N.zeros(ngal,'f') #flux
            self.errfflat24 = N.zeros(ngal,'f')
                   
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
                    (self.idflat24[i],self.imagexflat24[i],self.imageyflat24[i],self.raflat24[i],self.decflat24[i],self.fflat24[i],self.errfflat24[i])=(float(t[0]),float(t[1]),float(t[2]),float(t[3]),float(t[4]),float(t[5]),float(t[6]))
                    i=i+1
            input.close()
            outfile=open('xyflat24.dat','w')
            for i in range(len(self.imagex24)):
            #print self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i]
                    string="%6.2f %6.2f %8.1f %6.2f \n" % (self.imagex24[i],self.imagey24[i],self.f24[i],self.errf24[i])
                    outfile.write(string)
            outfile.close()
                                   
            #delta=3.#max allowed offset (in arcseconds) between matched sources
            #delta=delta/3600.
            x1=self.ra
            y1=self.dec
            x2=self.raflat24
            y2=self.decflat24
                                   
            self.matchediscsflat24 = N.zeros(len(x1), 'i')
            self.nmatchediscsflat24 = N.zeros(len(x1), 'i')
            self.matchflagediscsflat24 = N.ones(len(x1), 'i')
            for i in range(len(x1)):
                    (self.matchediscsflat24[i],self.matchflagediscsflat24[i],self.nmatchediscsflat24[i]) = findnearest(x1[i],y1[i],x2,y2,delta)

    def mastertable(self,name):
	    output=open('mastertable.fmt','w')
	    output.write('#1 - ediscsID \n')
	    output.write('#2 - ediscs RA \n')
	    output.write('#3 - ediscs Dec \n')
	    output.write('#4 - Halpha EW (A) \n')
	    output.write('#5 - err Halpha EW (A) \n')
	    output.write('#6 - Halpha SFR (Msol/yr, 1 mag extinct) \n')
	    output.write('#7 - Halpha err SFR (Msol/yr, 1 mag extinct) \n')
	    output.write('#8 - match flag Halpha: -1, 0, 1  \n')
	    output.write('#9 - GIM2D Type of HST image \n')
	    output.write('#10 - matchflag GIM2D Type of HST image \n')
	    output.write('#11 - HST visual morph \n')
	    output.write('#12 - HST visual morph matchflag\n')
	    output.write('#13 - 24 um total flux (uJy) \n')
	    output.write('#14 - error in 24 um total flux (uJy) \n')
	    output.write('#15 - matchflag 24 um \n')
	    output.write('#16-25 - isophotal mag and err: V,R,I,J,K \n')
	    output.write('#26-35 - 1" aper mag and err: V,R,I,J,K \n')
	    output.write('#36 - photoz memb flag \n')
	    output.write('#37 - spec membership flag: -1, 0, 1\n')
	    output.write('#38 - photoz + spec membership flag: -0, 1\n')
	    output.write('#39 - spec redshift\n')
	    output.write('#40 - spec type:1=abs-line; 2=abs-line+weak emission; 3=emission w/EW(OII)<25; 4=strong emission w/EW(OII)\n')
	    output.write('#41 - EW OII (A)\n')
	    output.write('#42 - match flag spectroscopy\n')
	    output.write('#43 - best photoz \n')
	    output.write('#44 - low photoz \n')
	    output.write('#45 - high photoz \n')
	    output.write('#46 - wmin: fraction of total exposure in near-IR\n')
	    output.write('#47 - pclust: probability of member based on photozs \n')
	    output.write('#48-62 - low, best, high rest-frame luminosity in U,B,V,R,I \n')
	    output.write('#63-71 - low, best, high rest-frame colors  U-B,B-V,U-V \n')
	    output.write('#72-77 - low, best, high rest-frame luminosity in J,K \n')
	    output.write('#78 - x position on geo corrected vlt image \n')
	    output.write('#79 - y position on geo corrected vlt image \n')
	    output.write('#80 - starflag \n')
	    output.write('#81 - K-band flux \n')
	    output.write('#82 - Star-forming flag based on Halpha criteria \n')
	    output.write('#83 - match flag for ediscs spectroscopy \n')
	    output.write('#84 - EW OII flag for ediscs spectroscopy \n')
	    output.write('#85 - number of optical ediscs galaxies within 2arcsec of 24um source \n')
	    output.write('#86 - matchflagediscsirac\n') 
	    output.write('#87-90 - iracf1 iracf2 iracf3 iracf4 (uJy) \n')
	    output.write('#91-94 - erriracf1 erriracf2 erriracf3 erriracf4 (uJy) \n') 
	    output.write('#95 - iracsexflag0 \n') 
	    output.write('#96 - iracsexflag1 \n') 
	    output.write('#97-100 - iracwch1 iracwch2 iracwch3 iracwch4 \n') 
	    output.write('#101 - iracwmin \n') 
	    output.write('#102 - nmatchediscsirac: number of optical sources w/in 2arcsec radius of irac source\n') 
	    output.write('#103 - 24 micron luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#104 - error in 24 micron luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#105 - Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#106 - error in Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    output.write('#107 - low, best, high rest-frame luminosity in U,B,V,R,I \n')
	    output.write('#63-71 - low, best, high rest-frame colors  U-B,B-V,U-V \n')
	    output.write('#72 - RA on 24um image\n')
	    output.write('#73 - Dec on 24um image\n')
	    output.write('#73 - flux80flag; is flux above 80% completeness limit? 1=yes, 0=no\n')


	    #output.write('#106 - error in Halpha luminosity (H0=70 km/s/Mpc)\n') 
	    output.close()
	    outfile=str(name)+'mastertable.dat'
	    output= open(outfile,'w')
	    outfile2=str(name)+'mastertable24.dat'
	    output2= open(outfile2,'w')
	    outfile3=str(name)+'AllMembers.dat'
	    output3=open(outfile3,'w')
	    outfile4=str(name)+'Members24.dat'
	    output4=open(outfile4,'w')
	    outfile5=str(name)+'ediscs24.v2.dat'#tables for ediscs collaboration
	    output5=open(outfile5,'w')
            output5.write('#ediscsID     matchflag24  f24(uJy) errf24(uJy) Vmag(1arc) errVmag photmembflag specmembflag photo+specmemb wmin spectype OIIEW HaFlag HaEW HaEWerr HaSFR HaSFRerr\n')

            output5.write('#matchflag24: -1=not on 24um image; 0=on 24um image but no 24um detection; 1=24um source \n')
            output5.write('#NOTE: f24 is meaningful only if matchflag24 > 0 \n')
            output5.write('#photmembflag: 0 = not a member; 1 = a member according to photoz \n')
	    output5.write('#specmembflag: -1 = no spec; 0 = spec but not a member; 1 = spec and a member \n')
            output5.write('#phot+specmemb: 0 = if photoz or spec says not a member; 1 = if spec or photoz says is a member \n')
            output5.write('#wmin = fraction of full near-IR exposure.  photozs are unreliable for wmin < 0.3\n')
	    output.write('#ediscsID  ra  dec  xcorr  ycorr  starflag  EWha  EWhaerr  SFR  SFRerr  matchflaghaediscs SF_flag hstgimtype  matchflagmorphgimediscs  hstvisnumtype  matchflagmorphvisediscs  matchflag24  f24  errf24  nmatchediscs24  misoV  misoeVapsim  misoR  misoeRapsim  misoI  misoeIapsim  misoJ  misoeJapsim  misoK  misoeKapsim   magV_1  mageVapsim  magR_1  mageRapsim  magI_1  mageIapsim  magJ_1  mageJapsim  magK_1  mageKapsim  membflag  newspecmatchflag defmembflag  specz  spectype  specEWOII matchflagspecediscs specEWOIIflag  bestz   lowz   highz   wmin   Pclust  LUlowzclust  LUbestzclust  LUhighzclust  LBlowzclust  LBbestzclust  LBhighzclust  LVlowzclust  LVbestzclust  LVhighzclust  LRlowzclust  LRbestzclust  LRhighzclust  LIlowzclust  LIbestzclust  LIhighzclust  LJlowzclust  LJbestzclust  LJhighzclust  LKlowzclust  LKbestzclust  LKhighzclust  fluxK  UBlowzclust  UBbestzclust  UBhighzclust  BVlowzclust  BVbestzclust  BVhighzclust  UVlowzclust  UVbestzclust  UVhighzclust  matchflagediscsirac iracf1 iracf2  iracf3 iracf4 erriracf1 erriracf2 erriracf3 erriracf4 iracsexflag0 iracsexflag1 iracwch1 iracwch2 iracwch3 iracwch4 iracwmin nmatchediscsirac L24 errL24 LHa errLHa  LVlow  LVbest  LVhigh  LRlow  LRbest  LRhigh UBlow  UBbest  UBhigh  BVlow  BVbest  BVhigh  UVlow  UVbest  UVhigh RA24 Dec24 flux80flag \n')
	    for i in range(len(self.ediscsID)):
		    try:
			    ii=self.matchediscsirac[i]
	                    ih=self.matchhaediscs[i]
			    i24=self.matchediscs24[i]
			    ispec=self.matchspecediscs[i]
#			    string="%s %12.6f %12.6f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %8.4e %8.4e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %12.6f %12.6f %5.2f %5.2f %5.2f %7.2f %7.2f %7.2f %7.2f %8.4e %7.2f %7.2f %7.2f %5.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n" % (self.ediscsID[i],self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EWha[ih],self.EWhaerr[ih],self.SFR[ih],self.SFRerr[ih],self.matchflaghaediscs[i],self.SF_flag[ih], self.hstgimtype[self.matchmorphgimediscs[i]],self.matchflagmorphgimediscs[i],self.hstvisnumtype[self.matchmorphvisediscs[i]],self.matchflagmorphvisediscs[i],self.matchflag24[i],self.f24[i24],self.errf24[i24],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i], self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i], self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.specz[ispec],self.spectype[ispec],self.specEWOII[ispec],self.matchflagspecediscs[i],self.specEWOIIflag[ispec],self.bestz[i], self.lowz[i], self.highz[i], self.wmin[i], self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[ii],self.iracf2[ii],self.iracf3[ii],self.iracf4[ii],self.erriracf1[ii],self.erriracf2[ii],self.erriracf3[ii],self.erriracf4[ii],self.iracsexflag0[ii],self.iracsexflag1[ii],self.iracwch1[ii],self.iracwch2[ii],self.iracwch3[ii],self.iracwch4[ii],self.iracwmin[ii],self.nmatchediscsirac[i],self.l24[i24],self.errl24[i24],self.lha[ih],self.lhaerr[ih],self.snr24[i24],self.imagex24[i24],self.imagey24[i24],self.fap1[i24],self.fap2[i24],self.fap3[i24],self.fap4[i24],self.fap5[i24],self.fap6[i24],self.fap7[i24],self.fap8[i24],self.fap9[i24],self.fap10[i24],self.errfap1[i24],self.errfap2[i24],self.errfap3[i24],self.errfap4[i24],self.errfap5[i24],self.errfap6[i24],self.errfap7[i24],self.errfap8[i24],self.errfap9[i24],self.errfap10[i24],self.LVlow[i],  self.LVbest[i],  self.LVhigh[i],  self.LRlow[i],  self.LRbest[i],  self.LRhigh[i], self.UBlow[i],  self.UBbest[i],  self.UBhigh[i],  self.BVlow[i],  self.BVbest[i],  self.BVhigh[i],  self.UVlow[i],  self.UVbest[i],  self.UVhigh[i],self.ra24[i24],self.dec24[i24],self.flux80flag[i],self.photmembflag[i],self.supermembflag[i],self.MR[i],self.MU[i],self.MV[i],self.MB[i],self.stellmass[i],self.MRbestz[i],self.MVbestz[i],self.MUbestz[i],self.redflag[i],self.L24[i],self.errL24[i],self.Lir[i],self.errLir[i],self.Lirbestz[i],self.errLirbestz[i],self.SFRir[i],self.SFRirerr[i])
			    string="%s %12.6f %12.6f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %8.4e %8.4e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.2f %7.2f %7.2f %7.2f %7.2f %7.4f %7.3f %7.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %12.6f %12.6f %5.2f \n" % (self.ediscsID[i],self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EWha[ih],self.EWhaerr[ih],self.SFR[ih],self.SFRerr[ih],self.matchflaghaediscs[i],self.SF_flag[ih], self.hstgimtype[self.matchmorphgimediscs[i]],self.matchflagmorphgimediscs[i],self.hstvisnumtype[self.matchmorphvisediscs[i]],self.matchflagmorphvisediscs[i],self.matchflag24[i],self.f24[i24],self.errf24[i24],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i], self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i], self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.specz[ispec],self.spectype[ispec],self.specEWOII[ispec],self.matchflagspecediscs[i],self.specEWOIIflag[ispec],self.bestz[i], self.lowz[i], self.highz[i], self.wmin[i], self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[ii],self.iracf2[ii],self.iracf3[ii],self.iracf4[ii],self.erriracf1[ii],self.erriracf2[ii],self.erriracf3[ii],self.erriracf4[ii],self.iracsexflag0[ii],self.iracsexflag1[ii],self.iracwch1[ii],self.iracwch2[ii],self.iracwch3[ii],self.iracwch4[ii],self.iracwmin[ii],self.nmatchediscsirac[i],self.l24[i24],self.errl24[i24],self.lha[ih],self.lhaerr[ih],self.snr24[i24],self.imagex24[i24],self.imagey24[i24],self.fap1[i24],self.fap2[i24],self.fap3[i24],self.fap4[i24],self.fap5[i24],self.fap6[i24],self.fap7[i24],self.fap8[i24],self.fap9[i24],self.fap10[i24],self.errfap1[i24],self.errfap2[i24],self.errfap3[i24],self.errfap4[i24],self.errfap5[i24],self.errfap6[i24],self.errfap7[i24],self.errfap8[i24],self.errfap9[i24],self.errfap10[i24],self.LVlow[i],  self.LVbest[i],  self.LVhigh[i],  self.LRlow[i],  self.LRbest[i],  self.LRhigh[i], self.UBlow[i],  self.UBbest[i],  self.UBhigh[i],  self.BVlow[i],  self.BVbest[i],  self.BVhigh[i],  self.UVlow[i],  self.UVbest[i],  self.UVhigh[i],self.ra24[i24],self.dec24[i24],self.flux80flag[i])
			    
			    output.write(string)
			    #print i,self.matchflaghaediscs[i],self.SFR[ih],self.SFRerr[ih],self.lha[ih],self.lhaerr[ih]
		    except IndexError:
			    print  "Index Error in mastertable",i,len(self.ediscsID)

		    #print "%s \n" %(self.ediscsID[i])
		    #print string

		    if self.matchflag24[i] > 0:
			    string5="%s %7.2f %8.2f %8.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n"%(self.ediscsID[i],self.matchflag24[i],self.f24[self.matchediscs24[i]],self.errf24[self.matchediscs24[i]],self.magV[i],self.mageVapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.wmin[i],self.spectype[self.matchspecediscs[i]],self.specEWOII[self.matchspecediscs[i]],self.matchflaghaediscs[i],self.EWha[self.matchhaediscs[i]],self.EWhaerr[self.matchhaediscs[i]],self.SFR[self.matchhaediscs[i]],self.SFRerr[self.matchhaediscs[i]])
		    else:
			    string5="%s %7.2f     0.00     0.00 %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f    0.00    0.00    0.00    0.00    0.00    0.00    0.00\n"%(self.ediscsID[i],self.matchflag24[i],self.magV[i],self.mageVapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.wmin[i])
		    output5.write(string5)
		    if self.matchflag24[i] > 0:
			    output2.write(string)

		    if self.defmembflag[i] > 0:
			    if self.matchflagmorphvisediscs[i] > 0.:
				    string2="%s %7.2f \n" % (self.ediscsID[i],self.hstvisnumtype[self.matchmorphvisediscs[i]])
				    output3.write(string2)
				    if self.matchflag24[i] > 0:
					    string2="%s %7.2f %7.2f %7.2f \n" % (self.ediscsID[i],self.hstvisnumtype[self.matchmorphvisediscs[i]],self.f24[self.matchediscs24[i]],self.errf24[self.matchediscs24[i]])
					    output4.write(string2)

		    i=i+1
	    output.close()
	    output2.close()
	    output3.close()
	    output4.close()
	    output5.close()

    def count(self):
	total = 0.
	photo = 0.
	total24 = 0.
	photo24 = 0.
	den = 0.
	num = 0.
	den24 = 0.
	num24 = 0.
	N = 0.
	N24 = 0.
	Nha = 0.
	totalha = 0.
	b = 0.
	bspec = 0.
	spec24 = 0.
	photoha = 0.
	photob = 0.
	
	self.complete = 0.
	self.contam = 0.
	N = 0.
	N24 = 0.
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .4:
		if self.newspecmatchflag[i] > 0:
		    total = total + 1.
		    if self.membflag[i] > 0:
			photo = photo + 1.
		if self.newspecmatchflag[i] > -1:
		    if self.membflag[i] > 0:
			den = den + 1.
			if self.newspecmatchflag[i] == 0:
			    num = num + 1.
		 
	complete = my.ratioerror(photo,total)
	contam = my.ratioerror(num,den)
	print complete, contam
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .4:
		if self.matchflag24[i] > 0:
		    if self.newspecmatchflag[i] > 0:
			total24 = total24 + 1.
			if self.membflag[i] > 0:
			    photo24 = photo24 + 1.
		    if self.newspecmatchflag[i] > -1:
			if self.membflag[i] > 0:
			    den24 = den24 + 1.
			    if self.newspecmatchflag[i] == 0:
				num24 = num24 + 1.
    	complete24 = my.ratioerror(photo24,total24)
	self.complete = complete24[0]
	contam24 = my.ratioerror(num24,den24)
	self.contam = contam[0]
	print  complete24, contam24
	for i in range(len(self.membflag)):
	    if self.wmin[i] > .4:
		if self.membflag[i] > 0:
		    N = N + 1.		    
		if self.matchflag24[i] > 0:
		    if self.membflag[i] > 0:
			N24 = N24 + 1.
		    if self.newspecmatchflag[i] > 0:
			spec24 = spec24 + 1
		if self.matchflagha[i] > 0:
		    if self.membflag[i] > 0:
			Nha = Nha + 1.
		    if self.newspecmatchflag[i] > 0:
			totalha = totalha + 1.
			if self.membflag[i] > 0:
			    photoha = photoha + 1
		    if self.matchflag24[i] > 0:
			if self.membflag[i] > 0:
			    b = b + 1.
			if self.newspecmatchflag[i] > 0:
			    bspec = bspec + 1.
			    if self.membflag[i] > 0:
				photob = photob + 1.
			    
			
	print total, N, photo
	print spec24, N24, photo24
	print totalha, Nha, photoha
	print bspec, b, photob

    def L24LIsub(self,mark):
	l24=[]
	li=[]
	LIbestzclust=self.LIbestzclust*1.e10/h**2#convert to solar luminosities & cosmol

	for i in range(len(self.L24)):
	    if (self.matchflag24[i] >0.) & (self.photmembflag[i] > 0.):
		l24.append(self.L24[i])
		li.append(LIbestzclust[i])
	l24=N.array(l24,'d')
	li=N.array(li,'d')
	pylab.plot(li,l24,mark,label=self.prefix,markersize=6.)

    def LIRLIsub(self,mark):
	l24=[]
	li=[]
	LIbestzclust=self.LIbestzclust*1.e10/h**2#convert to solar luminosities & cosmol
	Lir=self.L24/.09#approx conversion from Dale models

	for i in range(len(self.L24)):
	    if (self.matchflag24[i] >0.) & (self.photmembflag[i] > 0.):
		    #print "got one"
		    l24.append(Lir[i])
		    li.append(LIbestzclust[i])
	l24=N.array(l24,'d')
	li=N.array(li,'d')
	print self.prefix,"len(l24) = ",len(l24)
	s='CL'+str(self.idra)+'-'+str(self.iddec)
	pylab.plot(li,l24,mark,label=s,markersize=6.)
    def LIRhistsub(self,mark):
	lir=[]
	morph=[]
	morphbell=[]
	specflagbell=[]
	Ibell=[]
	VIbell=[]
	ngal=0
	specz=[]
	photz=[]
	dr=[]
	#Lir=self.L24/.09#approx conversion from Dale models
	for i in range(len(self.L24)):
	    flag=0
	    if self.matchflag24[i] > 0.1:
		    if self.supermembflag[i] > 0.1:
			    if (self.matchflagvistype[i] > 0.1) & (self.vistype[i] < 100):#if it has HST morphology
			    #if self.newspecmatchflag[i] > 0:
				    ngal += 1
				    for j in range(len(Lid)):
					    s=self.ediscsID[i]
					    if s.find(Lid[j]) > -1:
						    m=Lmorph[j]
						    for k in range(len(m)):
							    morphbell.append(m[k])
							    lir.append(self.Lir[i])
							    morph.append(self.vistype[i])
							    Ibell.append(self.magI[i])
							    VIbell.append((self.magV[i]-self.magI[i]))
							    dr.append(self.dr[i])
							    flag=1
							    photz.append(self.zcl)
							    if (self.newspecmatchflag[i] > 0.1):
								    specflagbell.append(1.)
								    specz.append(self.specz[i])
							    else:
								    specflagbell.append(0.)
								    specz.append(-99.)

				    if (flag < 1):
					    #print "No match for %s: HSTvistype = %5.1f, matchf24flag = %3.1f, photmembflag = %3.1f, matchflagvistype = %3.1f, newspecmatchflag = %3.1f"%(self.ediscsID[i],self.vistype[i],(self.matchflag24[i]),(self.photmembflag[i]),(self.matchflagvistype[i]),(self.newspecmatchflag[i]))
					    s="%s \n"%(self.ediscsID[i])
					    morphout.write(s)
	self.lirbell=N.array(lir,'d')
	self.vistypebell=N.array(morph,'f')
	self.morphbell=N.array(morphbell,'f')
	self.Ibell=N.array(Ibell,'f')
	self.VIbell=N.array(VIbell,'f')
	self.specflagbell=N.array(specflagbell,'f')
	self.speczbell=N.array(specz,'f')
	self.photzbell=N.array(photz,'f')
	self.drbell=N.array(dr,'f')
	print "length of arrays", len(self.morphbell),len(self.specflagbell),ngal
	return lir,morph

    def LIRhistsubspec(self,mark):
	lir=[]
	morph=[]
	Lir=self.L24/.09#approx conversion from Dale models
	for i in range(len(self.L24)):
	    if (int(self.matchflag24[i]) & int(self.matchflagvistype[i]) > 0) & (self.newspecmatchflag[i] > 0.) :
		lir.append(Lir[i])
		morph.append(self.vistype[i])
	return lir,morph

    def plotmorph(self):
	prefix=self.prefix
	mall=[]
	m24=[]
	for i in range(len(self.L24)):
	    if (int(self.photmembflag[i]) & int(self.matchflagvistype[i]) > 0.):
		mall.append(self.vistype[i])
		if self.matchflag24[i] > 0.:
		    m24.append(self.vistype[i])
	mall=N.array(mall,'f')
	m24=N.array(m24,'f')

	#xbins=N.array([-7,-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11,66,111],'f')
	xbins=N.array([-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11],'f')
	nmall=N.zeros(len(xbins),'f')
	nm24=N.zeros(len(xbins),'f')
	i=0
	for x in xbins:
	    for m in mall:
		if abs(m-x) < .1:
		    nmall[i]=nmall[i]+1.
	    for m in m24:
		if abs(m-x) < .1:
		    nm24[i]=nm24[i]+1.
	    i=i+1
	Nbin = len(xbins)
	ntot=int(N.sum(nmall))
	ntot24=int(N.sum(nm24))

	nmallerr=N.sqrt(nmall)/N.sum(nmall)*100.
	nmall=nmall/N.sum(nmall)*100.
	nm24err=N.sqrt(nm24)/N.sum(nm24)*100.
	nm24=nm24/N.sum(nm24)*100.
    #ind = N.arange(0,2*Nbin,2)  # the x locations for the groups
	ind = N.arange(Nbin)  # the x locations for the groups
	width = 0.35       # the width of the bars
	p1 = pylab.bar(ind, nmall, width, color='r', yerr=nmallerr)
	
	p2 = pylab.bar(ind+width, nm24, width, color='b', yerr=nm24err)
	
	pylab.text(.5,34.,prefix,fontsize=12.)
	pylab.xlim(-width,len(ind))
	pylab.ylim(0.,40.)
	ytick=N.arange(0.,40.,10.)
	#pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr','?','N/A') ,fontsize=10)
	#xbins=N.array([-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11],'f')
	pylab.xticks(ind+width, ('C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr') ,fontsize=12)
	
	pylab.legend( (p1[0], p2[0]), ('All ('+str(ntot)+' gal)', '24um ('+str(ntot24)+' gal)'), shadow=True)


    def plotpositions(self):
	prefix=self.prefix
	ra=[]
	dec=[]
	ra24=[]
	dec24=[]
	rao2=[]
	deco2=[]

	pra=[]
	pdec=[]
	pra24=[]
	pdec24=[]
	DA=my.DA(self.zcl,h)#kpc/arcsec
	for i in range(len(self.L24)):
	    if (int(self.supermembflag[i]) > 0.):
		    if (self.matchflag24[i] > -.5):
		    #dra=(self.ra[i]-self.rac)*N.cos(self.dec[i]*N.pi/180.)*60./self.r200arcmin
			    dra=(self.ra[i]-self.rac)*60.*N.cos(self.decc/180.*N.pi)#/self.r200arcmin
			    ddec=(self.dec[i]-self.decc)*60.#/self.r200arcmin
			    pra.append(dra)
			    pdec.append(ddec)
			    if (self.matchflag24[i] > 0.):
				    pra24.append(dra)
				    pdec24.append(ddec)

			    if (int(self.newspecmatchflag[i]) > 0):
				    ra.append(dra)
				    dec.append(ddec)

				    if (self.matchflag24[i] > 0.):
					    ra24.append(dra)
					    dec24.append(ddec)
				    if self.specEWOIIflag[i] > 0.1:#use OII instead of Halpha b/c not sure of overlap b/w Halpha and 24 micron image
					    rao2.append(dra)
					    deco2.append(ddec)


	#pylab.plot(pra,pdec,'wo',markeredgecolor='k',markersize=3)
	pylab.plot(pra24,pdec24,'wo',markeredgecolor='r',markersize=4)

	pylab.plot(ra,dec,'ko',markersize=6)
	pylab.plot(ra24,dec24,'ro',markersize=10)
	pylab.plot(rao2,deco2,'bo',markersize=6)
	r2=self.r200arcmin
	x=r2*N.arange(-1.,1+.01,.01)
	y=N.sqrt(r2**2-x**2)
	pylab.plot(x,y,color='0.4',ls='-')
	pylab.plot(x,-1.*y,color='0.4',ls='-')
	nr=4.9
	#pylab.axis('equal')
	ax=pylab.gca()
	name='CL'+str(self.idra)+'-'+str(self.iddec)

	s=name#+', z=%5.3f, $\sigma$=%3i'%(self.z,self.sigma)
	pylab.text(.1,.82,s,fontsize=14,transform=ax.transAxes)
	t=pylab.arange(-4,5,2)
	pylab.xticks(t)
	pylab.yticks(t)

	ax.axis('equal')
	pylab.axis([1.*nr,-1.*nr,-1.*nr,nr],'equal')

	#return ra,dec,ra24,dec24

    def colormag(self):
	    pylab.cla()
	    pylab.clf()
	    name='CL'+str(self.ra)+'-'+str(self.dec)

	    self.colormagallsub()
	    pylab.legend(loc='upper right')

	    pylab.ylabel(r'$\rm V-I$',fontsize=30)
	    pylab.xlabel(r'$\rm I$',fontsize=30)
	    title='V-I vs. I'+str(name)

	    #pylab.axis([19.,24.,0.,5,])
	    file=self.prefix+'colormag.eps'
	    pylab.savefig(file)

    def colormagallsub(self,dm):
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    #dm=2.5*N.log10(self.DL/cl1216.DL)

	    flux1=[]
	    errflux1=[]
	    flux2=[]
	    errflux2=[]	
	    flux1h=[]
	    errflux1h=[]
	    flux2h=[]
	    errflux2h=[]	
	    flux124=[]
	    errflux124=[]
	    flux224=[]
	    errflux224=[]	
	    
	    pflux1=[]
	    perrflux1=[]
	    pflux2=[]
	    perrflux2=[]	
	    pflux1h=[]
	    perrflux1h=[]
	    pflux2h=[]
	    perrflux2h=[]	
	    pflux124=[]
	    perrflux124=[]
	    pflux224=[]
	    perrflux224=[]	
	    n03=0#number w/photoz from zc-0.3 to zc+.3
	    n031=0#number w/photoz from zc-0.3 to z=1
	    n02=0#number w/photoz from zc+/-0.2
	    n03a=0#number w/photoz from zc-0.3  to a max zbest=1
	    s=self.prefix+'RedSeq24.dat'
	    johnfile=open(s,'w')
	    nredseq=0
	    nredseq24=0
	    for i in range(len(self.ra)):
	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newspecmatchflag[i],self.wmin[i]
		    if self.matchflag24[i] > 0.:#count how many 24 micron sources are within broad photoz cuts from cluster
			    if abs(self.zcl-self.bestz[i]) < 0.3:
				    n03=n03+1
				    if self.bestz[i] < 1:
					    n03a=n03a+1
			    if abs(self.zcl-self.bestz[i]) < 0.2:
				    n02=n02+1

			    if (self.bestz[i] > 0.3) & (self.bestz[i] < 1.):
				    n031=n031+1
			    
		    if (self.newspecmatchflag[i] > 0.1) & (self.matchflag24[i] > -1.):#spec memb and on 24um image
			    flux2.append(self.magV[i]-self.magI[i])
			    flux1.append(self.magI[i])
			    if self.specEWOIIflag[i] > 0.1:#use OII instead of Halpha b/c not sure of overlap b/w Halpha and 24 micron image
			    #if self.matchflagha[i] > 0:
				    #if self.SFflag[i] > 0:
				    flux2h.append(self.magV[i]-self.magI[i])
				    flux1h.append(self.magI[i])

			    I=self.magI[i]
			    VI=self.magV[i]-self.magI[i]
			    dmag=abs(-1.*0.09*I + dm - VI)
			    if (dmag < .3):
				    nredseq += 1

			    if self.matchflag24[i] > 0.1:
				    flux224.append(VI)
				    flux124.append(I)
			            #see if galaxy is on red sequence
				    if (dmag < 0.3):
					    s='%s %12.8f %12.8f %8.4f %3.1f \n'%(self.ediscsID[i],self.ra[i],self.dec[i],self.specz[i],self.spectype[i])
					    johnfile.write(s)
					    nredseq24 += 1
			    #print "spec+24 ",self.ediscsID[i],N.log10(self.L24[i]),self.vistype[i]
			    s='circle(%12.8f,%12.8f,4.0" \n'%(self.ra[i],self.dec[i])
			    #dsfile.write(s)
		    elif (self.supermembflag[i] > 0.):
			    pflux2.append(self.magV[i]-self.magI[i])
			    pflux1.append(self.magI[i])
			    if self.matchflagha[i] > 0:
				    if self.SFflag[i] > 0:
					    pflux2h.append(self.magV[i]-self.magI[i])
					    pflux1h.append(self.magI[i])
			    
			    if self.matchflag24[i] > 0:
				    pflux224.append(self.magV[i]-self.magI[i])
				    pflux124.append(self.magI[i])


	    johnfile.close()
	
	    flux1=N.array(flux1,'d')
	    flux2=N.array(flux2,'d')
	    flux1h=N.array(flux1h,'d')
	    flux2h=N.array(flux2h,'d')
	    flux124=N.array(flux124,'d')
	    flux224=N.array(flux224,'d')

	    self.downsize=min(flux124)-min(flux1)

	    pflux1=N.array(pflux1,'d')
	    pflux2=N.array(pflux2,'d')
	    pflux1h=N.array(pflux1h,'d')
	    pflux2h=N.array(pflux2h,'d')
	    pflux124=N.array(pflux124,'d')
	    pflux224=N.array(pflux224,'d')
		
	    #pylab.plot(flux1,flux2,'bo')
	    title='V-I vs. I'+str(name)
	    #pylab.title(name)

	    #pylab.plot(pflux1,pflux2, 'wo',markeredgecolor='k',markeredgewidth=1.,markersize=3,label="Photo-z Members")
	    
	    #pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="Photo-z + 24um")

	    #pylab.plot(pflux1h,pflux2h,'wo',markeredgecolor='b',mew=1.,markersize=12,label="Photo-z + Halpha")
	    
	    #pylab.plot(flux1,flux2, 'ko',markersize=6,label="Spec Members")
	    pylab.plot(flux1,flux2, 'ko',markersize=6,label="_nolegend_")

	    
	    pylab.plot(flux124,flux224,'ro',markersize=10,label=r"$\rm 24\mu m$")
	    #pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label=r"$\rm 24\mu m \ photoz$")
	    pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="_nolegend_")
	    pylab.plot(flux1h,flux2h,'bo',markersize=6,label=r"$\rm [OII]$")
	

	    print self.prefix, "z = ",self.zcl
	    print "Nspec+24, Nphot+24 = ",len(flux124),len(pflux124)
	    print "N 24 with abs(bestz-z) < 0.2 = ",n02
	    print "N 24 with abs(bestz-z) < 0.3 = ",n03
	    print "N 24 with abs(bestz-z) < 0.3, but capping at zbest < 1 = ",n03a
	    print "N 24 with bestz> 0.3 and bestz < 1 = ",n031

	#pylab.legend(loc='lower left')

	    #dm=2.5*N.log10(self.DL/cl1216.DL)
	    x=N.arange(19.8,24.4,1.)
	    y=-.09*(x-20)+dm-1.6
	    pylab.plot(x,y,'k-',label="_nolegend_")
	    y2=y-.3
	    pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    y2=y+.3
	    pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    
	    pylab.xticks(N.arange(19,24))
	    pylab.yticks(N.arange(1,4))
	    pylab.axis([18.5,23.5,0.1,3.9])

	    ax=pylab.gca()
	    s=name#+', z=%5.3f, Nsp=%2d, Nph=%2d'%(self.z,len(flux124),len(pflux124))
	    pylab.text(.2,.8,s,fontsize=14,transform=ax.transAxes)

	    self.nredseq = nredseq
	    self.nredseq24 = nredseq24
	    (a,b,c)=my.ratioerror(nredseq24,nredseq)
	    print 'Fraction of spectroscopic red seq members with 24um emission = %6.3f +/- %6.3f (%i, %i)'%(a,b, nredseq24,nredseq)
#	pylab.axis([0,3.5,14,22])
	#file=str(name)+'colormag.jpg'
	#pylab.savefig(file)



    def colormagallBVsub(self,dm):
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    #dm=2.5*N.log10(self.DL/cl1216.DL)

	    flux1=[]
	    errflux1=[]
	    flux2=[]
	    errflux2=[]	
	    flux1h=[]
	    errflux1h=[]
	    flux2h=[]
	    errflux2h=[]	
	    flux124=[]
	    errflux124=[]
	    flux224=[]
	    errflux224=[]	
	    
	    pflux1=[]
	    perrflux1=[]
	    pflux2=[]
	    perrflux2=[]	
	    pflux1h=[]
	    perrflux1h=[]
	    pflux2h=[]
	    perrflux2h=[]	
	    pflux124=[]
	    perrflux124=[]
	    pflux224=[]
	    perrflux224=[]	
	    n03=0#number w/photoz from zc-0.3 to zc+.3
	    n031=0#number w/photoz from zc-0.3 to z=1
	    n02=0#number w/photoz from zc+/-0.2
	    n03a=0#number w/photoz from zc-0.3  to a max zbest=1
	    nredseq=0
	    nredseq24=0
	    for i in range(len(self.ra)):
	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newspecmatchflag[i],self.wmin[i]
		    if (self.newspecmatchflag[i] > 0.1):#spec memb and on 24um image
			    flux2.append(self.MB[i]-self.MV[i])
			    flux1.append(self.MV[i])
			    if self.specEWOIIflag[i] > 0.1:#use OII instead of Halpha b/c not sure of overlap b/w Halpha and 24 micron image
			    #if self.matchflagha[i] > 0:
				    #if self.SFflag[i] > 0:
				    flux2h.append(self.MB[i]-self.MV[i])
				    flux1h.append(self.MV[i])

			    I=self.magI[i]
			    VI=self.magV[i]-self.magI[i]
			    dmag=abs(-1.*0.09*I + dm - VI)
			    if (dmag < .3):
				    flux224.append(self.MB[i]-self.MV[i])
				    flux124.append(self.MV[i])
			            #see if galaxy is on red sequence
		    elif (self.supermembflag[i] > 0.)  & (self.matchflag24[i] > -1.):
			    pflux2.append(self.MB[i]-self.MV[i])
			    pflux1.append(self.MV[i])
			    
			    I=self.magI[i]
			    VI=self.magV[i]-self.magI[i]
			    dmag=abs(-1.*0.09*I + dm - VI)
			    if (dmag < .3):
				    pflux224.append(self.MB[i]-self.MV[i])
				    pflux124.append(self.MV[i])



	
	    flux1=N.array(flux1,'d')
	    flux2=N.array(flux2,'d')
	    flux1h=N.array(flux1h,'d')
	    flux2h=N.array(flux2h,'d')
	    flux124=N.array(flux124,'d')
	    flux224=N.array(flux224,'d')


	    pflux1=N.array(pflux1,'d')
	    pflux2=N.array(pflux2,'d')
	    pflux1h=N.array(pflux1h,'d')
	    pflux2h=N.array(pflux2h,'d')
	    pflux124=N.array(pflux124,'d')
	    pflux224=N.array(pflux224,'d')
		
	    title='V-I vs. I'+str(name)
	    pylab.plot(flux1,flux2, 'bo',markersize=4,label="Blue Cloud")
	    pylab.plot(flux124,flux224,'ro',markersize=4,label="Red Seq")
	    #pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="_nolegend_")
	    #pylab.plot(flux1h,flux2h,'bo',markersize=6,label=r"$\rm [OII]$")
	
	    x=N.arange(-25.,-15.,1.)
	    y=-.022*(x+20)+.65#-.1*(self.zcl-.48) #cut b/w red and blue galaxies
	    pylab.plot(x,y,'k-',label="_nolegend_")
	    y2=y-.2
	    #pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    y2=y+.2
	    #pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    
	    pylab.xticks(N.arange(-24.,-17.,2))
	    #pylab.yticks(N.arange(1,4))
	    pylab.axis([-24.,-17.,0.1,1.2])

	    ax=pylab.gca()
	    s=name#+', z=%5.3f, Nsp=%2d, Nph=%2d'%(self.z,len(flux124),len(pflux124))
	    pylab.text(.2,.8,s,fontsize=14,transform=ax.transAxes)


    def colormagallUVsub(self,dm):
	    
	    flux1=N.compress((self.supermembflag > 0.1) & (abs(self.matchflag24) < .1), self.MV)
	    flux2=N.compress((self.supermembflag > 0.1) & (abs(self.matchflag24) < .1),self.MU)


	    flux2=flux2-flux1

	    flux124=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1), self.MV)
	    flux224=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1),self.MU)

	    flux224=flux224-flux124

#	    for i in range(len(self.ra)):
#	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newspecmatchflag[i],self.wmin[i]
#		    if ((self.supermembflag[i] > 0.1) &  (self.matchflag24[i] > -.1)) :#spec memb and on 24um image
#			    flux2.append(self.MU[i]-self.MV[i])
#			    flux1.append(self.MV[i])
#
#			    if (self.matchflag24[i] > .1):
#				    flux224.append(self.MU[i]-self.MV[i])
#				    flux124.append(self.MV[i])
#			            #see if galaxy is on red sequence
	

	    pylab.plot(flux1,flux2, 'k.',markersize=1,label="Blue Cloud",markeredgecolor='k',markerfacecolor='None')
	    pylab.plot(flux124,flux224,'ro',markersize=3,label="Red Seq",markeredgecolor='r',markerfacecolor='None')

	    if self.zcl < 0.6:
		    flux124=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1) & (self.Lir > 10.**lirminloz), self.MV)
		    flux224=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1) & (self.Lir > 10.**lirminloz),self.MU)
	    else:
		    flux124=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1) & (self.Lir > 10.**lirminhiz), self.MV)
		    flux224=N.compress((self.supermembflag > 0.1) & (self.matchflag24 > .1) & (self.Lir > 10.**lirminhiz),self.MU)
		    
	    flux224=flux224-flux124
	    pylab.plot(flux124,flux224,'ro',markersize=3,label="Red Seq",markeredgecolor='r',markerfacecolor='r')

	    if self.zcl < 0.6:
		    zmin=0.4
		    zmax=0.6
		    flux124=N.compress((self.supermembflag < 0.1) & (self.bestz>zmin) & (self.bestz < zmax)& (self.matchflag24 > .1) & (self.Lir > 10.**lirminloz), self.MVbestz)
		    flux224=N.compress((self.supermembflag < 0.1) & (self.bestz>zmin) & (self.bestz < zmax)& (self.matchflag24 > .1) & (self.Lir > 10.**lirminloz),self.MUbestz)
	    else:
		    zmin=0.62
		    zmax=0.8

		    flux124=N.compress((self.supermembflag < 0.1) & (self.bestz>zmin) & (self.bestz < zmax)& (self.matchflag24 > .1) & (self.Lirbestz > 10.**lirminhiz), self.MVbestz)
		    flux224=N.compress((self.supermembflag < 0.1) & (self.bestz>zmin) & (self.bestz < zmax)& (self.matchflag24 > .1) & (self.Lirbestz > 10.**lirminhiz),self.MUbestz)
		    
	    flux224=flux224-flux124
	    pylab.plot(flux124,flux224,'yx',markersize=3,label="Red Seq",markeredgecolor='y',markerfacecolor='y')
	
	    #pylab.axis([-24.,-17.,0.1,1.2])



    def colormagsaltsub(self,dm):#for salt proposal
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    #dm=2.5*N.log10(self.DL/cl1216.DL)

	    specmembra= self.ra[(self.newspecmatchflag > 0.1)]
	    specmembdec= self.dec[(self.newspecmatchflag > 0.1)]
	    specmembv= self.magV[(self.newspecmatchflag > 0.1)]
	    specmembr= self.magR[(self.newspecmatchflag > 0.1)]
	    specmembi= self.magI[(self.newspecmatchflag > 0.1)]
	    compv=[]
	    compr=[]
	    compi=[]
	    delta=3./3600.#3" offset to nearest neighbor
	    for i in range(len(specmembra)):
		    
		    imatch, matchflag,nmatch=findnearest(specmembra[i],specmembdec[i],self.ra,self.dec)
		    if nmatch > 0:
			    if nmatch < 2:
				    compv.append(float(self.magV[imatch]))
				    compi.append(float(self.magI[imatch]))
			    if nmatch > 1:
				    compv.append(float(self.magV[imatch[0]]))
				    compi.append(float(self.magI[imatch[0]]))
				    
	    pylab.plot(specmembv,specmembv-specmembi, 'ko',markersize=6,label="_nolegend_")
	    pylab.plot(compv,compv-compi,'bo',markersize=4,label=r"$\rm [OII]$")
	

	    #dm=2.5*N.log10(self.DL/cl1216.DL)
	    x=N.arange(19.8,24.4,1.)
	    y=-.09*(x-20)+dm-1.6
	    pylab.plot(x,y,'k-',label="_nolegend_")
	    y2=y-.3
	    pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    y2=y+.3
	    pylab.plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    
	    pylab.xticks(N.arange(19,24))
	    pylab.yticks(N.arange(1,4))
	    pylab.axis([18.5,23.5,0.1,3.9])

	    ax=pylab.gca()
	    s=name#+', z=%5.3f, Nsp=%2d, Nph=%2d'%(self.z,len(flux124),len(pflux124))
	    pylab.text(.2,.8,s,fontsize=14,transform=ax.transAxes)



    def colormagmorphallsub(self):
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    dm=2.5*N.log10(self.DL/cl1216.DL)

	    eI=[]
	    eVI=[]
	    sI=[]
	    sVI=[]
	    iI=[]
	    iVI=[]
	    pI=[]
	    pVI=[]
	    for i in range(len(self.magV)):
		    if self.newspecmatchflag[i] > 0:
			    if self.matchflagvistype[i] > 0:
				    if self.matchflag24[i] > 0.:
					    if (self.vistype[i] > -6.) & (self.vistype[i] < 0.):
						    eI.append(self.magI[i])
						    eVI.append(self.magV[i]-self.magI[i])
					    if (self.vistype[i] > 0.) & (self.vistype[i] < 9.):
						    sI.append(self.magI[i])
						    sVI.append(self.magV[i]-self.magI[i])
					    if (self.vistype[i] > 9.) & (self.vistype[i] < 60.):
						    iI.append(self.magI[i])
						    iVI.append(self.magV[i]-self.magI[i])
					    if (self.vistype[i] < -5.):
						    pI.append(self.magI[i])
						    pVI.append(self.magV[i]-self.magI[i])

	#xbins=N.array([-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11],'f')
	#pylab.xticks(ind+width, ('C', 'E', 'S0', 'Sa', '', 'Sb','','Sc','','Sd','','Sm','Im','Irr') ,fontsize=12)

		
	    #print "here's specflagbell ",self.specflagbell
	    #print "here's morphbell ",self.morphbell
	    #print "length of arrays ", len(self.specflagbell),len(self.morphbell),len(self.Ibell)
	    #morphbell=N.compress(self.specflagbell > 0.1,self.morphbell)
	    #Ibell=N.compress(self.specflagbell > 0.1,self.Ibell)
	    #VIbell=N.compress(self.specflagbell > 0.1,self.VIbell)

	    #eI=N.compress(morphbell< 0.1,self.Ibell)
	    #eVI=N.compress(morphbell< 0.1,self.VIbell)
	    #sI=N.compress(abs(morphbell- 1.) < 0.1,self.Ibell)
	    #sVI=N.compress(abs(morphbell-1.)< 0.1,self.VIbell)
	    #iI=N.compress(abs(morphbell-2.)< 0.1,self.Ibell)
	    #iVI=N.compress(abs(morphbell-2.)< 0.1,self.VIbell)
	    #pI=N.compress(abs(morphbell-3.)< 0.1,self.Ibell)
	    #pVI=N.compress(abs(morphbell-3.)< 0.1,self.VIbell)


	    pylab.plot(eI,eVI,'ro',markersize=8)
	    pylab.plot(sI,sVI,'bo',markersize=8)
	    pylab.plot(iI,iVI,'go',markersize=8)
	    pylab.plot(pI,pVI,'co',markersize=8)


	    dm=2.5*N.log10(self.DL/cl1216.DL)
	    x=N.arange(19.8,24.4,1.)
	    y=-.09*(x-20)+2.8+dm
	    pylab.plot(x,y,'k-',label="_nolegend_")
	    
	    pylab.xticks(N.arange(20,24))
	    pylab.yticks(N.arange(1,4))
	    pylab.axis([19.5,23.5,0.25,3.9])

	    s=name+', z=%5.3f'%(self.zcl)
	    pylab.text(19.8,3.2,s,fontsize=14)
#	pylab.axis([0,3.5,14,22])
	#file=str(name)+'colormag.jpg'
	#pylab.savefig(file)

    def colormagallsubspec(self,name,dm):
	flux1=[]
	errflux1=[]
	flux2=[]
	errflux2=[]	
	flux1h=[]
	errflux1h=[]
	flux2h=[]
	errflux2h=[]	
	flux124=[]
	errflux124=[]
	flux224=[]
	errflux224=[]	

	pflux1=[]
	perrflux1=[]
	pflux2=[]
	perrflux2=[]	
	pflux1h=[]
	perrflux1h=[]
	pflux2h=[]
	perrflux2h=[]	
	pflux124=[]
	perrflux124=[]
	pflux224=[]
	perrflux224=[]	

	for i in range(len(self.ra)):
	    #print name,' phot,spec,wmin= ',self.membflag[i],self.newspecmatchflag[i],self.wmin[i]
	    if self.newspecmatchflag[i] > 0:
		    flux2.append(self.magV[i]-self.magI[i])
		    flux1.append(self.magI[i])
		    if self.matchflagha[i] > 0:
			if self.SFflag[i] > 0:
			    flux2h.append(self.magV[i]-self.magI[i])
			    flux1h.append(self.magI[i])

		    if self.matchflag24[i] > 0:
			flux224.append(self.magV[i]-self.magI[i])
			flux124.append(self.magI[i])
			#print "spec+24 ",self.ediscsID[i],N.log10(self.L24[i]),self.vistype[i]
			#s='circle(%12.8f,%12.8f,4.0" \n'%(self.ra[i],self.dec[i])
			#dsfile.write(s)
	    elif (self.supermembflag[i] > 0.):
		    pflux2.append(self.magV[i]-self.magI[i])
		    pflux1.append(self.magI[i])
		    if self.matchflagha[i] > 0:
			if self.SFflag[i] > 0:
			    pflux2h.append(self.magV[i]-self.magI[i])
			    pflux1h.append(self.magI[i])
			    
		    if self.matchflag24[i] > 0:
			pflux224.append(self.magV[i]-self.magI[i])
			pflux124.append(self.magI[i])

	
	flux1=N.array(flux1,'d')
	flux2=N.array(flux2,'d')
	flux1h=N.array(flux1h,'d')
	flux2h=N.array(flux2h,'d')
	flux124=N.array(flux124,'d')
	flux224=N.array(flux224,'d')

	pflux1=N.array(pflux1,'d')
	pflux2=N.array(pflux2,'d')
	pflux1h=N.array(pflux1h,'d')
	pflux2h=N.array(pflux2h,'d')
	pflux124=N.array(pflux124,'d')
	pflux224=N.array(pflux224,'d')
		
	    #pylab.plot(flux1,flux2,'bo')
	title='V-I vs. I'+str(name)
	#pylab.title(title)

	#pylab.plot(pflux1,pflux2, 'wo',markeredgecolor='k',markeredgewidth=1.,markersize=3,label="Photo-z Members")

	#pylab.plot(pflux124,pflux224,'wo',markeredgecolor='r',mew=1.,markersize=5,label="Photo-z + 24um")
	#pylab.plot(pflux1h,pflux2h,'wo',markeredgecolor='b',mew=1.,markersize=3,label="Photo-z + Halpha")

	#pylab.plot(flux1,flux2, 'ko',markersize=16,label="Spec Members")


	pylab.plot(flux1,flux2, 'ko',markeredgecolor='k',markeredgewidth=2,markersize=8,label="Spec Members")
	pylab.plot(flux1h,flux2h,'ws',markersize=22,markeredgecolor='k',markeredgewidth=2,label="Spec +Halpha")
	pylab.plot(flux124,flux224,'k^',markersize=22,label="Spec + 24um")	



	print self.prefix
	print "Nspec+24, Nphot+24 = ",len(flux124),len(pflux124)
	#pylab.legend(loc='lower left')



	x=N.arange(20.2,23.4,1.)
	dm=5.*N.log10(self.DL/cl1216.DL)
	y=-.09*(x-20)+2.8+dm
	pylab.plot(x,y,'k-')
	
	pylab.xticks(N.arange(20,24),fontsize=20.)
	pylab.yticks(N.arange(1,4),fontsize=20.)
	pylab.axis([20.,23.5,0.5,3.8])

#	pylab.axis([0,3.5,14,22])
	#file=str(name)+'colormag.jpg'
	#pylab.savefig(file)


    def veldistallsub(self,dm):
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    print name
	    #dm=2.5*N.log10(self.DL/cl1216.DL)

	    zgalall=[]
	    zgalrs=[]
	    zgalbc=[]
	    zgal24=[]
	    for i in range(len(self.ra)):
		    if (self.newspecmatchflag[i] > 0.1) & (self.matchflag24[i] > -1.):#spec memb & on 24um image
			    zgalall.append(self.specz[i])

			    I=self.magI[i]
			    VI=self.magV[i]-self.magI[i]
			    dmag=(-1.*0.09*I + dm - VI)
			    if self.matchflag24[i] > 0.1:
				    zgal24.append(self.specz[i])
			    else:
				    if (abs(dmag) < .3):
					    zgalrs.append(self.specz[i])

				    if dmag > .3:
					    zgalbc.append(self.specz[i])

	    zgalall=N.array(zgalall,'f')-self.zcl
	    zgalrs=N.array(zgalrs,'f')-self.zcl
	    zgalbc=N.array(zgalbc,'f')-self.zcl
	    zgal24=N.array(zgal24,'f')-self.zcl

	    vcl=my.getvfromz(self.zcl)
	    dvall=N.array((my.getvfromz(zgalall+self.zcl)-vcl)/self.sigma,'d')
	    dvrs=N.array((my.getvfromz(zgalrs+self.zcl)-vcl)/self.sigma,'d')
	    dvbc=N.array((my.getvfromz(zgalbc+self.zcl)-vcl)/self.sigma,'d')
	    dv24=N.array((my.getvfromz(zgal24+self.zcl)-vcl)/self.sigma,'d')

	
	    #print 'zgalall = ',zgalall
	    #xmin=min(zgalall)
	    #xmax=max(zgalall)
	    xmin=-4
	    xmax=4
	    nbin=10
	    x=dvall
	    (xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	    #my.drawhistpylab(xbin,(ybin),'k-')

	    x=dvrs
	    (xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	    my.drawhistpylab(xbin,(ybin),'k-')

	    x=dvbc
	    (xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	    my.drawhistpylab(xbin,(ybin),'b-')

	    x=dv24
	    (xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	    my.drawhistpylab(xbin,(ybin),'r-')

	    pylab.axis([-4.,4.,0,13.])

	    ax=pylab.gca()
	    s=name#+', z=%5.3f, Nsp=%2d, Nph=%2d'%(self.z,len(flux124),len(pflux124))
	    pylab.text(.2,.8,s,fontsize=14,transform=ax.transAxes)


	    try:
		    weird=N.compress(abs(dvall) > 3.,dvall)
		    weirdz=N.compress(abs(dvall) > 3.,zgalall)+self.zcl
		    #print "spec members at > 3 sigma in cluster ",self.prefix
		    #print weird
		    #print 'zcl = ',self.zcl,' redshifts = ',weirdz
	    except:
		    print "NO spec members at > 3 sigma in cluster ",self.prefix
	    #print 'dv24 = ',dv24
	    #print 'v24 - vcl = ',my.getvfromz(zgal24)-vcl
	    return dvall,dvrs,dvbc,dv24

    def lirdistallsub(self):
	    xmin=9.
	    xmax=13.
	    nbin=25
	    lir=[]
	    errlir=[]
	    morphmemb=[]
	    morphmembflag=[]
	    drmemb=[]
	    for i in range(len(self.ra)):
		    if self.supermembflag[i] > 0.:
			    if self.matchflag24[i] > 0:
				    lir.append(self.Lir[i])
				    errlir.append(self.errLir[i])
				    morphmemb.append(self.vistype[i])
				    morphmembflag.append(self.matchflagvistype[i])
				    drmemb.append(self.dr[i])
	    #print lir
	    self.lirmemb=lir
	    self.errlirmemb=errlir
	    self.morphmemb=morphmemb
	    self.morphmembflag=morphmembflag
	    self.drmemb=drmemb
	    x=N.log10(lir)
	    (xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	    #ybin=N.clip(ybin,.001,(max(ybin)+1))
	    my.drawhistpylab(xbin,(ybin),'k-')
	    #min=N.log10(3.*errlir[1])
	    min=N.log10(self.Lir80)
	    #for i in range(len(lir)):
		#    print i,"errlir=",errlir[i],lir[i],min
	    pylab.axvline(x=min,color='0.4',ls=':')
	    pylab.axvline(x=lirmin,color='k',ls='--')
	    #pylab.axis([9.25,12.4,1.,18.])
	    pylab.axis([9.8,12.4,1.,15.])
	    ax=pylab.gca()
	    pylab.xticks(fontsize=14)
	    #ax.set_yscale('log')
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    s=name#+', z=%5.3f, N=%3i'%(self.z,sum(ybin))
	    pylab.text(.1,.8,s,fontsize=14,transform=ax.transAxes)

    def lirMrallsub(self):
	    xmin=9.
	    xmax=13.
	    nbin=25
	    lir=[]
	    errlir=[]
	    mr=[]
	    errmr=[]
	    lirs=[]
	    errlirs=[]
	    mrs=[]
	    errmrs=[]
	    s=self.prefix+'.irlf.dat'
	    output00=open(s,'w')
	    for i in range(len(self.ra)):
		    if self.supermembflag[i] > 0.:
			    if self.matchflag24[i] > 0:
				    lir.append(self.Lir[i])
				    errlir.append(self.errLir[i])
				    mr.append(self.MV[i])
				    string1 = "%s  %5.1f  %9.4f %5.1f %12.8f %12.8f %8.2f %8.2f %8.4e %8.4e %8.4e %8.4e %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n"%(self.ediscsID[i],self.newspecmatchflag[i],self.specz[i],self.spectype[i], self.ra[i],self.dec[i],self.vistype[i],self.matchflag24[i],self.Lir[i],self.f24c[i],self.errf24c[i],self.SFRir[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i], self.magI[i],self.mageIapsim[i], self.magJ[i],self.mageJapsim[i])
				    output00.write(string1)

				    if self.newspecmatchflag[i] > 0.:
					    lirs.append(self.Lir[i])
					    errlirs.append(self.errLir[i])
					    mrs.append(self.MV[i])
				    #errmr.append(self.errMRlo[i])
	    self.lirsmemb=N.array(lirs,'d')
	    self.lirmemb=N.array(lir,'d')
		
	    ax=pylab.gca()
	    ax.yaxis.tick_left()
			    
	    pylab.plot(mr,lir,'wo',markersize=6, label='photoz')
	    pylab.plot(mrs,lirs,'ko',markersize=8, label='spec')
	    #min=(3.*errlir[1])
	    min=self.Lir80
	    pylab.axhline(y=min,color='0.4',ls=':',label='_nolegend_')
	    pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_')
	    #pylab.axis([9.25,12.4,1.,18.])
	    #pylab.axis([-24.2,-16.2,5.e9,2.e13])
	    pylab.axis([-24.2,-16.2,6.e9,4.e12])
	    pylab.xticks(N.arange(-24,-16),(r'$-24$','',r'$-22$','',r'$-20$','',r'$-18$',''))
	    ax=pylab.gca()
	    ax.set_yscale('log')
	    name='CL'+str(self.idra)+'-'+str(self.iddec)
	    s=name#+', z=%5.3f, N=%3i'%(self.z,len(lir))
	    #s=name+', z=%5.3f'%(self.z)
	    pylab.text(.5,.8,s,fontsize=14,transform=ax.transAxes,horizontalalignment='center')



    def gotoitHa(self,catalogpath,prefix,ra,dec,speccat,morphidcat,mipsimage):
	    if mode < .1:
		    self.readoptIRphot(catalogpath+'photometry/cl'+str(ra)+'-'+str(dec)+'vrijk_optIRphot.v23.cat')
		    self.readphotoz(catalogpath+'photoz/v241/zcat.final.cl'+str(ra)+'-'+str(dec)+'vrijk.V2.4.1.dat')
		    self.read24('/home/rfinn/research/clusters/spitzer/MasterTables/'+str(prefix)+'mosaic_extract_final.tbl',mipsimage)
		    self.readirac(prefix)

		    self.readhafile(catalogpath+'hafiles/'+str(prefix)+'sfrtableapj.dat')
		    self.readgimmorphology(catalogpath+'morphology/HST/GIM2D/'+str(morphidcat),catalogpath+'morphology/HST/GIM2D/cl'+str(ra)+'-'+str(dec)+'_drz_sci.struct2.cat')
		    self.readvismorphology(catalogpath+'morphology/HST/Visual/cl'+str(ra)+'_'+str(dec)+'.f.cat')
		    self.readspectroscopy(catalogpath+'spectroscopy/'+str(speccat))
		    self.memb()
		    self.magnitudetofluxconversion(prefix,1.)
		    self.mastertable(fullprefix)
	    if (mode-1.) < .1:
		    self.mode1()

    def gotoit(self):
	    if mode < .1:
		    prefix=self.prefix
		    ra=self.idra
		    dec=self.iddec
		    cluster=str(ra)+'-'+str(dec)
		    s=catalogpath+'photometry/v23/cl'+cluster+'*optIRphot*'
		    files=glob.glob(s)
		    print s
		    print files
		    phot=files[0]
		    s=catalogpath+'photoz/v241/zcat.final.cl'+cluster+'*'
		    files=glob.glob(s)
		    photz=files[0]
		    self.readoptIRphot(phot)
		    #print "Length of ediscsID = ",len(self.ediscsID)
		    self.readphotoz(photz)
		    s='/home/rfinn/research/clusters/spitzer/MasterTables/'+str(prefix)+'mosaic_extract_final.tbl'
		    self.read24(s,self.mipsimage)
		    self.readirac(prefix)
		    s=catalogpath+'morphology/HST/GIM2D/cl'+cluster+'*'
		    files=glob.glob(s)
		    if len(files) > 0:
			    morphcat=files[0]
			    s=catalogpath+'morphology/HST/GIM2D/'+prefix+'*.id'
			    files=glob.glob(s)
			    #print s
			    #print files
			    morphid=files[0]

			    self.readgimmorphology(morphid,morphcat)
		    s=catalogpath+'morphology/HST/Visual/cl'+str(ra)+'_'+str(dec)+'.f.cat'
		    files=glob.glob(s)
		    if len(files) > 0:
			    print "READING HST VISUAL MORPHOLOGIES"
			    self.readvismorphology(catalogpath+'morphology/HST/Visual/cl'+str(ra)+'_'+str(dec)+'.f.cat')
			    
		    s=catalogpath+'spectroscopy/final/'+prefix+'*'
		    files=glob.glob(s)
		    speccat=files[0]
		    if speccat.find('v1') > -1:
			    self.readspectroscopy(speccat)
		    if speccat.find('v24') > -1:
			    self.readspectroscopyv2(speccat)
					    
		    self.memb()
		    self.magnitudetofluxconversion(prefix,0.)
		    self.mastertable(fullprefix)

	    if mode > 0.1:
		    self.mode1()
    def mode1(self):
	    self.readmaster()
	    #self.measurenearest()
	    self.LIRhistsub('ko')
	    #self.makeplots()

    def makeplots(self):
	#self.plothavs24(prefix)
	#self.plothavs24morph(prefix)
	self.count()
	self.colormag()
	#self.plotston(prefix)
	#self.hist(prefix)
	#self.morphhist(prefix)
	#self.plot24vsEWOII(prefix)

    def writeapjtbl(self,file):
	    #name z sigma R200 (Mpc) R200 (arcmin) Nspec24 Nspec+phot24 Nmorph24
	    #name='CL'+self.idra+'$-$'+self.iddec
	    name=EdiscsName[self.prefix]
	    errf80=ErrorF80[self.prefix]
	    self.n24spec=0
	    self.n24phot=0
	    self.n24spec2=0
	    self.n24phot2=0
	    for i in range(len(self.supermembflag)):
		    if self.supermembflag[i] > 0.:
			    if self.matchflag24[i] > 0:
				    if self.newspecmatchflag[i] > 0.1:
					    self.n24spec += 1
					    if self.dr[i] < 1.:
						    self.n24spec2 += 1
				    else:
					    self.n24phot += 1
					    if self.dr[i] < 1.:
						    self.n24phot2 += 1
	    self.n24memb=self.n24spec+self.n24phot
	    self.n24memb2=self.n24spec2+self.n24phot2
		    
	    #s='%s & %6.4f & %3i & %3i & %5.2f & %4.1f & %i  &%i & %i & %5.2f &%i & %i & %i\\\\ \n'%(name,self.zcl,self.sigma,self.f80,N.log10(self.Lir80),self.SFR80,self.n24spec,self.n24phot,self.n24memb,self.r200arcmin,self.n24spec2,self.n24phot2,self.n24memb2)

	    s='%s & %6.4f & %3i$^{+%i}_{-%i}$ & %3i  & %3i  & %3i$\pm$%s & %5.2f & %4.1f \\\\ \n'%(name,self.zcl,self.sigma,self.errsigmap,self.errsigmam,self.fmin,self.fmax,self.f80,errf80,N.log10(self.Lir80),self.SFR80)#
	    file.write(s)
	    return self.n24spec,self.n24phot,self.n24spec2,self.n24phot2

    def IRoptdist(self):
	    distall=[]
	    print "IRoptdist test", self.prefix
	    #print self.matchflag24 > .1
	    print len(self.ra24),len(self.matchflag24)
	    ra=[]
	    dec=[]
	    for i in range(len(self.matchflag24)):#create array of 24um ra and dec
		    if self.matchflag24[i] > .1:
			    ra.append(self.ra24[i])
			    dec.append(self.dec24[i])

	    ra=N.array(ra,'f')
	    dec=N.array(dec,'f')
	    #ra=self.ra24[self.matchflag24 > .1]#get ra for 24um sources
	    #dec=self.dec24[self.matchflag24 > .1]#get ra for 24um sources
	    for i in range(len(ra)):
		    d=pylab.sqrt((ra[i]-self.ra)**2+(dec[i]-self.dec)**2)*3600.#calc dist in arcsec to ediscs sources
		    d=d[d<6.]
		    for dist in d:
			    distall.append(dist)#keep only distances w/in 20"
	    #make histogram of all distances
	    nbin=40
	    xmin=0.
	    xmax=11.
	    (xbin,ybin,ybinerr) = my.horizontalhist(distall,xmin,xmax,nbin)
	    #ybin=N.clip(ybin,.001,(max(ybin)+1))
	    area=N.zeros(nbin,'f')
	    dx=xbin[1]-xbin[0]
	    for i in range(nbin):
		    area[i]=N.pi*((xbin[i]+dx)**2-(xbin[i])**2)
	    ybin=ybin/area#normalize counts by surface area
	    my.drawhistpylab(xbin,(ybin),'k-')
	    #pylab.hist(distall,nbin)
	    pylab.axis([0,5.7,0.,48])
	    ax=pylab.gca()
	    name='CL'+str(self.idra)+'-'+str(self.iddec)

	    s=name#+', z=%5.3f, $\sigma$=%3i'%(self.z,self.sigma)
	    pylab.text(.2,.82,s,fontsize=14,transform=ax.transAxes,horizontalalignment='left')

    def lirstellmass(self,color,lirlim):
	    mmass=minmass
	    if color.find('all') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.Lir)
	    elif color.find('red') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.Lir)
	    elif color.find('blue') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.Lir)
		    #mmass=minmassblue
	    x=pylab.array(x,'d')
	    y=pylab.array(y,'d')
	    pylab.plot(x,y,'ro')
	    if self.zcl > 0.6:
		    if color.find('blue') > -1:
			    x1=x[(y>10.**lirlim) & (x > minmasshizblue)]
			    y1=y[(y>10.**lirlim) & (x > minmasshizblue)]
		    else:
			    x1=x[(y>10.**lirlim) & (x > minmass)]
			    y1=y[(y>10.**lirlim) & (x > minmass)]
	    else:
		    if color.find('blue') > -1:
			    x1=x[(y>10.**lirlim) & (x > minmasslozblue)]
			    y1=y[(y>10.**lirlim) & (x > minmasslozblue)]
		    else:
			    x1=x[(y>10.**lirlim) & (x > minmassloz)]
			    y1=y[(y>10.**lirlim) & (x > minmassloz)]
	    x1=x1.tolist()
	    y1=y1.tolist()
	    return x1,y1
	    #pylab.plot(x[(y <= 10.**lirlim)],y[(y <= 10.**lirlim)],'r.')

    def ssfrstellmass(self,color,lirlim):
	    mmass=minmass
	    if color.find('all') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.SFRir)
	    elif color.find('red') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.SFRir)
	    elif color.find('blue') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.SFRir)
		    #mmass=minmassblue
	    x=pylab.array(x,'d')
	    y=pylab.array(y,'d')
	    pylab.plot(x,y/x,'ro')
	    if self.zcl > 0.6:
		    if color.find('blue') > -1:
			    x1=x[(y>10.**lirlim) & (x > minmasshizblue)]
			    y1=y[(y>10.**lirlim) & (x > minmasshizblue)]
		    else:
			    x1=x[(y>10.**lirlim) & (x > minmass)]
			    y1=y[(y>10.**lirlim) & (x > minmass)]
	    else:
		    if color.find('blue') > -1:
			    x1=x[(y>10.**lirlim) & (x > minmasslozblue)]
			    y1=y[(y>10.**lirlim) & (x > minmasslozblue)]
		    else:
			    x1=x[(y>10.**lirlim) & (x > minmassloz)]
			    y1=y[(y>10.**lirlim) & (x > minmassloz)]
	    x1=x1.tolist()
	    y1=y1.tolist()
	    return x1,y1
	    #pylab.plot(x[(y <= 10.**lirlim)],y[(y <= 10.**lirlim)],'r.')

    def lirstellmass2(self,color,lirlim):#for plotting lir vs Mv for red and blue
	    if color.find('all') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1),self.Lir)
	    elif color.find('red') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & (self.redflag),self.Lir)
	    elif color.find('blue') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 >0.1) & ~(self.redflag),self.Lir)

	    pylab.plot(x,y,'ro')

	    x=pylab.array(x,'d')
	    y=pylab.array(y,'d')
	    if self.zcl > 0.6:
		    x1=x[(y>10.**lirlim) & (x < Mvcuthz)]
		    y1=y[(y>10.**lirlim) & (x < Mvcuthz)]
	    else:
		    x1=x[(y>10.**lirlim) & (x < Mvcut)]
		    y1=y[(y>10.**lirlim) & (x < Mvcut)]
	    x1=x1.tolist()
	    y1=y1.tolist()
	    return x1,y1

    def fraclirstellmass(self,color,lirlim):
	    mmass=minmass
	    if color.find('all') > -1:
		    x=N.compress((self.supermembflag >0.1) &(self.matchflag24 >-0.1),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1),self.Lir)

	    elif color.find('red') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & (self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & (self.redflag),self.Lir)
	    elif color.find('blue') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & ~(self.redflag),self.stellmass)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & ~(self.redflag),self.Lir)
		    #mmass=minmassblue
	    x=pylab.array(x,'d')
	    y=pylab.array(y,'d')
	    pylab.plot(x,y,'ro')
	    if self.zcl > 0.6:
		    y1=(y>10.**lirminhiz)
	    else:
		    y1=(y>10.**lirminloz)

	    x1=x
	    x1=x1.tolist()
	    y1=y1.tolist()
	    return x1,y1


    def fraclirMv(self,color,lirlim):
	    mmass=minmass
	    if color.find('all') > -1:
		    x=N.compress((self.supermembflag >0.1) &(self.matchflag24 >-0.1),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1),self.Lir)

	    elif color.find('red') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & (self.redflag),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & (self.redflag),self.Lir)
	    elif color.find('blue') > -1:
		    x=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & ~(self.redflag),self.MV)
		    y=N.compress((self.supermembflag >0.1)&(self.matchflag24 > -0.1) & ~(self.redflag),self.Lir)
		    #mmass=minmassblue
	    x=pylab.array(x,'d')
	    y=pylab.array(y,'d')
	    #pylab.plot(x,y,'ro')
	    if self.zcl > 0.6:
		    y1=(y>10.**lirminhiz)
	    else:
		    y1=(y>10.**lirminloz)

	    x1=x
	    x1=x1.tolist()
	    y1=y1.tolist()
	    return x1,y1



def getgroups():
	#subroutine to get group members for bianca's groups based on (1) cluster image, (2) group redshift, and (3) group velocity dispersion
	lirall=[]
	zall=[]
	matchflag24all=[]
	OIIgroupall=[]
	mvall=[]
	ediscsIDall=[]
	raall=[]
	decall=[]
	clusters=[cl1040,cl1103,cl105411,cl1040,cl1301,cl1103,cl105412]
	groupz=[.7798,.7031,.6130,.6316,.3969,.6261,.7305]
	groupsigma=[259.,252.,227.,179.,391.,336.,182.]
	o2flag=[1,1,1,1,0,0,0,0]
	for i in range(len(clusters)):
		(z,lir,matchflag24,OIIgroup,mv,ediscsID,ra,dec)=getgroupsub(clusters[i],groupz[i],groupsigma[i],o2flag[i])
	#in addition to "clusters" cl1037 and cl1420
		lirall=lirall+lir
		zall=zall+z
		matchflag24all=matchflag24all+matchflag24
		OIIgroupall=OIIgroupall+OIIgroup
		mvall=mvall+mv#already making Mv cut in getgroupsub()
		ediscsIDall=ediscsIDall+ediscsID
		raall=raall+ra
		decall=decall+dec
	#print lirall
	lirall=N.array(lirall,'d')
	zall=N.array(zall,'f')
	matchflag24all=N.array(matchflag24all,'f')
	OIIgroupall=N.array(OIIgroupall,'f')
	mvall=N.array(mvall,'f')
	raall=N.array(raall,'f')
	decall=N.array(decall,'f')
	#want to measure fraction of LIRGS, keeping note of which groups have high OII fractions
	#so generate combined list of z, LIR, match24flag,OIIgroupflag
	z=N.compress(mvall < Mvcut,zall)
	print "average redshift of group members = %6.3f +/- %6.3f"%(N.average(z),pylab.std(z))
	lirmin=10.**10.95
	nlirg=len(N.compress((matchflag24all > .1)&(lirall>lirmin)&(mvall < Mvcut),zall))
	ntot=len(N.compress((mvall < Mvcut),zall))
	(a,b,c)=my.ratioerror(nlirg,ntot)
	print "Fraction of LIRGS in groups = %5.2f -%5.2f + %5.2f"%(a,b,c) 

	nlirg=len(N.compress((matchflag24all > .1)&(lirall>lirmin)&(OIIgroupall > .1)&(mvall < Mvcut),zall))

	ntot=len(N.compress((OIIgroupall > 0.1)&(mvall < Mvcut),zall))
	(a,b,c)=my.ratioerror(nlirg,ntot)
	print "Fraction of LIRGS in groups w/high OII fraction = %5.2f -%5.2f + %5.2f (%i/%i)"%(a,b,c,nlirg,ntot)
	index=N.arange(len(zall))
	ids=N.compress((matchflag24all > .1)&(lirall>lirmin)&(OIIgroupall > .1)&(mvall < Mvcut),index)
	IDs=[]
	for i in ids:
		IDs.append(ediscsIDall[i])
	print "Ediscs IDs of LIRGS in high OII groups are "
	print IDs

	outfile2=open('GroupLirgs.dat','w')
	outfile=open('GroupLirgs.reg','w')
	outfile.write("global color=red font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	for i in range(len(raall)):
		s="%s \n"%(ediscsIDall[i])
		outfile2.write(s)
		if i in ids:
			string1 = "circle(%12.8f, %12.8f, 2\")# color=blue \n"%(raall[i],decall[i])
		else:
			string1 = "circle(%12.8f, %12.8f, 2\")# \n"%(raall[i],decall[i])
		outfile.write(string1)
	outfile.close()
	outfile2.close()

	nlirg=len(N.compress((matchflag24all > .1)&(lirall>lirmin)&(OIIgroupall < .1)&(mvall < Mvcut),zall))
	ntot=len(N.compress((OIIgroupall < 0.1)&(mvall < Mvcut),zall))
	(a,b,c)=my.ratioerror(nlirg,ntot)
	print "Fraction of LIRGS in groups w/low OII fraction = %5.2f -%5.2f+ %5.2f (%i/%i)"%(a,b,c,nlirg,ntot) 



def getgroupsub(cluster,gz,gsigma,o2flag):
	lir=[]
	zspec=[]
	matchflag24=[]
	OIIgroup=[]#1 if high OII fraction
	mv=[]
	ediscsID=[]
	ra=[]
	dec=[]
	vcl=my.getvfromz(gz)
	vmax=vcl+3.*gsigma
	zmax=my.getzfromv(vmax)
	vmin=vcl-3.*gsigma
	zmin=my.getzfromv(vmin)
	gMvcut=calcMvcut(gz)
	for i in range(len(cluster.fieldz)):
		if cluster.fieldz[i] > zmin:
			if cluster.fieldz[i]< zmax:#then append Lir, z, matchflag24, OIIgroupflag
				if cluster.fieldMv[i] < gMvcut:
					matchflag24.append(cluster.fieldmatchflag24[i])
					lir.append(cluster.fieldLir[i])
					zspec.append(cluster.fieldz[i])
					OIIgroup.append(o2flag)
					mv.append(cluster.fieldMv[i])
					ediscsID.append(cluster.fieldediscsID[i])
					ra.append(cluster.fieldra[i])
					dec.append(cluster.fielddec[i])
	return zspec,lir,matchflag24,OIIgroup,mv,ediscsID,ra,dec


def plotcolormagall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.colormagallsub(4.403)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	prop=matplotlib.font_manager.FontProperties(size=12)
	pylab.legend(loc='lower left',numpoints=1,prop=prop,markerscale=0.6)
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,2)
	cl1354.colormagallsub(4.434)
	#pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.colormagallsub(4.485)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.colormagallsub(4.209)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,5)
	cl105411.colormagallsub(4.260)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	
	pylab.subplot(4,4,6)
	cl1227.colormagallsub(4.177)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('colormagalla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(4,4,7)
	#cl1103.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1353.colormagallsub(4.070)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1037.colormagallsub(4.053)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,9)
	cl1232.colormagallsub(3.980)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1411.colormagallsub(3.904)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	pylab.subplot(4,4,11)
	cl1420.colormagallsub(3.909)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('colormagallb.eps')

	#pylab.cla()
	#pylab.clf()
	pylab.subplot(4,4,12)
	cl1301.colormagallsub(3.889)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,13)
	cl1138.colormagallsub(3.845)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,14)
	cl1018.colormagallsub(3.855)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,15)
	cl1059.colormagallsub(3.826)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,16)
	cl1202.colormagallsub(3.742)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	ax=pylab.gca()
	pylab.text(-1.,-.5,r'$\rm I-band \ Magnitude$',fontsize=28,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.3,2.,r'$\rm V-I$',fontsize=28,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('colormagallc.eps')

def plotcolormagallBV():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.colormagallBVsub(4.403)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	prop=matplotlib.font_manager.FontProperties(size=12)
	pylab.legend(loc='lower left',numpoints=1,prop=prop,markerscale=0.6)
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,2)
	cl1354.colormagallBVsub(4.434)
	#pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.colormagallBVsub(4.485)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.colormagallBVsub(4.209)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,5)
	cl105411.colormagallBVsub(4.260)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	
	pylab.subplot(4,4,6)
	cl1227.colormagallBVsub(4.177)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('colormagalla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(4,4,7)
	#cl1103.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1353.colormagallBVsub(4.070)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1037.colormagallBVsub(4.053)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,9)
	cl1232.colormagallBVsub(3.980)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1411.colormagallBVsub(3.904)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	pylab.subplot(4,4,11)
	cl1420.colormagallBVsub(3.909)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('colormagallb.eps')

	#pylab.cla()
	#pylab.clf()
	pylab.subplot(4,4,12)
	cl1301.colormagallBVsub(3.889)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,13)
	cl1138.colormagallBVsub(3.845)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,14)
	cl1018.colormagallBVsub(3.855)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,15)
	cl1059.colormagallBVsub(3.826)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,16)
	cl1202.colormagallBVsub(3.742)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	ax=pylab.gca()
	pylab.text(-1.,-.5,r'$\rm V-band \ Magnitude$',fontsize=28,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.3,2.,r'$\rm B-V$',fontsize=28,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('colormagallBV.eps')

def plotcolormagallUV():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5

	xmin=-26.
	xmax=-15
	ymin=-1
	ymax=2.5
	pylab.cla()
	pylab.clf()

	pylab.subplot(2,1,1)
	x=N.compress((gems.hzwhmips > 1.1),gems.hzgMv)
	x2=N.compress((gems.hzwhmips > 1.1),gems.hzgMu)

	pylab.plot(x,x2-x,'.',markeredgecolor='c',markersize=1,markerfacecolor='None')
	#x=gems.gMv[pylab.where(gems.whmips < 1.1)]
	#x2=gems.gMu[pylab.where(gems.whmips < 1.1)]
	x=N.compress((gems.hzwhmips < 1.1),gems.hzgMv)
	x2=N.compress((gems.hzwhmips < 1.1),gems.hzgMu)
	pylab.plot(x,x2-x,'o',markeredgecolor='b',markerfacecolor='None',markersize=3)
	for i in range(6):
		if i == 0:
			cl=c1
		if i == 1:
			cl=c2
		if i == 2:
			cl=c3
		if i == 3:
			cl=c4
		if i == 4:
			cl=c5
		if i == 5:
			cl=c6

		cl.colormagallUVsub(4.403)

	xl=pylab.arange(xmin+1,xmax)
	yl=1.15-0.31*(0.79)-0.08*(xl+20)
	pylab.plot(xl,yl+.25,'k-',lw=1)
	pylab.plot(xl,yl,'k--',lw=1)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	pylab.axis([xmin,xmax,ymin,ymax])

	pylab.subplot(2,1,2)
	x=N.compress((gems.lzwhmips > 1.1),gems.lzgMv)
	x2=N.compress((gems.lzwhmips > 1.1),gems.lzgMu)
	pylab.plot(x,x2-x,'.',markeredgecolor='c',markersize=1,markerfacecolor='None')

	x=N.compress((gems.lzwhmips < 1.1),gems.lzgMv)
	x2=N.compress((gems.lzwhmips < 1.1),gems.lzgMu)
	pylab.plot(x,x2-x,'o',markeredgecolor='b',markerfacecolor='None',markersize=3)

	for i in range(6,16):
		if i == 6:
			cl=c7
		if i == 7:
			cl=c8
		if i == 8:
			cl=c9
		if i == 9:
			cl=c10
		if i == 10:
			cl=c11
		if i == 11:
			cl=c12
		if i == 12:
			cl=c13
		if i == 13:
			cl=c14
		if i == 14:
			cl=c15
		if i == 15:
			cl=c16
		if i == 16:
			cl=c17
		cl.colormagallUVsub(4.403)
	pylab.legend([r'$\rm GEMS$',r'$\rm GEMS \ IR$',r'$\rm EDisCS$',r'$\rm EDisCS \ IR$',r'$\rm EDi \ L_{IR}>min$'],loc='lower left',numpoints=1)

	xl=pylab.arange(xmin+1,xmax)
	yl=1.15-0.31*(0.6)-0.08*(xl+20)
	pylab.plot(xl,yl+.25,'k-',lw=1)
	pylab.plot(xl,yl,'k--',lw=1)

	pylab.axis([xmin,xmax,ymin,ymax])


	ax=pylab.gca()
	pylab.text(.5,-.2,r'$\rm M_V$',fontsize=28,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-.2,1.,r'$\rm M_U-M_V$',fontsize=28,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('colormagallUV.eps')

def plotcolormagsalt():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.95,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.colormagsaltsub(4.403)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	prop=matplotlib.font_manager.FontProperties(size=12)
	pylab.legend(loc='lower left',numpoints=1,prop=prop,markerscale=0.6)
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,2)
	cl1354.colormagsaltsub(4.434)
	#pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.colormagsaltsub(4.485)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.colormagsaltsub(4.209)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,5)
	cl105411.colormagsaltsub(4.260)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	
	pylab.subplot(4,4,6)
	cl1227.colormagsaltsub(4.177)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('colormagalla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(4,4,7)
	#cl1103.colormagallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1037.colormagallsub(4.053)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1232.colormagallsub(3.980)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,9)

	cl1138.colormagallsub(3.845)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1103.colormagallsub(3.742)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	ax=pylab.gca()
	pylab.text(-1.,-.5,r'$\rm I-band \ Magnitude$',fontsize=28,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.3,2.,r'$\rm V-I$',fontsize=28,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('colormagSALT.eps')


def plotveldistall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	(a,b,c,d)=cl1216.veldistallsub(4.403)

	dvall=a
	dvrs=b
	dvbc=c
	dv24=d

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,2)
	(a,b,c,d)=cl1354.veldistallsub(4.434)
	#pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	print 'dvall',len(dvall),len(a)
	dvall=num.append(dvall,a)
	print 'after appending ',len(dvall)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)

	pylab.subplot(4,4,3)
	(a,b,c,d)=cl105412.veldistallsub(4.485)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)



	pylab.subplot(4,4,4)
	(a,b,c,d)=cl1040.veldistallsub(4.209)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,5)
	(a,b,c,d)=cl105411.veldistallsub(4.260)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,6)
	(a,b,c,d)=cl1227.veldistallsub(4.177)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,7)
	(a,b,c,d)=cl1353.veldistallsub(4.070)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,8)
	(a,b,c,d)=cl1037.veldistallsub(4.053)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,9)
	(a,b,c,d)=cl1232.veldistallsub(3.980)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,10)
	(a,b,c,d)=cl1411.veldistallsub(3.904)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)



	pylab.subplot(4,4,11)
	(a,b,c,d)=cl1420.veldistallsub(3.909)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)

	pylab.subplot(4,4,12)
	(a,b,c,d)=cl1301.veldistallsub(3.889)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,13)
	(a,b,c,d)=cl1138.veldistallsub(3.845)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)



	pylab.subplot(4,4,14)
	(a,b,c,d)=cl1018.veldistallsub(3.855)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,15)
	(a,b,c,d)=cl1059.veldistallsub(3.826)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)


	pylab.subplot(4,4,16)
	(a,b,c,d)=cl1202.veldistallsub(3.742)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	dvall=num.append(dvall,a)
	dvrs=num.append(dvrs,b)
	dvbc=num.append(dvbc,c)
	dv24=num.append(dv24,d)

	ax=pylab.gca()
	pylab.text(-1.,-.5,r'$\rm \Delta v/\sigma$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.5,2.,r'$\rm N_{gal}$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('veldistall.eps')


	pylab.clf()
	pylab.cla()
	xmin=-3.
	xmax=3.
	nbin=24
	x=dvall
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	#my.drawhistpylab(xbin,(ybin),'k-')
	dbin=xbin[1]-xbin[0]

	x=dvrs
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybin=ybin/N.sum(ybin)
	my.drawhistpylab(xbin,(ybin),'k-')
	ybinerr=ybinerr/N.sum(ybin)
	#pylab.errorbar(xbin+dbin/2.,ybin,yerr=ybinerr,fmt=None,ecolor='k')
	

	x=dvbc
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybin=ybin/N.sum(ybin)
	my.drawhistpylab(xbin,(ybin),'b-')
	ybinerr=ybinerr/N.sum(ybin)
	#pylab.errorbar(xbin+dbin/2.,ybin,yerr=ybinerr,fmt=None,ecolor='b')
	
	x=dv24
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybin=ybin/N.sum(ybin)
	my.drawhistpylab(xbin,(ybin),'r-')
	ybinerr=ybinerr/N.sum(ybin)
	#pylab.errorbar(xbin+dbin/2.,ybin,yerr=ybinerr,fmt=None,ecolor='r')

	pylab.axis([-3.1,3.1,0.,.17])
	pylab.xlabel(r'$\Delta v/\sigma$',fontsize=24)
	pylab.ylabel(r'$\rm N_{gal}/\Sigma N_{gal}$',fontsize=24)
	ax=pylab.gca()
	pylab.text(.1,.9,'Red Sequence',color='k',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.1,.85,'Blue Cloud',color='b',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.1,.8,'IR Sources',color='r',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.savefig('dvsigma.eps')
	(a,b)=stats.ks_2samp(dv24,dvrs)
	print '2-sample KS test comparing vel dist of 24 w/RS D-value = %6.4f p-value=%6.4f'%(a,b)
	(a,b)=stats.ks_2samp(dv24,dvbc)
	print '2-sample KS test comparing vel dist of 24 w/BC D-value = %6.4f p-value=%6.4f'%(a,b)
	(a,b)=stats.ks_2samp(dvrs,dvbc)
	print '2-sample KS test comparing vel dist of RS w/BC D-value = %6.4f p-value=%6.4f'%(a,b)

	pylab.clf()
	pylab.cla()
	dv24=pylab.array(dv24,'f')
	dvbc=pylab.array(dvbc,'f')
	dvrs=pylab.array(dvrs,'f')

	(a,b)=my.cumulative(dv24[pylab.where(abs(dv24)<3)])
	pylab.plot(a,b,'r-')
	(a,b)=my.cumulative(dvbc[pylab.where(abs(dvbc) < 3)])
	#print len(a),len(b),a
	#print b
	try:
		pylab.plot(a,b,'b-')
	except:
		print "Error plotting cumulative distribution of blue cloud members"
	(a,b)=my.cumulative(dvrs[pylab.where(abs(dvrs)<3)])
	pylab.plot(a,b,'k-')
	pylab.axis([-3.1,3.1,-.01,1.01])
	pylab.xlabel(r'$\Delta v/\sigma$',fontsize=36)
	pylab.ylabel(r'$\rm Cumulative \ Distribution$',fontsize=36)
	ax=pylab.gca()
	pylab.text(.1,.9,'Red Sequence',color='k',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.1,.85,'Blue Cloud',color='b',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.1,.8,'IR Sources',color='r',fontsize=20,horizontalalignment='left',transform=ax.transAxes)
	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)


	pylab.savefig('dvsigmacum.eps')

def plotcolormagmorphall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=23.5
	xlaby=-1.
	ylabx=18.75
	ylaby=5.5
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(421)
	cl1216.colormagmorphallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(422)
	cl1354.colormagmorphallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(423)
	cl105412.colormagmorphallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(424)
	cl1040.colormagmorphallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(425)
	cl105411.colormagmorphallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(426)
	cl1227.colormagmorphallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(427)
	cl1037.colormagmorphallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(428)
	cl1232.colormagmorphallsub()
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.text(0.,-0.5,'I-band Magnitude',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-1.2,2.,'V-I',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,'I-band Magnitude',fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,'V-I',fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('colormagmorphall.eps')

def plotLIRLI():
    pylab.cla()
    pylab.clf()
    cl1216.LIRLIsub('bo')
    cl1354.LIRLIsub('ro')
    cl105412.LIRLIsub('go')
    cl1040.LIRLIsub('ko')
    cl105411.LIRLIsub('yo')
    cl1227.LIRLIsub('co')
    #cl1216.LIRLIsub('bo')
    pylab.xlabel(r'$\rm{L(I)/L_\odot}$',fontsize=24)
    pylab.ylabel(r'$\rm{L(IR)/L_\odot}$',fontsize=24)
    pylab.legend(loc='upper left')
    x=N.arange(1.,2.e12,1.e11)
    y=1.e11*N.ones(len(x),'d')
    pylab.plot(x,y,'k--')
    pylab.text(2.e8,1.1e11,'LIRGS',fontsize=18)
    pylab.axis([1.e8,1.e12,9.e9,1.e12])
    ax=pylab.gca()
    ax.set_yscale("log")
    ax=pylab.gca()
    ax.set_xscale("log")

    pylab.savefig('LIRLI.eps')





#cl0023 = ediscs24()
#cl0023.read24('/home/Corey/spitzer/cl0023/mips/24/r13809664/ch1/pbcd/output/mosaic_extract.tbl')
#cl0023.readha2file('/home/Corey/spitzer/cl0023/mips/24/r13809664/ch1/pbcd/output/cl0023sfrtable2.dat')





def plotlirdistall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=41.5
	xlaby=-.3
	ylabx=1.3
	ylaby=1.5
	xlabel='$\rm log_{10}(L_{IR}/L_\odot)$'
	ylabel='$\rm log_{10}(N_{gal})$'
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.lirdistallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	locs,labels=pylab.xticks()
	yloc=1.2*N.ones(len(locs),'f')
	sfrlab=[]
	for i in range(len(labels)):
		s="%3.1f" % (10.**float(locs[i])*bellconv*Lsol)
		s=str(s)
		print locs[i],s
		pylab.text(float(locs[i]),1.05,s,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
		     
	pylab.subplot(4,4,2)
	cl1354.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,5)
	cl105411.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,6)
	cl1227.lirdistallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirdistalla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(337)
	#cl1103.lirdistallsub()

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1353.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1037.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.text(-0.5,xlaby,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(-2.3,ylaby,r'$\rm log_{10}(N_{gal})$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	#pylab.savefig('lirdistalla.eps')
	#pylab.cla()
	#pylab.clf()

	pylab.subplot(4,4,9)
	cl1232.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1411.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,11)
	cl1420.lirdistallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirdistallb.eps')

	pylab.subplot(4,4,12)
	cl1301.lirdistallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,13)
	cl1138.lirdistallsub()
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,14)
	cl1018.lirdistallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,15)
	cl1059.lirdistallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,16)
	cl1202.lirdistallsub()
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	
	pylab.text(-1.,-.5,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=32,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.4,2.,r'$\rm log_{10}(N_{gal})$',fontsize=32,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,xlabel,fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,ylabel,fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('lirdistallb.eps')


def plotlirdistcomb():
	pylab.cla()
	pylab.clf()
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=46.
	xlaby=-.3
	ylabx=43.2
	ylaby=2.

	xlabx=1.-1
	xlaby=-.2
	ylabx=-.25-1
	ylaby=1.
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)



	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	lir=[]
	morph=[]
	morphflag=[]
	specflag=[]


	#morph=morph+list(c1.morphmemb)+list(c2.morphmemb)+list(c3.morphmemb)+list(c4.morphmemb)+list(c5.morphmemb)+list(c6.morphmemb)+list(c7.morphmemb)+list(c8.morphmemb)+list(c9.morphmemb)+list(c10.morphmemb)+list(c11.morphmemb)+list(c12.morphmemb)+list(c13.morphmemb)+list(c14.morphmemb)+list(c15.morphmemb)+list(c16.morphmemb)+list(c17.morphmemb)


	#morphflag=morphflag+list(c1.morphmembflag)+list(c2.morphmembflag)+list(c3.morphmembflag)+list(c4.morphmembflag)+list(c5.morphmembflag)+list(c6.morphmembflag)+list(c7.morphmembflag)+list(c8.morphmembflag)+list(c9.morphmembflag)+list(c10.morphmembflag)+list(c11.morphmembflag)+list(c12.morphmembflag)+list(c13.morphmembflag)+list(c14.morphmembflag)+list(c15.morphmembflag)+list(c16.morphmembflag)+list(c17.morphmembflag)


	#lir=lir+list(c1.lirbell)+list(c2.lirbell)+list(c3.lirbell)+list(c4.lirbell)+list(c5.lirbell)+list(c6.lirbell)+list(c7.lirbell)+list(c8.lirbell)+list(c9.lirbell)+list(c10.lirbell)+list(c11.lirbell)+list(c12.lirbell)+list(c13.lirbell)+list(c14.lirbell)+list(c15.lirbell)+list(c16.lirbell)+list(c17.lirbell)
	lir=lir+list(c1.lirbell)+list(c2.lirbell)+list(c3.lirbell)+list(c4.lirbell)+list(c5.lirbell)+list(c6.lirbell)+list(c7.lirbell)+list(c8.lirbell)+list(c9.lirbell)+list(c10.lirbell)+list(c11.lirbell)+list(c12.lirbell)+list(c13.lirbell)+list(c14.lirbell)+list(c15.lirbell)+list(c16.lirbell)

	#morph=morph+list(c1.morphbell)+list(c2.morphbell)+list(c3.morphbell)+list(c4.morphbell)+list(c5.morphbell)+list(c6.morphbell)+list(c7.morphbell)+list(c8.morphbell)+list(c9.morphbell)+list(c10.morphbell)+list(c11.morphbell)+list(c12.morphbell)+list(c13.morphbell)+list(c14.morphbell)+list(c15.morphbell)+list(c16.morphbell)+list(c17.morphbell)
	morph=morph+list(c1.morphbell)+list(c2.morphbell)+list(c3.morphbell)+list(c4.morphbell)+list(c5.morphbell)+list(c6.morphbell)+list(c7.morphbell)+list(c8.morphbell)+list(c9.morphbell)+list(c10.morphbell)+list(c11.morphbell)+list(c12.morphbell)+list(c13.morphbell)+list(c14.morphbell)+list(c15.morphbell)+list(c16.morphbell)

	specflag=specflag+list(c1.specflagbell)+list(c2.specflagbell)+list(c3.specflagbell)+list(c4.specflagbell)+list(c5.specflagbell)+list(c6.specflagbell)+list(c7.specflagbell)+list(c8.specflagbell)+list(c9.specflagbell)+list(c10.specflagbell)+list(c11.specflagbell)+list(c12.specflagbell)+list(c13.specflagbell)+list(c14.specflagbell)+list(c15.specflagbell)+list(c16.specflagbell)

	specz=[]
	specz=specz+list(c1.speczbell)+list(c2.speczbell)+list(c3.speczbell)+list(c4.speczbell)+list(c5.speczbell)+list(c6.speczbell)+list(c7.speczbell)+list(c8.speczbell)+list(c9.speczbell)+list(c10.speczbell)+list(c11.speczbell)+list(c12.speczbell)+list(c13.speczbell)+list(c14.speczbell)+list(c15.speczbell)+list(c16.speczbell)

	photz=[]
	photz=photz+list(c1.photzbell)+list(c2.photzbell)+list(c3.photzbell)+list(c4.photzbell)+list(c5.photzbell)+list(c6.photzbell)+list(c7.photzbell)+list(c8.photzbell)+list(c9.photzbell)+list(c10.photzbell)+list(c11.photzbell)+list(c12.photzbell)+list(c13.photzbell)+list(c14.photzbell)+list(c15.photzbell)+list(c16.photzbell)

	#print "here's lir",len(lir),len(morph)
	#for i in range(20):
	#	print i,'Lir morph =',lir[i],morph[i]
	#print "here's morph", morph


	lir=N.array(lir,'d')#
	morph=N.array(morph,'f')
	specflag = N.array(specflag,'f')

	#elir=[]#expand morphologies from 3 gems classifiers
	#emorph=[]
	#espec=[]
	#for i in range(len(lirhiz)):
	#	       m=morphhiz[i]
	#	       if pylab.average(m) < 0.:
	#		       continue
	#	       for j in range(len(m)):
	#		       elir.append(lirhiz[i])
	#		       emorph.append(m[j])
	#		       espec.append(specflag[i])
	#elir=N.array(elir,'d')
	#emorph=N.array(emorph,'f')
	#especflag=N.array(especflag,'f')
	#lir=elir
	#morph=emorph
	#specflag=especflag
	#for i in range(10):
	#	print i,lir[i],morph[i]

	print "number of galaxies for comparison w/bell, before lir cut = ",len(lir)
	print "the average redshift of these sources is = ",N.average(photz)
	#print "here's morph array",morph
	elir=N.compress((morph < 0.1),lir)
	slir=N.compress(abs(morph -1.) < .1,lir)
	ilir=N.compress(abs(morph -2.) < .1,lir)
	plir=N.compress(abs(morph -3.) < .1,lir)
	tot=N.sum(lir)
	ntot=len(lir)/3.#we have 3 classifiers in ediscs

	elirspec=N.compress((morph < 0.1) & (specflag > .1),lir)
	slirspec=N.compress((abs(morph -1.) < .1) & (specflag > .1),lir)
	ilirspec=N.compress((abs(morph -2.) < .1) & (specflag > .1),lir)
	plirspec=N.compress((abs(morph -3.) < .1) & (specflag > .1),lir)
	lirspec=N.compress((specflag > .1),lir)
	spectot=N.sum(lirspec)
	ntotspec=len(lirspec)/3.#we have 3 classifiers in ediscs
	print "for spec and phot members:"
	print "Fraction of Lir from Spirals = %5.3f"%(N.sum(slir)/tot)
	print "Fraction of Lir from E/S0s = %5.3f"%(N.sum(elir)/tot)
	print "Fraction of Lir from Irr = %5.3f"%(N.sum(ilir)/tot)
	print "Fraction of Lir from Pec/Comp = %5.3f"%(N.sum(plir)/tot)

	print "Fraction of gal from Spirals = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(slir)/3.,ntot))
	print "Fraction of gal from E/S0s = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(elir)/3.,ntot))
	print "Fraction of gal from Irr = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(ilir)/3.,ntot))
	print "Fraction of gal from Pec/Comp = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(plir)/3.,ntot))


	(a,b)=stats.stats.ks_2samp(lirspec,lir)
	print 'KS - test comparing lir for spec and all members:  D value = %8.4f, p-value = %8.4f'%(a,b)
	print "for spec members only: # of spec members = ", ntotspec
	print "Fraction of Lir from Spirals = %5.3f"%(N.sum(slirspec)/spectot)
	print "Fraction of Lir from E/S0s = %5.3f"%(N.sum(elirspec)/spectot)
	print "Fraction of Lir from Irr = %5.3f"%(N.sum(ilirspec)/spectot)
	print "Fraction of Lir from Pec/Comp = %5.3f"%(N.sum(plirspec)/spectot)

	print "Fraction of gal from Spirals = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(slirspec)/3.,ntotspec))
	print "Fraction of gal from E/S0s = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(elirspec)/3.,ntotspec))
	print "Fraction of gal from Irr = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(ilirspec)/3.,ntotspec))
	print "Fraction of gal from Pec/Comp = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(plirspec)/3.,ntotspec))

	#print "here's elir array",elir
	#print "here's slir array",slir
	#print "here's ilir array",ilir



#xbins=N.array([-7,-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11,66,111],'f')
#pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', '','Sb','','Sc','','Sd','','Sm','Im','Irr','?','N/A') ,fontsize=18)
	pylab.subplot(221)
	lirdistcombsub(lir,elir,'All','E/S0')
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	xaxissfrlog(ax)

	pylab.subplot(222)
	lirdistcombsub(lir,slir,'All','Sa-Sm')
 	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	xaxissfrlog(ax)

	pylab.subplot(223)
	lirdistcombsub(lir,ilir,'All','Irr/Comp')
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(224)
	lirdistcombsub(lir,plir,'All','Pec/Int')
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	#pylab.subplot(224)
	#ax=pylab.gca()
	##ax.set_xticklabels(([]))
	##ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.text(xlabx,xlaby,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(xlabx,2.1,r'$\rm SFR \ (M_{\odot}/yr)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(ylabx,ylaby,r'$\rm log_{10}(N_{gal}/log_{10}L_{IR})$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.savefig('lirdistcomb.eps')



	pylab.cla()
	pylab.clf()
	pylab.subplot(221)
	lirdistcombsub(lirspec,elirspec,'All','E/S0')
	ax=pylab.gca()
	ax.set_xticklabels(([]))

	pylab.subplot(222)
	lirdistcombsub(lirspec,slirspec,'All','Sa-Sm')
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(223)
	lirdistcombsub(lirspec,ilirspec,'All','Irr/Comp')
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(224)
	lirdistcombsub(lirspec,plirspec,'All','Pec/Int')
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	#pylab.subplot(224)
	#ax=pylab.gca()
	##ax.set_xticklabels(([]))
	##ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.text(xlabx,xlaby,r'$\rm log_{10}(L_{IR}/L\odot)$',fontsize=28,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(xlabx,2.2,r'$\rm SFR \ (M_{\odot}/yr)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(ylabx,ylaby,r'$\rm log_{10}(N_{gal}/log_{10}L_{IR})$',fontsize=28,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.savefig('lirdistcombspec.eps')

	
	#compare distribution of LIR between GEMS and our sample
	pylab.cla()
	pylab.clf()
	lirmin=10.75#for lowz gems sample
	lmin=10.**lirmin
	lirloz=[]
	lirloz=lirloz+list(c7.lirmemb)+list(c8.lirmemb)+list(c9.lirmemb)+list(c10.lirmemb)+list(c11.lirmemb)+list(c12.lirmemb)+list(c13.lirmemb)+list(c14.lirmemb)+list(c15.lirmemb)+list(c16.lirmemb)#all cluster members in high-z clusters
	lirlz=N.array(lirloz,'d')
	edilirlz=N.compress(lirlz > lmin,lirlz)

	lzgemslir=N.compress((gems.lzLir > lmin) & (gems.lzwhmips <1.1),gems.lzLir)
	lzgemsvtype=gems.lzgvtype[pylab.where((gems.lzLir > lmin)& (gems.lzwhmips <1.1))]
	print "Number of GEMS and EDisCS galaxies for loz Lir comparison, after lir cur =  ",len(lzgemslir)/3.,len(edilirlz)
	print "Average Lir of GEMS (Lir > lirmin) = ",pylab.average(lzgemslir),pylab.log10(pylab.average(lzgemslir))



	lirdistcombsub2(lzgemslir,edilirlz,'GEMS','EDisCS')
	pylab.axis([10.3,12.7,0.,3.5])
	pylab.xlabel(r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=36)
	#pylab.ylabel(r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L/$\rm \Sigma n_{gal}$)',fontsize=24)
	#pylab.ylabel(r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L)',fontsize=24)
	pylab.ylabel(r'$\rm log_{10}(N_{gal}/log_{10}L_{IR})$',fontsize=36)
	ax=pylab.gca()
	pylab.text(.1,.9,r'$\rm 0.42 < z < 0.6, log_{10}(L_{IR}/L_\odot) > 10.75$',fontsize=28,color='k',horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.5,1.1,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,color='k',horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(.15,.8,r'$\rm GEMS$',fontsize=28,color='0.5',horizontalalignment='left',transform=ax.transAxes)

	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)
	xaxissfrlog(ax)
	pylab.savefig('GemsEdiscsLirComplz.eps')
	(a,b)=stats.ks_2samp(lzgemslir,edilirlz)
	print '2-sample KS test comparing Lir dist of z< 0.6 cluster and field D-value = %6.4f p-value=%6.4f'%(a,b)



	pylab.cla()
	pylab.clf()
	lirmin=10.95#for highz gems sample
	lmin=10.**lirmin
	lirhiz=[]
	lirhiz=lirhiz+list(c1.lirmemb)+list(c2.lirmemb)+list(c3.lirmemb)+list(c4.lirmemb)+list(c5.lirmemb)+list(c6.lirmemb)#+list(c7.lirmemb)#all cluster members in high-z clusters
	lirhiz=N.array(lirhiz,'d')
	edilir=N.compress(lirhiz > lmin,lirhiz)

	gemslir=N.compress((gems.hzLir > lmin) & (gems.hzwhmips <1.1),gems.hzLir)
	gemsvtype=gems.hzgvtype[pylab.where((gems.hzLir > lmin)& (gems.hzwhmips <1.1))]
	print "Number of GEMS and EDisCS galaxies for Lir comparison, after lir cur =  ",len(gemslir)/3.,len(edilir)
	print "Average Lir of GEMS (Lir > lirmin) = ",pylab.average(gemslir),pylab.log10(pylab.average(gemslir))



	lirdistcombsub2(gemslir,edilir,'GEMS','EDisCS')
	pylab.axis([10.3,12.7,0.,3.5])
	pylab.xlabel(r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=36)
	#pylab.ylabel(r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L/$\rm \Sigma n_{gal}$)',fontsize=24)
	#pylab.ylabel(r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L)',fontsize=24)
	pylab.ylabel(r'$\rm log_{10}(N_{gal}/log_{10}L_{IR})$',fontsize=36)
	ax=pylab.gca()
	pylab.text(.1,.9,r'$\rm 0.6 < z < 0.8, log_{10}(L_{IR}/L_\odot) > 10.95$',fontsize=28,color='k',horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.5,1.1,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,color='k',horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(.15,.8,r'$\rm GEMS$',fontsize=28,color='0.5',horizontalalignment='left',transform=ax.transAxes)

	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)
	xaxissfrlog(ax)

	(a,b)=stats.ks_2samp(gemslir,edilir)
	print '2-sample KS test comparing Lir dist of z> 0.6 cluster and field D-value = %6.4f p-value=%6.4f'%(a,b)

	pylab.savefig('GemsEdiscsLirComp.eps')




	pylab.cla()
	pylab.clf()
	lirmin=10.95
	lmin=10.**lirmin
	lirhiz=[]
	lirhiz=lirhiz+list(c1.lirbell)+list(c2.lirbell)+list(c3.lirbell)+list(c4.lirbell)+list(c5.lirbell)+list(c6.lirbell)#+list(c7.lirbell)#all cluster members in high-z clusters that also have HST morphologies
	morphhiz=[]
	morphhiz=morphhiz+list(c1.morphbell)+list(c2.morphbell)+list(c3.morphbell)+list(c4.morphbell)+list(c5.morphbell)+list(c6.morphbell)#+list(c7.morphbell)

	#print lirhiz
	lirhiz=pylab.array(lirhiz,'d')
	morphhiz=N.array(morphhiz,'f')

	#elir=[]#expand morphologies from 3 gems classifiers
	#emorph=[]
	#for i in range(len(lirhiz)):
	#	       m=morphhiz[i]
	#	       if pylab.average(m) < 0.:
	#		       continue
	#	       for j in range(len(m)):
	#		       elir.append(lirhiz[i])
	#		       emorph.append(m[j])

	#elir=pylab.array(elir,'d')
	#emorph=pylab.array(emorph,'f')

	
	edilir=N.compress(lirhiz > lmin,lirhiz)
	edimorph=N.compress(lirhiz > lmin,morphhiz)
	print 'length of edilir, edimorph = ',len(edilir),len(edimorph)
	for i in range(10):
		print i,edilir[i],edimorph[i]
	elir=N.compress((edimorph < 0.1),edilir)
	slir=N.compress(abs(edimorph -1.) < .1,edilir)
	ilir=N.compress(abs(edimorph -2.) < .1,edilir)
	plir=N.compress(abs(edimorph -3.) < .1,edilir)
	tot=N.sum(lir)



	#elir=N.compress(elir > lmin,elir)
	#slir=N.compress(slir > lmin,slir)
	#ilir=N.compress(ilir > lmin,ilir)
	#plir=N.compress(plir > lmin,plir)
	tot=N.sum(elir)+N.sum(slir)+N.sum(ilir)+N.sum(plir)
	ntot=len(elir)+len(slir)+len(ilir)+len(plir)
	ntot=ntot/3.
	print "EDisCS After lirmin cut and HST morph:"
	print "Fraction of Lir from Spirals = %5.3f"%(N.sum(slir)/tot)
	print "Fraction of Lir from E/S0s = %5.3f"%(N.sum(elir)/tot)
	print "Fraction of Lir from Irr = %5.3f"%(N.sum(ilir)/tot)
	print "Fraction of Lir from Pec/Comp = %5.3f"%(N.sum(plir)/tot)
	print "Fraction of gal from Spirals = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(slir)/3.,ntot))
	print "Fraction of gal from E/S0s = %5.3f %5.3f+ %5.3f"%(my.ratioerror(len(elir)/3.,ntot))
	print "Fraction of gal from Irr = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(ilir)/3.,ntot))
	print "Fraction of gal from Pec/Comp = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(plir)/3.,ntot))
	print "total number of galaxies w/morph and lir>lirmin = ",ntot
	#print "here's elir array",elir
	#print "here's slir array",slir
	#print "here's ilir array",ilir


#0=E,1=S0,2=Sa,3=Sbc,4=Sdm,5=Irr,6=Pec/Int,7=comp,X=no image
	print "len of gemsvtype and gemslir = ",len(gemsvtype),len(gemslir)
	glir=[]#expand morphologies from 3 gems classifiers
	gmorph=[]
	for i in range(len(gemslir)):
		       m=gemsvtype[i]
		       #print 'm = ',m
		       if pylab.average(m) < 0.:
			       continue
		       for j in range(len(m)):
			       glir.append(gemslir[i])
			       gmorph.append(m[j])

	glir=pylab.array(glir,'d')
	gmorph=pylab.array(gmorph,'f')
	gelir=N.compress(gmorph < 1.9,glir)
	flag=(gmorph > 1.9) & (gmorph < 4.1)
	gslir=N.compress(flag,glir)
	gilir=N.compress((abs(gmorph -5.) < 0.1) | (abs(gmorph -7.) < 0.1),glir)
	gplir=N.compress((abs(gmorph -6.) < 0.1),glir)
	gtot=N.sum(gelir)+N.sum(gslir)+N.sum(gilir)+N.sum(gplir)
	ngtot=len(gelir)+len(gslir)+len(gilir)+len(gplir)
	ngtot=ngtot/3.

	print "GEMS After lirmin cut:"
	print "Fraction of Lir from Spirals = %5.3f"%(N.sum(gslir)/gtot)
	print "Fraction of Lir from E/S0s = %5.3f"%(N.sum(gelir)/gtot)
	print "Fraction of Lir from Irr = %5.3f"%(N.sum(gilir)/gtot)
	print "Fraction of Lir from Pec/Comp = %5.3f"%(N.sum(gplir)/gtot)

	print "Fraction of gal from Spirals = %5.3f %5.3f + %5.3f"%(my.ratioerror(len(gslir)/3.,ngtot))
	print "Fraction of gal from E/S0s = %5.3f %5.3f +%5.3f"%(my.ratioerror(len(gelir)/3.,ngtot))
	print "Fraction of gal from Irr = %5.3f %5.3f +%5.3f"%(my.ratioerror(len(gilir)/3.,ngtot))
	print "Fraction of gal from Pec/Comp = %5.3f %5.3f +%5.3f"%(my.ratioerror(len(gplir)/3.,ngtot))

	#print "Fraction of galaxies from Spirals = %5.3f"%(len(gslir)/ngtot)
	#print "Fraction of galaxies from E/S0s = %5.3f"%(len(gelir)/ngtot)
	#print "Fraction of galaxies from Irr = %5.3f"%(len(gilir)/ngtot)
	#print "Fraction of galaxies from Pec/Comp = %5.3f"%(len(gplir)/ngtot)

	print "total number of galaxies w/morph and lir>lirmin = ",ngtot

	scale=1.*ngtot/(1.*ntot)#ngal(gems)/ngal(ediscs)

#xbins=N.array([-7,-6,-5,-2,1,2,3,4,5,6,7,8,9,10,11,66,111],'f')
#pylab.xticks(ind+width, ('St', 'C', 'E', 'S0', 'Sa', '','Sb','','Sc','','Sd','','Sm','Im','Irr','?','N/A') ,fontsize=18)
	pylab.subplot(221)
	lirdistcombsub3(gelir,elir,'GEMS E/S0','E/S0',scale)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.text(.7,.85,'GEMS',fontsize=18,color='0.4',transform=ax.transAxes)
	pylab.text(.7,.75,'EDisCS',fontsize=18,color='k',transform=ax.transAxes)

	xaxissfrlog(ax)

	pylab.subplot(222)
	lirdistcombsub3(gslir,slir,'GEMS Sa-Sm','Sa-Sm',scale)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	xaxissfrlog(ax)

	pylab.subplot(223)
	lirdistcombsub3(gilir,ilir,'GEMS Irr/Comp','Irr/Comp',scale)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(224)
	lirdistcombsub3(gplir,plir,'GEMS Pec/Int','Pec/Int',scale)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	#pylab.subplot(224)
	#ax=pylab.gca()
	##ax.set_xticklabels(([]))
	##ax=pylab.gca()
	#ax.set_yticklabels(([]))



	pylab.text(xlabx,xlaby,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(xlabx,2.1,r'$\rm SFR \ (M_\odot/yr)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(ylabx,ylaby,r'$\rm log_{10}(N_{gal}/log_{10}L_{IR})$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.savefig('GEMSlirdistcomb.eps')




def lirdistcombsub(x,x2,label1,label2):#label 1 = All for morph plot
	greyscale='0.5'
	xmin=10.25
	xmax=13.
	nbin=10
	x=N.log10(x)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	dbin=xbin[1]-xbin[0]
	ybinerrm=pylab.zeros(len(ybinerr),'d')
	ybinerrp=pylab.zeros(len(ybinerr),'d')
	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	ybin=ybin/3.#account for classifications by 3 authors
	yberrm=yberrm/3.#account for classifications by 3 authors
	yberrp=yberrp/3.#account for classifications by 3 authors

	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			else:
				ybinerrm[i]=999.
			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		else:
			ybin[i]=-10.
	#print "here's ybin ",ybin
	#print "dbin = ",dbin
	#ybin=ybin/(dbin)
	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	pylab.plot(xp,yp,c=greyscale,ls='-')
	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm
	yerrall[1]=ybinerrp

	#print ybinerrm
	#print ybinerrp

	pylab.errorbar(xbin+dbin/2.,ybin,yerr=yerrall,fmt=None,ecolor='0.5')


	x=N.log10(x2)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybinerrm=pylab.zeros(len(ybinerr),'d')
	ybinerrp=pylab.zeros(len(ybinerr),'d')
	ybin=ybin/3.#account for classifications by 3 authors
	yberrm=yberrm/3.#account for classifications by 3 authors
	yberrp=yberrp/3.#account for classifications by 3 authors

	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			else:
				ybinerrm[i]=999.

			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		else:
			ybin[i]=-10.

#	for i in range(len(ybin)):
#		if ybin[i] > 0.:#

#			ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((ybin[i]-ybinerr[i])/dbin)
#			ybinerrp[i]=pylab.log10((ybin[i]+ybinerr[i])/dbin)-pylab.log10(ybin[i]/dbin)
#			if (ybin[i] - ybinerr[i]) < 0.1:
#				ybinerrm[i]=pylab.log10(ybin[i]/dbin)
#			ybin[i]=N.log10(ybin[i]/dbin)

#		else:
#			ybin[i]=-.1
	#print "here's ybin ",ybin
	my.drawhistpylab(xbin,(ybin),'k-')

	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm
	yerrall[1]=ybinerrp
	    
	pylab.errorbar(xbin+dbin/2.,ybin,yerr=yerrall,fmt=None,ecolor='k')

	pylab.axis([9.8,12.4,0.,3.4])
	ax=pylab.gca()
	#s='All'
	s=label1
	pylab.text(.7,.9,s,fontsize=18,color=greyscale,transform=ax.transAxes)
	pylab.text(.7,.8,label2,fontsize=18,color='k',transform=ax.transAxes)

	#ax=pylab.gca()
	#ax.set_yscale('log')

def lirdistcombsub2(x,x2,label1,label2):#label 1 = All for morph plot
	greyscale='0.5'
	xmin=10.5
	xmax=13.
	nbin=15
	#run ks test comparing two populations
	ngem=len(x)
	nedi=len(x2)
	scale=1.*ngem/3./(1.*nedi)
	x=N.log10(x)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)#bin gems counts

	dbin=xbin[1]-xbin[0]
	ybinerrm=pylab.zeros(len(ybinerr),'d')
	ybinerrp=pylab.zeros(len(ybinerr),'d')
	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	ybin=ybin/3.#correct for 3 gems classifiers 
	yberrm=yberrm/3.#correct for 3 gems classifiers 
	yberrp=yberrp/3.#correct for 3 gems classifiers 

	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				if (yberrm[i]/dbin) < 1.:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
					print 'got here ybinerrm = ',ybinerrm[i]
				else:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)


			else:
				ybinerrm[i]=999.
			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		elif ybin[i] < .01:
			ybin[i]=-10.
			ybinerr[i]=0.
			ybinerrm[i]=0.
			ybinerrp[i]=0.

	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	#ybin=ybin/3.#correct for 3 edi classifiers 
	#ybinerr=ybinerr/3.#correct for 3 edi classifiers 

	pylab.plot(xp,yp,c=greyscale,ls='-',label=label1)
	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm
	yerrall[1]=ybinerrp
	#yerrall[0]=yberrm
	#yerrall[1]=yberrp
	pylab.errorbar(xbin+dbin/2.,ybin,yerr=yerrall,fmt=None,ecolor='0.5',label='_nolegend_')
	
	#ediscs data
	x=N.log10(x2)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybinerrm=pylab.zeros(len(ybinerr),'f')
	ybinerrp=pylab.zeros(len(ybinerr),'f')
	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				if (yberrm[i]/dbin) < 1.:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)+pylab.log10((yberrm[i])/dbin)
					print 'got here ybinerrm = ',ybinerrm[i]
				else:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)

			else:
				ybinerrm[i]=999.

			#ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		else:
			ybin[i]=-10.
			ybinerr[i]=0.
	#print "here's ybin ",ybin
	#for i in range(len(xbin)):
	#	print i, 'EDisCS',xbin[i],ybin[i],ybinerrm[i],ybinerrp[i]

	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	#pylab.plot(xp,yp,c='k',ls='-')
	pylab.plot(xp,yp+pylab.log10(scale),c='k',ls='-',label=label2)
	#pylab.plot(xp,yp*scale,c='k',ls='--')
	#yerrall=zip(ybinerrm,ybinerrp)
	#yerrall=pylab.array(yerrall)
	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm#+pylab.log10(scale)
	yerrall[1]=ybinerrp#+pylab.log10(scale)

	pylab.errorbar(xbin+dbin/2.,ybin+pylab.log10(scale),yerr=yerrall,fmt=None,ecolor='k',label='_nolengend_')
	#my.drawhistpylab(xbin,(ybin),'k-')

	#pylab.axis([9.8,12.4,0.,3.])
	ax=pylab.gca()
	#s='All'
	s=label1
	#pylab.text(.75,.8,s,fontsize=20,color=greyscale,transform=ax.transAxes)
	#pylab.text(.75,.75,label2,fontsize=20,color='k',transform=ax.transAxes)

	location=.7,.7
	pylab.legend(loc=location)


	#ax=pylab.gca()
	#ax.set_yscale('log')

def lirdistcombsub3(x,x2,label1,label2,scale):#label 1 = All for morph plot
	greyscale='0.5'
	xmin=10.9
	xmax=13.
	nbin=10
	x=N.log10(x)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	#ybin=ybin/3.#correct for 3 gems classifiers 
	#ybinerr=ybinerr/3.#correct for 3 gems classifiers 
	dbin=xbin[1]-xbin[0]
	ybinerrm=pylab.zeros(len(ybinerr),'d')
	ybinerrp=pylab.zeros(len(ybinerr),'d')
	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	ybin=ybin/3.#correct for 3 gems classifiers 
	yberrm=yberrm/3.#correct for 3 gems classifiers 
	yberrp=yberrp/3.#correct for 3 gems classifiers 

	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				if (yberrm[i]/dbin) < 1.:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)+pylab.log10((yberrm[i])/dbin)
					print 'got here ybinerrm = ',ybinerrm[i]
				else:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			else:
				ybinerrm[i]=999.
				print 'Error, yberrm < 0.',ybin[i],dbin,yberrm[i]

			#ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		else:
			ybin[i]=-10.
			ybinerr[i]=0.

	#print "here's ybin ",ybin
	#print "dbin = ",dbin
	#ybin=ybin/(dbin)
	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	pylab.plot(xp,yp,c=greyscale,ls='-')
	#print 'HEEYYYYYYYY'
	#print ybin
	#print ybinerrm
	#print ybinerrp
	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm
	yerrall[1]=ybinerrp
	    
	pylab.errorbar(xbin+dbin/2.,ybin,yerr=yerrall,fmt=None,ecolor='0.5')


	x=N.log10(x2)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybin=ybin/3.#correct for 3 gems classifiers 
	ybinerr=ybinerr/3.#correct for 3 gems classifiers 

	(yberrm,yberrp)=poisson.Poisson_lim(ybin)
	ybin=ybin*scale#scale by ngal(gems)/ngal(ediscs)
	ybinerr=ybinerr*scale#scale by ngal(gems)/ngal(ediscs)
	ybinerrm=pylab.zeros(len(ybinerr),'d')
	ybinerrp=pylab.zeros(len(ybinerr),'d')
	yberrm=yberrm*scale
	yberrp=yberrp*scale

	for i in range(len(ybin)):
		if ybin[i] > 0.:
			if yberrm[i] > 0.:
				if (yberrm[i]/dbin) < 1.:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)+pylab.log10((yberrm[i])/dbin)
					print 'got here ybinerrm = ',ybinerrm[i]
				else:
					ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)


			else:
				ybinerrm[i]=999.
			#ybinerrm[i]=pylab.log10(ybin[i]/dbin)-pylab.log10((yberrm[i])/dbin)
			ybinerrp[i]=pylab.log10((yberrp[i])/dbin)-pylab.log10(ybin[i]/dbin)
			ybin[i]=N.log10(ybin[i]/dbin)
		else:
			ybin[i]=-10.
			ybinerr[i]=0.

	#print "here's ybin ",ybin
	my.drawhistpylab(xbin,(ybin),'k-')

	yerrall=pylab.zeros([2,len(ybinerrm)],'f')
	yerrall[0]=ybinerrm
	yerrall[1]=ybinerrp
	    
	pylab.errorbar(xbin+dbin/2.,ybin,yerr=yerrall,fmt=None,ecolor='k')

	#pylab.axis([9.8,12.4,0.,3.])
	pylab.axis([10.3,12.9,0.,2.9])
	ax=pylab.gca()
	#s='All'
	s=label1
	#pylab.text(.1,.85,s,fontsize=18,color=greyscale,transform=ax.transAxes)
	pylab.text(.08,.85,label2,fontsize=18,color='k',transform=ax.transAxes)

	#ax=pylab.gca()
	#ax.set_yscale('log')


def plotlirdistz():
	pylab.cla()
	pylab.clf()
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)
	matplotlib.rc('legend',numpoints=3,fontsize=20,markerscale=1)

	xlabx=46.
	xlaby=-.3
	ylabx=43.2
	ylaby=2.

	pylab.cla()
	pylab.clf()

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	lir=[]
	lirloz=[]
	lirhiz=[]
	drloz=[]
	drhiz=[]

	lir=lir+list(c1.lirmemb)+list(c2.lirmemb)+list(c3.lirmemb)+list(c4.lirmemb)+list(c5.lirmemb)+list(c6.lirmemb)+list(c7.lirmemb)+list(c8.lirmemb)+list(c9.lirmemb)+list(c10.lirmemb)+list(c11.lirmemb)+list(c12.lirmemb)+list(c13.lirmemb)+list(c14.lirmemb)+list(c15.lirmemb)+list(c16.lirmemb)#+list(c17.lirmemb)

	lirloz=lirloz+list(c8.lirmemb)+list(c9.lirmemb)+list(c10.lirmemb)+list(c11.lirmemb)+list(c12.lirmemb)+list(c13.lirmemb)+list(c14.lirmemb)+list(c15.lirmemb)+list(c16.lirmemb)#+list(c17.lirmemb)

	lirhiz=lirhiz+list(c1.lirmemb)+list(c2.lirmemb)+list(c3.lirmemb)+list(c4.lirmemb)+list(c5.lirmemb)+list(c6.lirmemb)+list(c7.lirmemb)

	drloz=drloz+list(c8.drmemb)+list(c9.drmemb)+list(c10.drmemb)+list(c11.drmemb)+list(c12.drmemb)+list(c13.drmemb)+list(c14.drmemb)+list(c15.drmemb)+list(c16.drmemb)#+list(c17.drmemb)

	drhiz=drhiz+list(c1.drmemb)+list(c2.drmemb)+list(c3.drmemb)+list(c4.drmemb)+list(c5.drmemb)+list(c6.drmemb)+list(c7.drmemb)
	nclloz=10.
	nclhiz=7.

	lir=N.array(lir,'d')#
	lirloz=N.array(lirloz,'d')#
	lirhiz=N.array(lirhiz,'d')#
	drloz=N.array(drloz,'d')#
	drhiz=N.array(drhiz,'d')#

	lirloz=N.compress((drloz < 1.) & (lirloz > 10.**lirmin/3.*5.),lirloz)
	lirhiz=N.compress((drhiz < 1.) & (lirhiz > 10.**lirmin/3.*5.),lirhiz)

	x=lir
	x2=lirloz
	x3=lirhiz
	greyscale='0.3'
	xmin=9
	xmax=13.
	nbin=40
	x=N.log10(x)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	#errybin=N.sqrt(ybin)/len(x)
	#ybinerr=ybinerr/len(x)
	#ybin=ybin/len(x)#normalize by number of galaxies
	ybinerr=ybinerr/nclloz#len(x)
	ybin=ybin/nclloz#len(x)#normalize by number of galaxies
	dbin=xbin[1]-xbin[0]
	#for i in range(len(ybin)):
	#	if ybin[i] > 0.:
	#		ybin[i]=N.log10(ybin[i])
	#	else:
	#		ybin[i]=-99.
	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	#pylab.plot(xp,yp,c=greyscale,ls='-',label='All')
	#pylab.errorbar(xbin,ybin,ybinerr,fmt='ko')

	x=N.log10(x2)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybinerr=ybinerr/len(x)
	ybin=ybin/len(x)
	#for i in range(len(ybin)):
	#	if ybin[i] > 0.:
	#		ybin[i]=N.log10(ybin[i])
	#	else:
	#		ybin[i]=-99.
	#print "here's ybin ",ybin

	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'k-')
	pylab.plot(xp,yp,c='k',ls='-',label=r'$z < 0.6$')
	pylab.errorbar(xbin+dbin/2.,ybin,ybinerr,fmt=None,ecolor='k',label='_nolegend_')

	x=N.log10(x3)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)
	ybinerr=ybinerr/len(x)
	ybin=ybin/len(x)
	#for i in range(len(ybin)):
	#	if ybin[i] > 0.:
	#		ybin[i]=N.log10(ybin[i])
	#	else:
	#		ybin[i]=-99.
	#print "here's ybin ",ybin

	(xp,yp)=my.drawhistpylabnoplot(xbin,(ybin),'b-')
	pylab.plot(xp,yp,c='0.5',ls='-',label=r'$z > 0.6$')
	pylab.errorbar(xbin+dbin/2.,ybin,ybinerr,fmt=None,ecolor='0.5',label='_nolegend_')
	#my.drawhistpylab(xbin,(ybin),'b-',label=r'$z > 0.6$')

	pylab.legend(loc='best')
	#pylab.axis([9.8,12.6,-2.6,.2])
	pylab.axvline(x=N.log10(10.**lirmin*5./3),color='k',ls='--')
	#pylab.plot(xp-.4,yp,c='0.7',ls='-',label=r'$z > 0.6$')
	pylab.axis([10.6,12.6,-.01,.45])
	ax=pylab.gca()
	s='All'
	#pylab.text(.7,.9,s,fontsize=18,color=greyscale,transform=ax.transAxes)
	#pylab.text(.7,.8,label2,fontsize=18,color='k',transform=ax.transAxes)




	xlabx=0.5
	xlaby=-.1
	ylabx=-.1
	ylaby=0.5

	pylab.text(xlabx,xlaby,r'log$_{10}$(L$\rm _{IR}$/L$_\odot$)',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(ylabx,ylaby,r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L)',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.text(ylabx,ylaby,r'Fraction of Galaxies per bin',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.savefig('lirdistz.eps')


def calcmatchstats():
	nmatch=0.
	nomatch=0.
	nmultimatch=0.

	if mode < 0.1:
		clusters=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16]
		outfile100=open('matchstats.dat','w')
		for cl in clusters:
			s="%i %i %i\n"%(cl.nedi24match,cl.nedi24nomatch,cl.nedi24multimatch)
			outfile100.write(s)
		outfile100.close()
		for cl in clusters:
			nmatch=cl.nedi24match+nmatch
			nomatch=nomatch+cl.nedi24nomatch
			nmultimatch=nmultimatch+cl.nedi24multimatch

	if mode > 0.1:
		infile100=open('matchstats.dat','r')
		for line in infile100:
			t=line.split()
			nmatch=nmatch+float(t[0])
			nomatch=nomatch+float(t[1])
			nmultimatch=nmultimatch+float(t[2])
		infile100.close()
	ntot=nmatch+nomatch
	(a,bu,bd)=my.ratioerror(nmatch,ntot)
	print "Number of 24um sources w/edi match: %i/%i = %5.2f +%5.2f -%5.2f\n"%(nmatch,ntot,a,bu,bd)
	(a,bu,bd)=my.ratioerror(nomatch,ntot)
	print "Number of 24um sources w/no edi match: %i/%i = %5.2f +%5.2f -%5.2f\n"%(nomatch,ntot,a,bu,bd)
	(a,bu,bd)=my.ratioerror(nmultimatch,ntot)
	print "Number of 24um sources w/multi edi matches: %i/%i = %5.2f +%5.2f -%5.2f\n"%(nmultimatch,ntot,a,bu,bd)

def IRoptdist(): #histogram of distance b/w IR and opt sources

	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	clusters=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16]
	for i in range(1,17):
		pylab.subplot(4,4,i)
		clusters[i-1].IRoptdist()
		if i< 13:
			ax=pylab.gca()
			ax.set_xticklabels(([]))
		f=1.*i/4.-pylab.floor(1.*i/4.)
		if abs(f-0.25) > .01:
			ax=pylab.gca()
			ax.set_yticklabels(([]))
	ax=pylab.gca()
	pylab.text(-1.,-.5,r'$\rm Distance from IR Galaxy (arcsec)$',fontsize=32,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.5,2.,r'$\rm N_{gal}/Square Arcsecond$',fontsize=32,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.xticks(fontsize=24)
	pylab.yticks(fontsize=24)

	pylab.savefig('IRoptdist.eps')
	
	

def calcIRfraco2():

	spec=[]
	spec=spec+list(c1.newspecmatchflag)+list(c2.newspecmatchflag)+list(c3.newspecmatchflag)+list(c4.newspecmatchflag)+list(c5.newspecmatchflag)+list(c6.newspecmatchflag)+list(c7.newspecmatchflag)+list(c8.newspecmatchflag)+list(c9.newspecmatchflag)+list(c10.newspecmatchflag)+list(c11.newspecmatchflag)+list(c12.newspecmatchflag)+list(c13.newspecmatchflag)+list(c14.newspecmatchflag)+list(c15.newspecmatchflag)+list(c16.newspecmatchflag)

	o2=[]
	o2=o2+list(c1.specEWOIIflag)+list(c2.specEWOIIflag)+list(c3.specEWOIIflag)+list(c4.specEWOIIflag)+list(c5.specEWOIIflag)+list(c6.specEWOIIflag)+list(c7.specEWOIIflag)+list(c8.specEWOIIflag)+list(c9.specEWOIIflag)+list(c10.specEWOIIflag)+list(c11.specEWOIIflag)+list(c12.specEWOIIflag)+list(c13.specEWOIIflag)+list(c14.specEWOIIflag)+list(c15.specEWOIIflag)+list(c16.specEWOIIflag)

	f24=[]
	f24=f24+list(c1.matchflag24)+list(c2.matchflag24)+list(c3.matchflag24)+list(c4.matchflag24)+list(c5.matchflag24)+list(c6.matchflag24)+list(c7.matchflag24)+list(c8.matchflag24)+list(c9.matchflag24)+list(c10.matchflag24)+list(c11.matchflag24)+list(c12.matchflag24)+list(c13.matchflag24)+list(c14.matchflag24)+list(c15.matchflag24)+list(c16.matchflag24)


	lir=[]
	lir=lir+list(c1.Lir)+list(c2.Lir)+list(c3.Lir)+list(c4.Lir)+list(c5.Lir)+list(c6.Lir)+list(c7.Lir)+list(c8.Lir)+list(c9.Lir)+list(c10.Lir)+list(c11.Lir)+list(c12.Lir)+list(c13.Lir)+list(c14.Lir)+list(c15.Lir)+list(c16.Lir)
	spec=N.array(spec,'d')
	o2=N.array(o2,'d')
	f24=N.array(f24,'d')
	lir=N.array(lir,'d')

	sumlirspec=N.compress((spec > 0.1) & (f24 > 0.1),lir)
	print "Number of spec members with 24um emission = ",len(sumlirspec)
	sumlirspec=N.sum(sumlirspec)
	sumlirspeco2=(N.compress((spec > 0.1) & (f24 > 0.1) & (o2 > 0.1),lir))
	print "Number of spec members with 24um emission = ",len(sumlirspeco2)
	sumlirspeco2=N.sum(sumlirspeco2)
	print "IR spec members: Fraction of tot IR from galaxies w/o2 emission = %5.2f"%(sumlirspeco2/sumlirspec)


def plotsigmaz():
	matplotlib.rc('legend',numpoints=1,fontsize=20,markerscale=1)
	ncl=16
	downsize=N.zeros(ncl,'f')
	z=N.zeros(ncl,'f')
	sigma=N.zeros(ncl,'f')
	errsigmap=N.zeros(ncl,'f')
	errsigmam=N.zeros(ncl,'f')
	nlir=N.zeros(ncl,'f')
	ntot=N.zeros(ncl,'f')
	nlirred=N.zeros(ncl,'f')
	ntotred=N.zeros(ncl,'f')
	nlirblue=N.zeros(ncl,'f')
	ntotblue=N.zeros(ncl,'f')

	sumlir=N.zeros(ncl,'f')
	sumlirerr=N.zeros(ncl,'f')
	drall=N.zeros(ncl,'f')

	nlirs=N.zeros(ncl,'f')
	ntots=N.zeros(ncl,'f')

	nlirsred=N.zeros(ncl,'f')
	ntotsred=N.zeros(ncl,'f')

	nlirsblue=N.zeros(ncl,'f')
	ntotsblue=N.zeros(ncl,'f')

	sumlirs=N.zeros(ncl,'f')
	sumlirerrs=N.zeros(ncl,'f')
	#lirmin=11.2
	lirmin=10.95
	lmin=10.**lirmin
	for i in range(ncl):
		if i == 0:
			cl=c1
		if i == 1:
			cl=c2
		if i == 2:
			cl=c3
		if i == 3:
			cl=c4
		if i == 4:
			cl=c5
		if i == 5:
			cl=c6
		if i == 6:
			cl=c7
		if i == 7:
			cl=c8
		if i == 8:
			cl=c9
		if i == 9:
			cl=c10
		if i == 10:
			cl=c11
		if i == 11:
			cl=c12
		if i == 12:
			cl=c13
		if i == 13:
			cl=c14
		if i == 14:
			cl=c15
		if i == 15:
			cl=c16
		if i == 16:
			cl=c17
		z[i]=(cl.zcl)
		sigma[i]=(cl.sigma)
		errsigmap[i]=cl.errsigmap
		errsigmam[i]=cl.errsigmam
		downsize[i]=cl.downsize
		n=0
		nred=0
		nblue=0
		nt=0
		ntred=0
		ntblue=0
		sum=0
		sumerr=0
		sumred=0
		sumblue=0
		nts=0
		ntsred=0
		ntsblue=0
		ns=0
		nsred=0
		nsblue=0
		sums=0
		sumsred=0
		sumsblue=0
		sumerrs=0
		for j in range(len(cl.Lir)):
			if cl.matchflag24[j] > -1.:#galaxy on 24 micron image
				#if cl.magI[j] < 23.:#galaxy brighter than spec cut
				if cl.MV[j] < cl.Mvcut:#galaxy brighter than spec cut
					if cl.dr[j] <= 1.:#galaxy w/in R200
						if cl.defmembflag[j] > 0:
							if cl.photmembflag[j] > 0:
								nt += 1
								if (cl.redflag[j]):
									ntred +=1
								else:
									ntblue += 1

								if cl.matchflag24[j] > 0:
									if cl.Lir[j] > lmin:
										if i == 15:
 											print j,cl.Lir[j]/1.e11,sum/1.e11
										sum += cl.Lir[j]
										sumerr += cl.errLir[j]**2
										n += 1
										if (cl.redflag[j]):
											nred +=1
										else:
											nblue += 1

						if cl.newspecmatchflag[j] > 0.:
							nts += 1
							if (cl.redflag[j]):
								ntsred +=1
							else:
								ntsblue += 1
							if cl.matchflag24[j] > 0:
								if cl.Lir[j] > lmin:
									sums += cl.Lir[j]
									sumerrs += cl.errLir[j]**2
									ns += 1

									if (cl.redflag[j]):
										nsred +=1
									else:
										nsblue += 1


		sumlir[i]=sum
		sumlirerr[i]=N.sqrt(sumerr)
		nlir[i]=n
		ntot[i]=nt
		nlirred[i]=nred
		ntotred[i]=ntred
		nlirblue[i]=nblue
		ntotblue[i]=ntblue
		sumlirs[i]=sums
		sumlirerrs[i]=N.sqrt(sumerrs)
		nlirs[i]=ns
		ntots[i]=nts
		nlirsred[i]=nsred
		ntotsred[i]=ntsred
		nlirsblue[i]=nsblue
		ntotsblue[i]=ntsblue
	#x=x+list(c1.z)+list(c2.z)+list(c3.z)+list(c4.z)+list(c5.z)+list(c6.z)+list(c7.z)+list(c8.z)+list(c9.z)+list(c10.z)+list(c11.z)+list(c12.z)+list(c13.z)+list(c14.z)+list(c15.z)+list(c16.z)+list(c17.z)

	#y=y+list(c1.sigma)+list(c2.sigma)+list(c3.sigma)+list(c4.sigma)+list(c5.sigma)+list(c6.sigma)+list(c7.sigma)+list(c8.sigma)+list(c9.sigma)+list(c10.sigma)+list(c11.sigma)+list(c12.sigma)+list(c13.sigma)+list(c14.sigma)+list(c15.sigma)+list(c16.sigma)+list(c17.sigma)
	output10=open('IntClusterProp.dat','w')
	output10.write('#z   sigma(km/s)  N_IR  N_tot  Sum(L_IR) (L_sun)  Err_sumlir')
	output11=open('IntClusterPropSpec.dat','w')
	output11.write('#z   sigma(km/s)  N_IR  N_tot  Sum(L_IR) (L_sun)  Err_sumlir')
	for i in range(len(ntot)):
		#print "sfrmass %i %6.4f sumlir=%5.2f nlir=%5.2f avelir=%5.2f ntot=%5.2f "%(i,z[i],sumlir[i]/1.e11,nlir[i],sumlir[i]/1.e11/nlir[i],ntot[i])
		#print "sfrmass spec sumlirs=%5.2f nlirs=%5.2f avelirs=%5.2f ntots=%5.2f "%(sumlirs[i]/1.e11,nlirs[i],sumlirs[i]/1.e11/nlirs[i],ntots[i])

		s="%6.4f %i %i %i %5.3e %5.3e\n"%(z[i],sigma[i],nlir[i],ntot[i],sumlir[i],sumlirerr[i])
		output10.write(s)
		s="%6.4f %i %i %i %5.3e %5.3e\n"%(z[i],sigma[i],nlirs[i],ntots[i],sumlirs[i],sumlirerrs[i])
		output11.write(s)
	output10.close()
	output11.close()
	pylab.cla()
	pylab.clf()
	#pylab.plot(z,sigma,'ko',markersize=16)
	sigerr=N.zeros([2,len(sigma)],'f')
	sigerr[0]=errsigmam
	sigerr[1]=errsigmam
	pylab.errorbar(z,sigma,yerr=sigerr,fmt='ko',markersize=16)
	pylab.axvline(x=0.6,color='k',ls='--')
	pylab.axis([0.36,0.84,120.,1250.])
	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)
	pylab.xlabel(r'$\rm z$',fontsize=36)
	pylab.ylabel(r'$\rm \sigma (km/s)$',fontsize=36)
	pylab.savefig('sigmaz.eps')


	pylab.cla()
	pylab.clf()
	pylab.plot(z,ntot,'ko',markersize=10)

	xmin=0.4
	xmax=0.81
	nbin=64
	gMvcut=calcMvcut(gems.z)
	x=N.compress((gems.gMv < gMvcut),gems.z)
	glir=N.compress((gems.gMv < gMvcut),gems.Lir)
	(xbin,ybin,ybinerr) = my.horizontalhist(x,xmin,xmax,nbin)

	my.drawhistpylab2(xbin,(ybin),'-','0.5')
	pylab.xlabel(r'$\rm z$',fontsize=36)
	pylab.ylabel(r'$\rm N_{gal}$',fontsize=36)
	pylab.axis([.38,.82,0.,200.])
	ax=pylab.gca()
	pylab.text(.15,.9,r'$\rm EDisCS$',fontsize=28,color='k',horizontalalignment='left',transform=ax.transAxes)
	pylab.text(.15,.8,r'$\rm GEMS$',fontsize=28,color='0.5',horizontalalignment='left',transform=ax.transAxes)
	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)


	matplotlib.rc('xtick',labelsize=24)
	matplotlib.rc('ytick',labelsize=24)


	pylab.savefig('zdist.eps')

	pylab.cla()
	pylab.clf()
	pylab.plot(x,glir,'ko')
	pylab.xlabel(r'$M_V$')
	pylab.ylabel(r'$L_{IR}$')
	ax=pylab.gca()
	ax.set_yscale('log')
	pylab.savefig('gemslir.eps')

	#calculate average fraction of LIRGS w/in r200
	mnlirg=N.sum(N.compress(z > 0.6,nlirs))
	mntot=N.sum(N.compress(z > 0.6,ntots))
	(a,b,c)=my.ratioerror(mnlirg,mntot)
	print "Fraction of LIRGS w/in R200 for z > 0.6 clusters = %5.2f + %5.2f - %5.2f (SPEC memb only)"%(a,b,c)

	mnlirg=N.sum(N.compress(z < 0.6,nlirs))
	mntot=N.sum(N.compress(z < 0.6,ntots))
	(a,b,c)=my.ratioerror(mnlirg,mntot)
	print "Fraction of LIRGS w/in R200 for z < 0.6 clusters = %5.2f + %5.2f - %5.2f  (SPEC memb only)"%(a,b,c)

	mnlirg=N.sum(N.compress(z > 0.6,nlir))
	mntot=N.sum(N.compress(z > 0.6,ntot))
	(a,b,c)=my.ratioerror(mnlirg,mntot)
	print "Fraction of LIRGS w/in R200 for z > 0.6 clusters = %5.2f + %5.2f - %5.2f (All memb)"%(a,b,c)

	mnlirg=N.sum(N.compress(z < 0.6,nlir))
	mntot=N.sum(N.compress(z < 0.6,ntot))
	(a,b,c)=my.ratioerror(mnlirg,mntot)
	print "Fraction of LIRGS w/in R200 for z < 0.6 clusters = %5.2f + %5.2f - %5.2f  (all memb)"%(a,b,c)


	mnlirg=N.sum(N.compress(z < 1.,nlir))
	mntot=N.sum(N.compress(z < 1.,ntot))
	(a,b,c)=my.ratioerror(mnlirg,mntot)
	print "Fraction of LIRGS w/in R200 for all clusters = %5.2f + %5.2f - %5.2f  (all memb)"%(a,b,c)


	#pylab.cla()
	#pylab.clf()
	#pylab.plot(z,nlir/ntot,'wo',label="Photoz-z",markersize=10)
	#pylab.plot(z,nlirs/ntots,'ko',label="Spec",markersize=10)
	#pylab.legend(loc='upper left')
	#pylab.axis([0.36,0.84,0.,.45])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	#pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=32)
	#pylab.savefig('fracIRz.eps')



	#pylab.cla()
	#pylab.clf()
	#pylab.plot(z,sumlir/ntot,'wo',label="Photo-z",markersize=10)
	#pylab.plot(z,sumlirs/ntots,'ko',label="Spec",markersize=10)
	#pylab.legend(loc='upper left')
	#pylab.axis([0.36,0.84,1.e9,2.e11])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	#ax=pylab.gca()
	#ax.set_yscale('log')
	#pylab.savefig('sumLIRz.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.plot(sigma,nlir/ntot,'wo',label='Photoz',markersize=10)
	#pylab.plot(sigma,nlirs/ntots,'ko',label='Spec',markersize=10)
	#pylab.legend(loc='upper left')
	#pylab.axis([120.,1250.,0.,.45])
	#pylab.xlabel(r'$\rm \sigma~(km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=32)
	#pylab.savefig('fracIRsigma.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.plot(sigma,sumlir/ntot,'wo',label='Photoz',markersize=10)
	#pylab.plot(sigma,sumlirs/ntots,'ko',label='Spec',markersize=10)
	#pylab.legend(loc='upper left')
	#pylab.axis([120.,1250.,1.e9,2.e11])
	#pylab.xlabel(r'$\rm \sigma~ (km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	#ax=pylab.gca()
	#ax.set_yscale('log')
	#pylab.savefig('sumLIRsigma.eps')



	fieldir=N.compress((gems.Lir > lmin) & (gems.gMv < Mvcut),gems.Lir)

	fieldirflag=N.compress((gems.gMv < Mvcut), gems.Lir)#returs array of zeros and ones, zero = Lir < lirmin; one = Lir > lirmin
	fieldirflag=(fieldirflag>lmin)#returs array of zeros and ones, zero = Lir < lirmin; one = Lir > lirmin

	fieldirflagblue=(N.compress((gems.gMv < Mvcut) & (gems.redflag < 0.1), gems.Lir) > lmin)
	fieldirflagred=(N.compress((gems.gMv < Mvcut) & (gems.redflag), gems.Lir) > lmin)

#create field sample 
	fieldirz=N.compress((gems.Lir > lmin) & (gems.gMv < Mvcut),gems.z)
	fieldirzred=N.compress((gems.Lir > lmin) & (gems.gMv < Mvcut) & (gems.redflag),gems.z)
	fieldirzblue=N.compress((gems.Lir > lmin) & (gems.gMv < Mvcut) & (gems.redflag < 0.1),gems.z)
	ngems=len(N.compress((gems.gMv < Mvcut),gems.Lir))
	ngemsz=len(N.compress((gems.gMv < Mvcut),gems.z))
	fieldz=(N.compress((gems.gMv < Mvcut),gems.z))
	fieldirall=(N.compress((gems.gMv < Mvcut),gems.Lir))

	fieldzred=(N.compress((gems.gMv < Mvcut) & (gems.redflag),gems.z))
	fieldirallred=(N.compress((gems.gMv < Mvcut) & (gems.redflag),gems.Lir))

	fieldzblue=(N.compress((gems.gMv < Mvcut) & (gems.redflag < 0.1),gems.z))
	fieldirallblue=(N.compress((gems.gMv < Mvcut) & (gems.redflag < 0.1),gems.Lir))

	nbin=5
	fieldNir=pylab.zeros(nbin,'f')
	fieldNirblue=pylab.zeros(nbin,'f')

	fieldNirred=pylab.zeros(nbin,'f')
	fieldSumir=pylab.zeros(nbin,'f')
	fieldSumirerr=pylab.zeros(nbin,'f')
	fieldNtot=pylab.zeros(nbin,'f')
	fieldNtotred=pylab.zeros(nbin,'f')
	fieldNtotblue=pylab.zeros(nbin,'f')
	fieldzbin=pylab.zeros(nbin,'f')
	fieldzbinred=pylab.zeros(nbin,'f')
	fieldzbinblue=pylab.zeros(nbin,'f')
	fieldflag=pylab.zeros(nbin,'f')

	nx=len(fieldz)
	y1=N.take(fieldirall,N.argsort(fieldz))#sort y according to x rankings
	y2=N.take(fieldirflag,N.argsort(fieldz))#sort y according to x rankings
	y3=N.take(fieldz,N.argsort(fieldz))#sort y according to x rankings
	t=len(y2)
	for i in range(nbin):
		nmin=i*int(float(nx)/float(nbin))
		nmax=(i+1)*int(float(nx)/float(nbin))
		sumir=0
		nir=0
		for j in range(nmin,nmax):
			#print "sfrmassz", j,t
			if y2[j] > 0.1:
				sumir=sumir+y1[j]
				nir=nir+1

		fieldNir[i]=nir
		fieldSumir[i]=sumir
		fieldSumirerr[i]=pylab.std(y1[nmin:nmax])
		fieldzbin[i]=pylab.average(y3[nmin:nmax])
		fieldNtot[i]=len(y3[nmin:nmax])

	nx=len(fieldzred)
	y1=N.take(fieldirallred,N.argsort(fieldzred))#sort y according to x rankings
	y2=N.take(fieldirflagred,N.argsort(fieldzred))#sort y according to x rankings
	y3=N.take(fieldzred,N.argsort(fieldzred))#sort y according to x rankings
	t=len(y2)
	for i in range(nbin):
		nmin=i*int(float(nx)/float(nbin))
		nmax=(i+1)*int(float(nx)/float(nbin))
		sumir=0
		nir=0
		for j in range(nmin,nmax):
			#print "sfrmassz", j,t
			if y2[j] > 0.1:
				sumir=sumir+y1[j]
				nir=nir+1

		fieldNirred[i]=nir
		fieldzbinred[i]=pylab.average(y3[nmin:nmax])
		fieldNtotred[i]=len(y3[nmin:nmax])

	nx=len(fieldzblue)
	y1=N.take(fieldirallblue,N.argsort(fieldzblue))#sort y according to x rankings
	y2=N.take(fieldirflagblue,N.argsort(fieldzblue))#sort y according to x rankings
	y3=N.take(fieldzblue,N.argsort(fieldzblue))#sort y according to x rankings
	t=len(y2)
	for i in range(nbin):
		nmin=i*int(float(nx)/float(nbin))
		nmax=(i+1)*int(float(nx)/float(nbin))
		sumir=0
		nir=0
		for j in range(nmin,nmax):
			#print "sfrmassz", j,t
			if y2[j] > 0.1:
				sumir=sumir+y1[j]
				nir += 1

		fieldNirblue[i]=nir
		fieldzbinblue[i]=pylab.average(y3[nmin:nmax])
		fieldNtotblue[i]=len(y3[nmin:nmax])



	#for right side plots, vs sigma
	(fieldIRfrac,errlo,errhigh)=my.ratioerror(1.*len(fieldir),(1.*ngems))
	
	fieldLirave=(N.average(fieldir))
	fieldLiraveerr=pylab.std(fieldir)/N.sqrt(1.0*len(fieldir))

	fieldLirNtot=N.sum(fieldir)/N.sqrt(1.0*ngems)
	fieldLirNtoterr=pylab.std(fieldirall)/N.sqrt(1.0*ngems)



	#BLUE and RED fractions separately
	pylab.cla()#multpanel version of 4 plots above
	pylab.clf()
	pylab.subplot(2,1,1)
	(a,b,c)=my.ratioerror(nlirblue,ntotblue)
	print 'BLUE PHOT',nlirblue,ntotblue,a
	yerrall=pylab.zeros([2,len(z)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	#pylab.errorbar(z,a,yerr=yerrall,lw=1,fmt=None,ecolor='k')
	pylab.plot(z,a,'wo',mec='k',label='_nolegend_',markersize=6)

	(xbin,y1bin,y2bin)=sffcolorsub(z,nlirblue,ntotblue)
	#print xbin,y1bin,y2bin
	(a,b,c)=my.ratioerror(y1bin,y2bin)
	yerrall=pylab.zeros([2,len(xbin)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(xbin,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(xbin,a,'wo',mec='k',label=r"$\rm All\ Blue$",markersize=18)



	(a,b,c)=my.ratioerror(nlirsblue,ntotsblue)
	print 'BLUE SPEC',nlirsblue,ntotsblue,a
	yerrall=pylab.zeros([2,len(z)],'f')
	yerrall[0]=b
	yerrall[1]=c
	#pylab.errorbar(z,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(z,a,'ko',label='_nolegend_',markersize=6)

	(xbin,y1bin,y2bin)=sffcolorsub(z,nlirsblue,ntotsblue)
	#print xbin,y1bin,y2bin
	(a,b,c)=my.ratioerror(y1bin,y2bin)
	yerrall=pylab.zeros([2,len(xbin)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(xbin,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(xbin,a,'ko',mec='k',label=r"$\rm Spec \ Blue$",markersize=18)


	(a,b,c)=my.ratioerror(fieldNirblue,fieldNtotblue)
	yerrall=pylab.zeros([2,len(fieldNirblue)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(fieldzbinblue,a,yerr=yerrall,fmt=None,ecolor='0.5')
	pylab.plot(fieldzbinblue,a,'^',color='0.5',label=r"$\rm GEMS\ Blue$",markersize=14)


	pylab.axis([0.36,0.84,-0.03,1.03])
	pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=26)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	prop=matplotlib.font_manager.FontProperties(size=14)
	pylab.legend(loc='upper left',numpoints=1,prop=prop, markerscale=0.6)

	pylab.subplot(2,1,2)
	(a,b,c)=my.ratioerror(nlirred,ntotred)
	print 'BLUE PHOT',nlirred,ntotred,a
	yerrall=pylab.zeros([2,len(z)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	#pylab.errorbar(z,a,yerr=yerrall,fmt=None,ecolor='0.5')
	pylab.plot(z,a,'wo',mec='k',label='_nolegend_',markersize=6)


	(xbin,y1bin,y2bin)=sffcolorsub(z,nlirred,ntotred)
	(a,b,c)=my.ratioerror(y1bin,y2bin)
	yerrall=pylab.zeros([2,len(xbin)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(xbin,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(xbin,a,'wo',mec='k',label=r"$\rm All\ Red$",markersize=18)

	(a,b,c)=my.ratioerror(nlirsred,ntotsred)
	print 'BLUE PHOT',nlirsred,ntotsred,a
	yerrall=pylab.zeros([2,len(z)],'f')
	yerrall[0]=b
	yerrall[1]=c
	#pylab.errorbar(z,a,yerr=yerrall,fmt=None,ecolor='0.5')
	pylab.plot(z,a,'ko',label='_nolegend_',markersize=6)

	(xbin,y1bin,y2bin)=sffcolorsub(z,nlirsred,ntotsred)
	(a,b,c)=my.ratioerror(y1bin,y2bin)
	yerrall=pylab.zeros([2,len(xbin)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(xbin,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(xbin,a,'ko',mec='k',label=r"$\rm Spec \ Red$",markersize=18)


	(a,b,c)=my.ratioerror(fieldNirred,fieldNtotred)
	yerrall=pylab.zeros([2,len(fieldNirred)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(fieldzbinred,a,yerr=yerrall,fmt=None,ecolor='0.5')
	pylab.plot(fieldzbinred,a,'^',color='0.5',label=r"$\rm GEMS \ Red$",markersize=14)

	prop=matplotlib.font_manager.FontProperties(size=14)
	pylab.legend(loc='upper left',numpoints=1,prop=prop, markerscale=0.6)
	ax=pylab.gca()
	#leg=ax.legend(fancybox=True)
	#leg.get_frame().set_alpha(0.5)

	#pylab.legend(loc='upper left')
	pylab.axis([0.36,0.84,-0.03,1.03])
	pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=26)
	pylab.xlabel(r'$\rm z$',fontsize=26)
	pylab.savefig('sffzcolor.eps')



	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)

	pylab.cla()#multpanel version of 4 plots above
	pylab.clf()
	pylab.subplot(321)
	(a,b,c)=my.ratioerror(nlir,ntot)
	yerrall=pylab.zeros([2,len(z)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(z,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(z,nlir/ntot,'wo',label=r"$\rm All$",markersize=13)
	(a,b,c)=my.ratioerror(nlirs,ntots)
	yerrall=pylab.zeros([2,len(z)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(z,a,yerr=yerrall,fmt=None,ecolor='k')
	pylab.plot(z,nlirs/ntots,'ko',label=r"$\rm Spec$",markersize=10)
	(a,b,c)=my.ratioerror(fieldNir,fieldNtot)
	yerrall=pylab.zeros([2,len(fieldNir)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(fieldzbin,a,yerr=yerrall,fmt=None,ecolor='0.5')
	pylab.plot(fieldzbin,fieldNir/fieldNtot,'^',color='0.5',label=r"$\rm GEMS$",markersize=10)
	prop=matplotlib.font_manager.FontProperties(size=14)
	pylab.legend(numpoints=1,prop=prop, markerscale=0.6)
	ax=pylab.gca()
	#leg=ax.legend(fancybox=True)
	#leg.get_frame().set_alpha(0.5)

	#pylab.legend(loc='upper left')
	pylab.axis([0.36,0.84,-0.01,.85])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=26)
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	
	(a,b,c)=my.ratioerror(N.sum(nlirs),N.sum(ntots))
	print "Average number of LIRGS w/in R200 for entire sample = %5.2f + %5.2f -%5.2f"%(a,b,c)
	#ax.set_yticklabels(([]))

	#pylab.savefig('fracIRz.eps')

	#pylab.cla()
	#pylab.clf()


	pylab.subplot(323)
	#pylab.plot(z,sumlir/ntot,'wo',label="Photo-z",markersize=10)
	#pylab.plot(z,sumlirs/ntots,'ko',label="Spec",markersize=10)
	pylab.errorbar(z,sumlir/nlir,yerr=sumlirerr/nlir,fmt=None,ecolor='k')

	pylab.plot(z,sumlir/nlir,'wo',label="All",markersize=13)

	pylab.errorbar(z,sumlirs/nlirs,yerr=sumlirerrs/nlirs,fmt=None,ecolor='k')

	pylab.plot(z,sumlirs/nlirs,'ko',label="Spec",markersize=10)
	#pylab.legend(loc='upper left')
	#pylab.axis([0.36,0.84,1.e9,2.e11])


	pylab.plot(fieldzbin,fieldSumir/fieldNir,'^',color='0.5',label="GEMS",markersize=10)
	pylab.errorbar(fieldzbin,fieldSumir/fieldNir,yerr=fieldSumirerr/fieldNir,fmt=None,ecolor='0.5')

	pylab.axis([0.36,0.84,4.e10,2.e12])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	pylab.ylabel(r'$\rm \Sigma L_{IR}/N_{IR}(L_\odot)$',fontsize=26)
	ax=pylab.gca()
	ax.set_yscale('log')
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	#pylab.savefig('sumLIRz.eps')

	pylab.subplot(322)
	(a,b,c)=my.ratioerror(nlir,ntot)
	yerrall=pylab.zeros([2,len(z)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(sigma,a,yerr=yerrall,fmt=None,ecolor='k')

	pylab.plot(sigma,nlir/ntot,'wo',label=r'$\rm \ All$',markersize=13)

	(a,b,c)=my.ratioerror(nlirs,ntots)
	yerrall=pylab.zeros([2,len(z)],'f')
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(sigma,a,yerr=yerrall,fmt=None,ecolor='k')

	pylab.plot(sigma,nlirs/ntots,'ko',label=r'$\rm \ Spec$',markersize=10)

	pylab.axhline(y=fieldIRfrac,color='0.5',ls='-',label=r'$\rm \ GEMS$')
	pylab.axhline(y=fieldIRfrac-errlo,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldIRfrac+errhigh,color='0.5',ls='--',label='_nolegend_')


	pylab.legend(loc='upper right',numpoints=2, prop=prop, markerscale=0.6)
	ax=pylab.gca()
	#leg=ax.legend(fancybox=True)
	#leg.get_frame().set_alpha(0.5)
	pylab.axis([120.,1250.,-0.01,.85])
	#pylab.xlabel(r'$\rm \sigma~(km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=32)
	#pylab.savefig('fracIRsigma.eps')
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(324)
	#pylab.plot(sigma,sumlir/ntot,'wo',label='Photoz',markersize=10)
	#pylab.plot(sigma,sumlirs/ntots,'ko',label='Spec',markersize=10)

	pylab.errorbar(sigma,sumlir/nlir,yerr=sumlirerr/nlir,fmt=None,ecolor='k')
	pylab.plot(sigma,sumlir/nlir,'wo',label='Photoz',markersize=13)

	#9/18/09 commenting line below to try to get plotting to finish
	#pylab.errorbar(sigma,sum,yerr=yerrall,fmt=None,ecolor='k')

	pylab.plot(sigma,sumlirs/nlirs,'ko',label='Spec',markersize=10)
	#pylab.legend(loc='upper left')

	pylab.axhline(y=fieldLirave,color='0.5',ls='-',label=r'$\rm GEMS$')
	pylab.axhline(y=fieldLirave-fieldLiraveerr,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldLirave+fieldLiraveerr,color='0.5',ls='--',label='_nolegend_')


	pylab.axis([120.,1250.,4.e10,2.e12])
	#pylab.xlabel(r'$\rm \sigma~ (km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=24)
	ax=pylab.gca()
	ax.set_yscale('log')
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(325)
	#pylab.plot(z,sumlir/ntot,'wo',label="Photo-z",markersize=10)
	#pylab.plot(z,sumlirs/ntots,'ko',label="Spec",markersize=10)
	pylab.errorbar(z,sumlir/ntot,yerr=sumlirerr/ntot,fmt=None,ecolor='k')
	pylab.plot(z,sumlir/ntot,'wo',label="Photo-z",markersize=13)

	pylab.errorbar(z,sumlirs/ntots,yerr=sumlirerrs/ntots,fmt=None,ecolor='k')

	pylab.plot(z,sumlirs/ntots,'ko',label="Spec",markersize=10)

	pylab.plot(fieldzbin,fieldSumir/fieldNtot,'^',color='0.5',label="GEMS",markersize=10)
	pylab.errorbar(fieldzbin,fieldSumir/fieldNtot,yerr=fieldSumirerr/fieldNtot,fmt=None,ecolor='0.5')

	#pylab.legend(loc='upper left')
	#pylab.axis([0.36,0.84,1.e9,2.e11])
	pylab.axis([0.36,0.84,1.e9,2.e11])
	pylab.xlabel(r'$\rm z$',fontsize=26)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	pylab.ylabel(r'$\rm \Sigma L_{IR}/N_{tot}(L_\odot)$',fontsize=26)
	ax=pylab.gca()
	ax.set_yscale('log')
	#pylab.savefig('sumLIRz.eps')



	pylab.subplot(326)
	#pylab.plot(sigma,sumlir/ntot,'wo',label='Photoz',markersize=10)
	#pylab.plot(sigma,sumlirs/ntots,'ko',label='Spec',markersize=10)
	pylab.errorbar(sigma,sumlir/ntot,yerr=sumlirerr/ntot,fmt=None,ecolor='k')

	pylab.plot(sigma,sumlir/ntot,'wo',label='Photoz',markersize=13)

	pylab.errorbar(sigma,sumlirs/ntots,yerr=sumlirerrs/ntots,fmt=None,ecolor='k')

	pylab.plot(sigma,sumlirs/ntots,'ko',label='Spec',markersize=10)
	#pylab.legend(loc='upper left')


	pylab.axhline(y=fieldLirNtot,color='0.5',ls='-',label=r'$\rm GEMS$')
	pylab.axhline(y=fieldLirNtot-fieldLirNtoterr,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldLirNtot+fieldLirNtoterr,color='0.5',ls='--',label='_nolegend_')


	pylab.axis([120.,1250.,1.e9,2.e11])
	pylab.xlabel(r'$\rm \sigma (km/s)$',fontsize=26)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=24)
	ax=pylab.gca()
	ax.set_yscale('log')
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	ax.set_yticklabels(([]))



	pylab.savefig('sfrmassz.eps')

	matplotlib.rc('figure',figsize='9, 5')
	pylab.cla()
	pylab.clf()
	#pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.01,hspace=0.001)
	#pylab.subplot(121)
	pylab.plot(z,downsize,'ko',markersize=10)
	a=N.compress(z< .6,z)
	b=N.compress(z< .6,downsize)

	pylab.plot(a,b,'k^',markersize=10)

	pylab.axis([0.38,0.81,-0.05,1.4])
	pylab.xlabel(r'$\rm z$',fontsize=36)
	pylab.ylabel(r'$\rm \Delta I$',fontsize=36)
		
	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)
   
	pylab.savefig('downsizez.eps')

	pylab.cla()
	pylab.clf()

	#pylab.subplot(122)
	pylab.plot(sigma,downsize,'wo',markersize=10)
	a=N.compress(z< .6,sigma)
	b=N.compress(z< .6,downsize)
	pylab.plot(a,b,'ko',markersize=10)
	pylab.xlabel(r'$\rm \sigma \ (km/s)$',fontsize=36)
	pylab.ylabel(r'$\rm \Delta I$',fontsize=36)

	pylab.axis([200.,1200.,-0.05,1.4])

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.xticks(fontsize=30)
	pylab.yticks(fontsize=30)

	pylab.savefig('downsizesigma.eps')



	#pylab.cla()
	#pylab.clf()
	#pylab.plot(z,nlirs/ntots,'ko')
	#pylab.axis([0.36,0.84,0.,.45])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	#pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=32)
	#pylab.savefig('fracIRzspec.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.plot(z,sumlirs/ntots,'ko')
	#pylab.axis([0.36,0.84,1.e9,2.e11])
	#pylab.xlabel(r'$\rm z$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	#ax=pylab.gca()
	#ax.set_yscale('log')
	#pylab.savefig('sumLIRzspec.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.plot(sigma,nlirs/ntots,'ko')
	#pylab.axis([120.,1250.,0.,.45])
	#pylab.xlabel(r'$\rm \sigma~(km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm N_{IR}/N_{tot}$',fontsize=32)
	#pylab.savefig('fracIRsigmaspec.eps')


	#pylab.cla()
	#pylab.clf()
	#pylab.plot(sigma,sumlirs/ntots,'ko')
	#pylab.axis([120.,1250.,1.e9,2.e11])
	#pylab.xlabel(r'$\rm \sigma ~(km/s)$',fontsize=32)
	#pylab.ylabel(r'$\rm \Sigma L_{IR}/(L_\odot)/N_{tot}$',fontsize=32)
	#ax=pylab.gca()
	#ax.set_yscale('log')
	#pylab.savefig('sumLIRsigmaspec.eps')


	lozdown=N.compress(z < 0.6, downsize)
	hizdown=N.compress(z > 0.6, downsize)
	a=N.average(lozdown)
	b=pylab.std(lozdown)/N.sqrt(float(len(lozdown)))
	print "Average I-band offset to brightest IR gal for z < 0.6= %5.2f +/- %5.2f"%(a,b)
	a=N.average(hizdown)
	b=pylab.std(hizdown)/N.sqrt(float(len(hizdown)))
	print "Average I-band offset to brightest IR gal for z > 0.6 = %5.2f +/- %5.2f"%(a,b)

	print "SPEC cut: downsize vs z"
	x=z
	y=downsize
	my.dospear(x,y)

	print "SPEC cut: downsize vs sigma"
	x=sigma
	y=downsize
	my.dospear(x,y)


	print "PHOTO-Z cut: frac IR vs z"
	x=z
	y=nlir/ntot
	my.dospear(x,y)

	print "SPEC-Z cut: frac IR vs z"
	x=z
	y=nlirs/ntots
	my.dospear(x,y)


	print "PHOTO-Z cut: frac IR vs sigma"
	x=sigma
	y=nlir/ntot
	my.dospear(x,y)

	print "SPEC-Z cut: frac IR vs sigma"
	x=sigma
	y=nlirs/ntots
	my.dospear(x,y)
	#now for sum LIR
	print "PHOTO-Z cut: sum IR vs z"
	x=z
	y=sumlir/ntot
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR vs z"
	x=z
	y=sumlirs/ntots
	my.dospear(x,y)

	print "PHOTO-Z cut: sum IR vs sigma"
	x=sigma
	y=sumlir/ntot
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR vs sigma"
	x=sigma
	y=sumlirs/ntots
	my.dospear(x,y)

	print "PHOTO-Z cut: sum IR/NIR vs z"
	x=z
	y=sumlir/nlir

	my.dospear(x,y)

	print "SPEC-Z cut: sum IR/NIR vs z"
	x=z
	y=sumlirs/nlirs
	#for i in range(len(x)):
	#	print x[i],y[i]
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR/NIR vs z (non-zero points only)"
	r=sumlirs/nlirs
	x=N.compress(nlirs > 0.,z)
	y=N.compress(nlirs > 0.,r)
	my.dospear(x,y)

	print "PHOTO-Z cut: sum IR/NIR vs sigma"
	x=sigma
	y=sumlir/nlir
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR/NIR vs sigma"
	x=sigma
	y=sumlirs/nlirs
	my.dospear(x,y)


	print "PHOTO-Z cut: sum IR/Ntot vs sigma"
	x=sigma
	y=sumlir/ntot
	my.dospear(x,y)


	print "SPEC-Z cut: sum IR/Ntot vs z"
	x=z
	y=sumlirs/ntots
	#for i in range(len(x)):
	#	print x[i],y[i]
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR/Ntot vs z (non-zero points only)"
	r=sumlirs/ntots
	x=N.compress(nlirs > 0.,z)
	y=N.compress(nlirs > 0.,r)
	my.dospear(x,y)

	print "PHOTO-Z cut: sum IR/Ntot vs sigma"
	x=sigma
	y=sumlir/ntot
	my.dospear(x,y)

	print "SPEC-Z cut: sum IR/Ntot vs sigma"
	x=sigma
	y=sumlirs/ntots
	my.dospear(x,y)


def plotlirMrall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	matplotlib.rc('figure',figsize=(11,11))

	xlabx=41.5
	xlaby=-.3
	ylabx=1.3
	ylaby=1.5
	ylabel='$\rm (L_{IR}/L_\odot)$'
	xlabel='$\rm M_V$'
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.14, right=.94,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.lirMrallsub()
	prop=matplotlib.font_manager.FontProperties(size=12)
	pylab.legend(loc='lower left',numpoints=1,markerscale=0.6,prop=prop)


	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	locs,labels=pylab.xticks()


	#yloc=1.2*N.ones(len(locs),'f')
	#sfrlab=[]
	#for i in range(len(labels)):
	#	s="%3.1f"%(10.**float(locs[i])*Lsol*4.5e-44)
	#	s=str(s)
	#	print locs[i],s
	#	pylab.text(float(locs[i]),1.05,s,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
		     
	pylab.subplot(4,4,2)
	cl1354.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	yaxissfr(ax)

	pylab.subplot(4,4,5)
	cl105411.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,6)
	cl1227.lirMrallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirMralla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(337)
	#cl1103.lirMrallsub()

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1353.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1037.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	yaxissfr(ax)

	#pylab.text(-0.5,xlaby,r'$\rm M_R$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(-2.3,ylaby,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	#pylab.savefig('lirMralla.eps')
	#pylab.cla()
	#pylab.clf()

	pylab.subplot(4,4,9)
	cl1232.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1411.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,11)
	cl1420.lirMrallsub()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirMrallb.eps')

	pylab.subplot(4,4,12)
	cl1301.lirMrallsub()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	yaxissfr(ax)

	pylab.subplot(4,4,13)
	cl1138.lirMrallsub()
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,14)
	cl1018.lirMrallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,15)
	cl1059.lirMrallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,16)
	cl1202.lirMrallsub()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))
	yaxissfr(ax)
	
	pylab.text(-1.,-0.4,r'$\rm M_V$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.5,2.,r'$\rm log_{10}(L_{IR}/L_\odot)$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.text(1.2,2.,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,verticalalignment='center',rotation=270,transform=ax.transAxes)

	#pylab.text(labx,xlaby,xlabel,fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,ylabel,fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('lirMralla.eps')


def sffcolorsub(x,y1,y2):
	xbin= N.zeros(4,'f')
	y1bin=N.zeros(4,'f')
	y2bin=N.zeros(4,'f')
	for i in range(4):
		j1=i*4
		j2=(i+1)*4
		xbin[i]=pylab.average(x[j1:j2])
		y1bin[i]=pylab.sum(y1[j1:j2])
		y2bin[i]=pylab.sum(y2[j1:j2])
	
	return xbin,y1bin,y2bin

def plotlirpositionsall():
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)
	matplotlib.rc('figure',figsize=(9,10))

	xlabx=41.5
	xlaby=-.3
	ylabx=1.3
	ylaby=1.5
	ylabel='$\rm log_{10}(L_{IR}/L_\odot)$'
	xlabel='$\rm M_R$'
	pylab.cla()
	pylab.clf()
	#pylab.figure(figsize=(9,10))
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(4,4,1)
	cl1216.plotpositions()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	locs,labels=pylab.xticks()


	#yloc=1.2*N.ones(len(locs),'f')
	#sfrlab=[]
	#for i in range(len(labels)):
	#	s="%3.1f"%(10.**float(locs[i])*Lsol*4.5e-44)
	#	s=str(s)
	#	print locs[i],s
	#	pylab.text(float(locs[i]),1.05,s,horizontalalignment='center',transform=ax.transAxes,fontsize=14)
		     
	pylab.subplot(4,4,2)
	cl1354.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,3)
	cl105412.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,4)
	cl1040.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,5)
	cl105411.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,6)
	cl1227.plotpositions()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirMralla.eps')

	#pylab.cla()
	#pylab.clf()
	#pylab.subplot(337)
	#cl1103.plotpositions()

	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,7)
	cl1353.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,8)
	cl1037.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.text(-0.5,xlaby,r'$\rm \Delta RA/R_{200}$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	#pylab.text(-2.3,ylaby,r'$\rm \Delta Dec/R_{200}$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	#pylab.savefig('lirpositionsa.eps')
	#pylab.cla()
	#pylab.clf()

	pylab.subplot(4,4,9)
	cl1232.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))

	pylab.subplot(4,4,10)
	cl1411.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,11)
	cl1420.plotpositions()
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	#pylab.savefig('lirMrallb.eps')

	pylab.subplot(4,4,12)
	cl1301.plotpositions()

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,13)
	cl1138.plotpositions()
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	#ax.set_yticklabels(([]))


	pylab.subplot(4,4,14)
	cl1018.plotpositions()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	pylab.subplot(4,4,15)
	cl1059.plotpositions()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))


	pylab.subplot(4,4,16)
	cl1202.plotpositions()
	#ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax=pylab.gca()
	ax.set_yticklabels(([]))

	
	pylab.text(-1.,-0.5,r'$\rm \Delta RA (min)$',fontsize=24,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-3.4,2.,r'$\rm \Delta Dec (min)$',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	#pylab.text(xlabx,xlaby,xlabel,fontsize=24,horizontalalignment='center')
	#pylab.text(ylabx,ylaby,ylabel,fontsize=24,verticalalignment='center',rotation='vertical')


	#pylab.subplot(326)
	#pylab.legend(loc='upper right')
	pylab.savefig('lirpositionsa.eps')



def plotsffdrcomb():
	pylab.cla()
	pylab.clf()
	matplotlib.rc('xtick',labelsize=18)
	matplotlib.rc('ytick',labelsize=18)

	xlabx=46.
	xlaby=-.3
	ylabx=43.2
	ylaby=2.

	xlabx=.5
	xlaby=-.2
	ylabx=-.2-1
	ylaby=1.
	xaxmax=1.9
	pylab.cla()
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.95,bottom=.15,top=0.9,wspace=0.001,hspace=0.001)
	pylab.subplot(221)


	ax=pylab.gca()
	ax.set_xticklabels(([]))
	#ax=pylab.gca()
	#ax.set_yticklabels(([]))
	lir=[]

	lirflag=[]
	dr=[]
	mr=[]
	mv=[]
	mb=[]
	photflag=[]
	photflag2=[]
	specflag=[]
	oldspecflag=[]

	mvcutall=[]
	print "total number of field galaxies w/24 micron detections is =",len(fieldlir)


	lir=lir+list(c1.Lir)+list(c2.Lir)+list(c3.Lir)+list(c4.Lir)+list(c5.Lir)+list(c6.Lir)+list(c7.Lir)+list(c8.Lir)+list(c9.Lir)+list(c10.Lir)+list(c11.Lir)+list(c12.Lir)+list(c13.Lir)+list(c14.Lir)+list(c15.Lir)+list(c16.Lir)#+list(c17.Lir)

	hizflag=N.ones(len(list(c1.Lir)+list(c2.Lir)+list(c3.Lir)+list(c4.Lir)+list(c5.Lir)+list(c6.Lir)),'i')
	lozflag=N.zeros(len(list(c7.Lir)+list(c8.Lir)+list(c9.Lir)+list(c10.Lir)+list(c11.Lir)+list(c12.Lir)+list(c13.Lir)+list(c14.Lir)+list(c15.Lir)+list(c16.Lir)),'i')#+list(c17.Lir)),'i')

	hizflag=list(hizflag)+list(lozflag)

	hizflag=N.array(hizflag,'i')


	dr=dr+list(c1.dr)+list(c2.dr)+list(c3.dr)+list(c4.dr)+list(c5.dr)+list(c6.dr)+list(c7.dr)+list(c8.dr)+list(c9.dr)+list(c10.dr)+list(c11.dr)+list(c12.dr)+list(c13.dr)+list(c14.dr)+list(c15.dr)+list(c16.dr)##+list(c17.dr)

	mr=mr+list(c1.MR)+list(c2.MR)+list(c3.MR)+list(c4.MR)+list(c5.MR)+list(c6.MR)+list(c7.MR)+list(c8.MR)+list(c9.MR)+list(c10.MR)+list(c11.MR)+list(c12.MR)+list(c13.MR)+list(c14.MR)+list(c15.MR)+list(c16.MR)#+list(c17.MR)

	mv=mv+list(c1.MV)+list(c2.MV)+list(c3.MV)+list(c4.MV)+list(c5.MV)+list(c6.MV)+list(c7.MV)+list(c8.MV)+list(c9.MV)+list(c10.MV)+list(c11.MV)+list(c12.MV)+list(c13.MV)+list(c14.MV)+list(c15.MV)+list(c16.MV)#+list(c17.MV)

	mvcutall=list(calcMvcut(c1.zcl)*N.ones(len(c1.MV)))+list(calcMvcut(c2.zcl)*N.ones(len(c2.MV)))+list(calcMvcut(c3.zcl)*N.ones(len(c3.MV)))+list(calcMvcut(c4.zcl)*N.ones(len(c4.MV)))+list(calcMvcut(c5.zcl)*N.ones(len(c5.MV)))+list(calcMvcut(c6.zcl)*N.ones(len(c6.MV)))+list(calcMvcut(c7.zcl)*N.ones(len(c7.MV)))+list(calcMvcut(c8.zcl)*N.ones(len(c8.MV)))+list(calcMvcut(c9.zcl)*N.ones(len(c9.MV)))+list(calcMvcut(c10.zcl)*N.ones(len(c10.MV)))+list(calcMvcut(c11.zcl)*N.ones(len(c11.MV)))+list(calcMvcut(c12.zcl)*N.ones(len(c12.MV)))+list(calcMvcut(c13.zcl)*N.ones(len(c13.MV)))+list(calcMvcut(c14.zcl)*N.ones(len(c14.MV)))+list(calcMvcut(c15.zcl)*N.ones(len(c15.MV)))+list(calcMvcut(c16.zcl)*N.ones(len(c16.MV)))

	mb=mb+list(c1.MB)+list(c2.MB)+list(c3.MB)+list(c4.MB)+list(c5.MB)+list(c6.MB)+list(c7.MB)+list(c8.MB)+list(c9.MB)+list(c10.MB)+list(c11.MB)+list(c12.MB)+list(c13.MB)+list(c14.MB)+list(c15.MB)+list(c16.MB)#+list(c17.MV)

	mbv=N.array(mb,'f')-N.array(mv,'f')

	lirflag=lirflag+list(c1.matchflag24)+list(c2.matchflag24)+list(c3.matchflag24)+list(c4.matchflag24)+list(c5.matchflag24)+list(c6.matchflag24)+list(c7.matchflag24)+list(c8.matchflag24)+list(c9.matchflag24)+list(c10.matchflag24)+list(c11.matchflag24)+list(c12.matchflag24)+list(c13.matchflag24)+list(c14.matchflag24)+list(c15.matchflag24)+list(c16.matchflag24)#+list(c17.matchflag24)

	photflag=photflag+list(c1.defmembflag)+list(c2.defmembflag)+list(c3.defmembflag)+list(c4.defmembflag)+list(c5.defmembflag)+list(c6.defmembflag)+list(c7.defmembflag)+list(c8.defmembflag)+list(c9.defmembflag)+list(c10.defmembflag)+list(c11.defmembflag)+list(c12.defmembflag)+list(c13.defmembflag)+list(c14.defmembflag)+list(c15.defmembflag)+list(c16.defmembflag)#+list(c17.defmembflag)
	photflag2=photflag2+list(c1.photmembflag)+list(c2.photmembflag)+list(c3.photmembflag)+list(c4.photmembflag)+list(c5.photmembflag)+list(c6.photmembflag)+list(c7.photmembflag)+list(c8.photmembflag)+list(c9.photmembflag)+list(c10.photmembflag)+list(c11.photmembflag)+list(c12.photmembflag)+list(c13.photmembflag)+list(c14.photmembflag)+list(c15.photmembflag)+list(c16.photmembflag)#+list(c17.photmembflag)

	specflag=specflag+list(c1.newspecmatchflag)+list(c2.newspecmatchflag)+list(c3.newspecmatchflag)+list(c4.newspecmatchflag)+list(c5.newspecmatchflag)+list(c6.newspecmatchflag)+list(c7.newspecmatchflag)+list(c8.newspecmatchflag)+list(c9.newspecmatchflag)+list(c10.newspecmatchflag)+list(c11.newspecmatchflag)+list(c12.newspecmatchflag)+list(c13.newspecmatchflag)+list(c14.newspecmatchflag)+list(c15.newspecmatchflag)+list(c16.newspecmatchflag)#+list(c17.newspecmatchflag)

	oldspecflag=oldspecflag+list(c1.matchflagspecediscs)+list(c2.matchflagspecediscs)+list(c3.matchflagspecediscs)+list(c4.matchflagspecediscs)+list(c5.matchflagspecediscs)+list(c6.matchflagspecediscs)+list(c7.matchflagspecediscs)+list(c8.matchflagspecediscs)+list(c9.matchflagspecediscs)+list(c10.matchflagspecediscs)+list(c11.matchflagspecediscs)+list(c12.matchflagspecediscs)+list(c13.matchflagspecediscs)+list(c14.matchflagspecediscs)+list(c15.matchflagspecediscs)+list(c16.matchflagspecediscs)#+list(c17.matchflagspecediscs)

	print "len lir,dr,lirflag = ",len(lir),len(dr),len(lirflag)
	#lir=N.compress((lirflag > -0.5),lir)
	#dr=N.compress((lirflag > -0.5),dr)
	#lirflag=N.compress((lirflag > -.5),lirflag)
	lir1=[]
	dr1=[]
	lirflag1=[]
	lir1s=[]
	dr1s=[]
	lirflag1s=[]
	photflag1=[]

	lir1hiz=[]
	dr1hiz=[]
	lirflag1hiz=[]
	lir1shiz=[]
	dr1shiz=[]
	lirflag1shiz=[]
	photflag1hiz=[]
	bv=[]
	bvs=[]
	nspecmemb=0.
	nspecmemb24=0.
	ncompl=0.#nspecmemb+photmemb=0
	ncompl24=0.#nspecmemb+photmemb24=0
	ncontam=0.#specnonmemb-photmemb
	ncontam24=0.#specnonmemb-photmemb24
	nnonmemb=0.
	n24=0.
	nphotmemb=0.
	hizflagall=[]
	hizflagalls=[]
	mv1=[]
	mv1s=[]
	for i in range(len(lir)):
		if (lirflag[i] > -0.5):#galaxy is on 24um image
                        if (mv[i] < mvcutall[i]):
       			#print i,n24,photflag[i],specflag[i],oldspecflag[i]
				if (photflag[i] > 0.5) & (photflag2[i] > 0.5):
					lir1.append(lir[i])
					dr1.append(dr[i])
					lirflag1.append(lirflag[i])
					nphotmemb += 1.
					mv1.append(mv[i])
					bv.append(mbv[i])
					if hizflag[i]:
						hizflagall.append(1)
					else:
						hizflagall.append(0)
				if specflag[i] > 0.:
					lir1s.append(lir[i])
					dr1s.append(dr[i])
					lirflag1s.append(lirflag[i])
					mv1s.append(mv[i])
					bvs.append(mbv[i])
					if hizflag[i]:
						hizflagalls.append(1)
					else:
						hizflagalls.append(0)

					if specflag[i] > 0.1:
						nspecmemb += 1.
						if lirflag[i] > 0.:
							nspecmemb24 += 1.
						if photflag[i] > 0.1:
							ncompl+=1
							if lirflag[i] > 0:
								ncompl24 +=1.
				if abs(specflag[i]) < 0.1:
					nnonmemb += 1
				#print "got a spec nonmember!",nnonmemb,photflag[i]
					if photflag[i] > 0.1:
						ncontam += 1
						if lirflag[i] > 0:
							ncontam24 +=1.

	lir1=N.array(lir1,'d')
	dr1=N.array(dr1,'d')
	lirflag1=N.array(lirflag1,'f')#photoz sample
	lir1s=N.array(lir1s,'d')#spectroscopic
	dr1s=N.array(dr1s,'d')
	lirflag1s=N.array(lirflag1s,'f')
	hizflagall=N.array(hizflagall,'i')
	hizflagalls=N.array(hizflagalls,'i')
	mv1=N.array(mv1,'d')
	mv1s=N.array(mv1s,'d')
	bv=N.array(bv,'d')
	bvs=N.array(bvs,'d')

	print "after first cut, len(lir1) = ",len(lir1),len(lir1s),len(hizflagall),N.sum(lirflag1),N.sum(lirflag1s)

	#lminloz=10.**lirminloz
	#lminhiz=10.**lirmin
	lminloz=10.**10.75#based on flux limit of CL1353, GEMS goes deeper
	lminhiz=10.**10.95#based on flux limit of GEMS galaxies (83 uJy) at z=0.8
	#lminhiz=10.**10.85#based on flux limit of GEMS galaxies (83 uJy) at z=0.8
	for i in range(len(lirflag1)):

		if (hizflagall[i] > .1):
			if (lirflag1[i] > 0.1) & (lir1[i] > lminhiz):
				lirflag1[i]=1.
			else:
				lirflag1[i]=0.
		if (hizflagall[i] < .1):
			if (lirflag1[i] > 0.1) & (lir1[i] > lminloz):
				lirflag1[i]=1
			else:
				lirflag1[i]=0.


	for i in range(len(lirflag1s)):

		if (hizflagalls[i] > .1):
			if (lirflag1s[i] > 0.1) & (lir1s[i] > lminhiz):
				lirflag1s[i]=1.
			else:
				lirflag1[i]=0.
		if (hizflagalls[i] < .1):
			if (lirflag1s[i] > 0.1) & (lir1s[i] > lminloz):
				lirflag1s[i]=1
			else:
				lirflag1s[i]=0.

	bvall=N.compress((lir1>lminhiz) & (mv1 < Mvcut),bv)
	drall=N.compress((lir1>lminhiz) & (mv1 < Mvcut),dr1)
	bvsall=N.compress((lir1s>lminhiz) & (mv1s < Mvcut),bvs)
	drsall=N.compress((lir1s>lminhiz) & (mv1s < Mvcut),dr1s)

	lir1lz=N.compress((hizflagall < 1) & (mv1 < Mvcut),lir1)
	lir1hz=N.compress((hizflagall > 0.1)  & (mv1 < Mvcut) ,lir1)

	dr1lz=N.compress((hizflagall < 1)  & (mv1 < Mvcut),dr1)
	dr1hz=N.compress((hizflagall>0.1)  & (mv1 < Mvcut),dr1)

	lirflag1lz=N.compress((hizflagall < 1)  & (mv1 < Mvcut),lirflag1)
	lirflag1hz=N.compress((hizflagall>0.1)  & (mv1 < Mvcut),lirflag1)

	lir1slz=N.compress((hizflagalls < 1)  & (mv1s < Mvcut),lir1s)
	lir1shz=N.compress((hizflagalls > 0.1)& (mv1s < Mvcut),lir1s)

	dr1slz=N.compress((hizflagalls < 1)  & (mv1s < Mvcut),dr1s)
	dr1shz=N.compress((hizflagalls > 0.1) & (mv1s < Mvcut),dr1s)

	lirflag1slz=N.compress((hizflagalls < 1)  & (mv1s < Mvcut),lirflag1s)
	lirflag1shz=N.compress((hizflagalls>0.1)  & (mv1s < Mvcut),lirflag1s)

	print "length of arrays = ",len(lir1lz),len(lir1hz),len(dr1lz),len(dr1hz),N.sum(lirflag1lz),N.sum(lirflag1hz)



	fieldirloz=N.compress((gems.lzLir > lminloz) & (gems.lzgMv < gems.lzMvcut),gems.lzLir)
	fieldirhiz=N.compress((gems.hzLir > lminhiz) & (gems.hzgMv < gems.hzMvcut),gems.hzLir)
	fieldloz=(N.average(fieldirloz))
	fieldhiz=(N.average(fieldirhiz))
	#fieldztotloz=N.compress(fieldz < 0.6,fieldztot)
	#fieldztothiz=N.compress(fieldz > 0.6,fieldztot)
	ngemslz=len(N.compress((gems.lzgMv < gems.lzMvcut),gems.lzLir))
	ngemshz=len(N.compress((gems.hzgMv < gems.hzMvcut),gems.hzLir))
	(fieldsffracloz,fielderrloz,fielderrlozup)=my.ratioerror(len(fieldirloz),ngemslz)
	print "number of loz gems ir gals and total = ",len(fieldirloz),ngemslz
	print fieldsffracloz,fielderrloz,fielderrlozup
	(fieldsffrachiz,fielderrhiz,fielderrhizup)=my.ratioerror(len(fieldirhiz),ngemshz)
	print "number of hiz gems ir gals and total = ",len(fieldirhiz),ngemshz

	nbin=5
	print "total number of 24 um sources = ",n24
	print "total number of phot members + 24 um detections = ",N.sum(lirflag1)
	print "total number of 24 um detections w/spectra = ",N.sum(lirflag1s)


	(a,b,c)=my.ratioerror(ncompl,nspecmemb)
	print "Overall completeness (nspecmemb+photmemb)/nspecmembtotal: %i / %i = %5.2f-%5.2f+%5.2f "%(ncompl,nspecmemb,a,b,c)
	(a,b,c)=my.ratioerror(ncompl24,nspecmemb24)
	print "Completeness of 24um (nspecmemb+photmemb+24)/nspecmembtotal24: %i / %i = %5.2f-%5.2f+%5.2f "%(ncompl24,nspecmemb24,a,b,c)

	(a,b,c)=my.ratioerror(ncontam,(ncompl+ncontam))
	print "Overall contamination (nspecnonmemb+photmemb)/(nspecmemb+ncontam): %i / %i = %5.2f+/-%5.2f "%(ncontam,(ncompl+ncontam),a,b)
	(a,b,c)=my.ratioerror(ncontam24,(ncompl24+ncontam24))
	print "Contamination of 24um (nspecmemb+photmemb+24)/(nspecmemb24+ncontam24): %i / %i = %5.2f+/-%5.2f+ %5.2f "%(ncontam24,(ncompl24+ncontam24),a,b,c)


	#pylab.subplot(221)
	pylab.subplot(111)
	pylab.cla()
	pylab.clf()

	(xbin,ybin,ybinerr)=my.biniterr(dr1lz,lirflag1lz,nbin)

	#pylab.plot(xbin,ybin,'k',markersize=8)
	pylab.errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='k',label='_nolegend')
	#pylab.scatter(xbin,ybin,color='w',marker='o',s=20,label='All Members')
	pylab.plot(xbin,ybin,'wo',markersize=8,label=r'$\rm \ All \ Members$')


	(xbin,ybin,ybinerr)=my.biniterr(dr1slz,lirflag1slz,nbin)
	#pylab.scatter(xbin,ybin,c='k',marker='o',s=20,label='Spec Members')
	pylab.plot(xbin,ybin,'ko',markersize=8,label=r'$\rm \ Spec \ Members$')
	pylab.errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='k',label='_nolegend_')
	pylab.axhline(y=fieldsffracloz,color='0.5',ls='-',label=r'$\rm \ GEMS$')
	pylab.axhline(y=fieldsffracloz-fielderrloz,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldsffracloz+fielderrlozup,color='0.5',ls='--',label='_nolegend_')


	pylab.legend(loc='upper right',numpoints=2)
	pylab.axis([0.,xaxmax,.0,.6])
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#pylab.xlabel('$\rm dr/R_{200}$',fontsize=24)
	#ax.set_yscale('log')
	#pylab.text(0.15,.75,r'GEMS',fontsize=18,color='0.5',horizontalalignment='left',transform=ax.transAxes)
	location=.5,.7

	#pylab.legend(loc=location)


	ax=pylab.gca()
	#pylab.text(0.5,.85,r'$0.43 < z < 0.6$, $\rm log_{10}(L_{IR}/L_\odot) > 10.75$',fontsize=16,horizontalalignment='center',transform=ax.transAxes)

	pylab.ylabel(r'$\rm SF \ Fraction$',fontsize=32)
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)
	pylab.savefig('sffdrloz.eps')
	#pylab.subplot(222)
	pylab.cla()
	pylab.clf()
	(xbin,ybin,ybinerr)=my.biniterr(dr1hz,lirflag1hz,nbin)

	#pylab.plot(xbin,ybin,'k',markersize=8)
	pylab.errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='k',label='_nolegend_')
	pylab.plot(xbin,ybin,'wo',markersize=8,label=r'$\rm \ All \ Members$')

	(xbin,ybin,ybinerr)=my.biniterr(dr1shz,lirflag1shz,nbin)
	pylab.plot(xbin,ybin,'ko',markersize=8,label=r'$\rm \ Spec \ Members$')
	pylab.errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='k',label='_nolegend_')

	pylab.axhline(y=fieldsffrachiz,color='0.5',ls='-',label=r'$\rm \ GEMS$')
	pylab.axhline(y=fieldsffrachiz-fielderrhiz,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldsffrachiz+fielderrhizup,color='0.5',ls='--',label='_nolegend_')



	pylab.legend(loc='upper left',numpoints=2)
	#pylab.xlabel('$\rm dr/R_{200}$',fontsize=24)
	#ax.set_yscale('log')
	#pylab.ylabel(r'$\rm SF \ Fraction$',fontsize=24)
	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax.set_yticklabels(([]))
	pylab.axis([0.,xaxmax,.0,.6])
	ax=pylab.gca()
	#pylab.text(0.15,.75,r'GEMS',fontsize=18,color='0.5',horizontalalignment='left',transform=ax.transAxes)
	#pylab.text(0.5,.85,r'$0.6 < z < 0.8$, $\rm log_{10}(L_{IR}/L_\odot) > 10.95$',fontsize=16,horizontalalignment='center',transform=ax.transAxes)


	pylab.ylabel(r'$\rm SF \ Fraction$',fontsize=32)
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)
	pylab.savefig('sffdrhiz.eps')


	#location=.7,.7
	#pylab.legend(loc=location)

#	#use gems galaxies instead of ediscs field
#	(fieldsffracloz,fielderrloz,fielderrlozup)=my.ratioerror(nfieldlirloz,ntotfieldloz)
#	(fieldsffrachiz,fielderrhiz,fielderrhizup)=my.ratioerror(nfieldlirhiz,ntotfieldhiz)

#	pylab.subplot(221)
#	pylab.axhline(y=fieldsffracloz,color='b',ls='-')
#	pylab.axhline(y=fieldsffracloz-fielderrloz,color='b',ls='--')
#	pylab.axhline(y=fieldsffracloz+fielderrlozup,color='b',ls='--')

#	pylab.subplot(222)
#	pylab.axhline(y=fieldsffrachiz,color='b',ls='-')
#	pylab.axhline(y=fieldsffrachiz-fielderrhiz,color='b',ls='--')
#	pylab.axhline(y=fieldsffrachiz+fielderrhizup,color='b',ls='--')




	#pylab.subplot(223)
	#pylab.cla()
	pylab.cla()
	pylab.clf()

	lir=N.compress(lirflag1lz > 0.1,lir1lz)
	dr=N.compress(lirflag1lz > 0.1,dr1lz)
	print "number of z < 0.6 galaxies above lirminloz = ",len(lir),len(dr)
	#lir=N.log10(lir)
	(xbin,ybin,ybinerr)=my.biniterr(dr,lir,nbin)

	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k',label='_nolegend_')
	pylab.plot(xbin,ybin,'wo',markersize=8,label='_nolegend_')

	lir=N.compress(lirflag1slz > 0.1,lir1slz)
	dr=N.compress(lirflag1slz > 0.1,dr1slz)
	#lir=N.log10(lir)
	(xbin,ybin,ybinerr)=my.biniterr(dr,lir,nbin)
	pylab.plot(xbin,ybin,'ko',markersize=8,label='_nolegend_')
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k',label='_nolegend_')
	#pylab.plot(dr,lir,'k.')

	#x=N.arange(0.,2.,0.1)
	#y=field*N.ones(len(x),'d')
	
	std=pylab.std(fieldirloz)/N.sqrt(1.*len(fieldlirloz))
	pylab.axhline(y=fieldloz,color='0.5',ls='-',label=r'$\rm GEMS$')
	pylab.axhline(y=fieldloz+std,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldloz-std,color='0.5',ls='--',label='_nolegend_')


	#location=.7,.7
	#pylab.legend(loc=location)
#	#add ediscs field in too
#	pylab.axhline(y=fieldlirlozave,color='b',ls='-')
#	pylab.axhline(y=fieldlirlozave+fieldlirlozstd,color='b',ls='--')
#	pylab.axhline(y=fieldlirlozave-fieldlirlozstd,color='b',ls='--')


	pylab.ylabel(r'$\rm {\bar{L}_{IR}/L_\odot}$',fontsize=32)	
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)

	ax=pylab.gca()
	ax.set_yscale('log')

	pylab.axis([0.,xaxmax,6.e10,5.e11])

	ax=pylab.gca()
	#pylab.text(0.,-.3,r'$\rm r/R_{200}$',fontsize=32,horizontalalignment='center',transform=ax.transAxes)

	pylab.text(1.2,.5,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,horizontalalignment='center',verticalalignment='center',rotation=270,transform=ax.transAxes)

	yaxissfr(ax)

	pylab.savefig('avelirdrloz.eps')


	#pylab.text(0.15,.75,r'GEMS',fontsize=18,color='0.5',horizontalalignment='left',transform=ax.transAxes)
	#pylab.subplot(224)


	pylab.cla()
	pylab.clf()

	lir=N.compress(lirflag1hz > 0.1,lir1hz)
	dr=N.compress(lirflag1hz > 0.1,dr1hz)
	#lir=N.log10(lir)
	print "number of z > 0.6 galaxies above lirminloz = ",len(lir),len(dr)
	(xbin,ybin,ybinerr)=my.biniterr(dr,lir,nbin)
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k',label='_nolegend_')
	pylab.plot(xbin,ybin,'wo',markersize=8,label='_nolegend_')
	#pylab.plot(dr,lir,'k.')

	lir=N.compress(lirflag1shz > 0.1,lir1shz)
	dr=N.compress(lirflag1shz > 0.1,dr1shz)
	#lir=N.log10(lir)
	(xbin,ybin,ybinerr)=my.biniterr(dr,lir,nbin)
	pylab.plot(xbin,ybin,'ko',markersize=8,label='_nolegend_')
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k',label='_nolegend_')


	#x=N.arange(0.,2.,0.1)
	#y=field*N.ones(len(x),'d')
	
	std=pylab.std(fieldirhiz)/N.sqrt(1.*len(fieldlirhiz))
	pylab.axhline(y=fieldhiz,color='0.5',ls='-',label=r'$\rm GEMS$')
	pylab.axhline(y=fieldhiz+std,color='0.5',ls='--',label='_nolegend_')
	pylab.axhline(y=fieldhiz-std,color='0.5',ls='--',label='_nolegend_')


	#pylab.legend(loc=location)

#	#add ediscs field in too
#	pylab.axhline(y=fieldlirhizave,color='b',ls='-')
#	pylab.axhline(y=fieldlirhizave+fieldlirhizstd,color='b',ls='--')
#	pylab.axhline(y=fieldlirhizave-fieldlirhizstd,color='b',ls='--')


	#pylab.axhline(y=gems.avelir,color='r',ls='-')
	#pylab.axhline(y=gems.avelir+gems.avelirerr,color='r',ls='--')
	#pylab.axhline(y=gems.avelir-gems.avelirerr,color='r',ls='--')
	#pylab.plot(x,y,'k-')
	#for i in range(len(x)):
	#	print x[i],y[i]



	ax=pylab.gca()
	ax.set_yscale('log')

	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	#ax.set_yticklabels(([]))
	#ax.yaxis.tick_left()

	pylab.axis([0.,xaxmax,6.e10,5.e11])
	pylab.ylabel(r'$\rm {\bar{L}_{IR}/L_\odot}$',fontsize=32)	
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)

	ax=pylab.gca()
	#pylab.text(0.,-.3,r'$\rm r/R_{200}$',fontsize=32,horizontalalignment='center',transform=ax.transAxes)

	pylab.text(1.2,.5,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,horizontalalignment='center',verticalalignment='center',rotation=270,transform=ax.transAxes)

	yaxissfr(ax)

	#pylab.text(0.15,.75,r'GEMS',fontsize=18,color='0.5',horizontalalignment='left',transform=ax.transAxes)
	#pylab.xlabel(r'$\rm dr/R_{200}$',fontsize=24)

	#pylab.subplot(224)

	#ax.set_yticklabels(([]))
	
	#pylab.text(ylabx,ylaby,r'log$_{10}$(n$\rm _{gal}$/log$_{10}$L)',fontsize=24,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('avelirdrhiz.eps')

	#pylab.savefig('sffdrcomb.eps')

def plotbvdr():
	lirmin=10.**10.95#based on flux limit of GEMS galaxies (83 uJy) at z=0.8
	nbin=5
	mv=[]
	mb=[]
	dr=[]
	mr=[]
	bv=[]
	ub=[]
	lirflag=[]
	photflag=[]
	specflag=[]
	lir=[]
	mv=mv+list(c1.MV)+list(c2.MV)+list(c3.MV)+list(c4.MV)+list(c5.MV)+list(c6.MV)+list(c7.MV)+list(c8.MV)+list(c9.MV)+list(c10.MV)+list(c11.MV)+list(c12.MV)+list(c13.MV)+list(c14.MV)+list(c15.MV)+list(c16.MV)#+list(c17.MV)
	mvcutall=list(calcMvcut(c1.zcl)*N.ones(len(c1.MV)))+list(calcMvcut(c2.zcl)*N.ones(len(c2.MV)))+list(calcMvcut(c3.zcl)*N.ones(len(c3.MV)))+list(calcMvcut(c4.zcl)*N.ones(len(c4.MV)))+list(calcMvcut(c5.zcl)*N.ones(len(c5.MV)))+list(calcMvcut(c6.zcl)*N.ones(len(c6.MV)))+list(calcMvcut(c7.zcl)*N.ones(len(c7.MV)))+list(calcMvcut(c8.zcl)*N.ones(len(c8.MV)))+list(calcMvcut(c9.zcl)*N.ones(len(c9.MV)))+list(calcMvcut(c10.zcl)*N.ones(len(c10.MV)))+list(calcMvcut(c11.zcl)*N.ones(len(c11.MV)))+list(calcMvcut(c12.zcl)*N.ones(len(c12.MV)))+list(calcMvcut(c13.zcl)*N.ones(len(c13.MV)))+list(calcMvcut(c14.zcl)*N.ones(len(c14.MV)))+list(calcMvcut(c15.zcl)*N.ones(len(c15.MV)))+list(calcMvcut(c16.zcl)*N.ones(len(c16.MV)))
	mb=mb+list(c1.MB)+list(c2.MB)+list(c3.MB)+list(c4.MB)+list(c5.MB)+list(c6.MB)+list(c7.MB)+list(c8.MB)+list(c9.MB)+list(c10.MB)+list(c11.MB)+list(c12.MB)+list(c13.MB)+list(c14.MB)+list(c15.MB)+list(c16.MB)#+list(c17.MV)

	dr=dr+list(c1.dr)+list(c2.dr)+list(c3.dr)+list(c4.dr)+list(c5.dr)+list(c6.dr)+list(c7.dr)+list(c8.dr)+list(c9.dr)+list(c10.dr)+list(c11.dr)+list(c12.dr)+list(c13.dr)+list(c14.dr)+list(c15.dr)+list(c16.dr)##+list(c17.dr)

	mr=mr+list(c1.MR)+list(c2.MR)+list(c3.MR)+list(c4.MR)+list(c5.MR)+list(c6.MR)+list(c7.MR)+list(c8.MR)+list(c9.MR)+list(c10.MR)+list(c11.MR)+list(c12.MR)+list(c13.MR)+list(c14.MR)+list(c15.MR)+list(c16.MR)#+list(c17.MR)

	#mv=mv+list(c1.MV)+list(c2.MV)+list(c3.MV)+list(c4.MV)+list(c5.MV)+list(c6.MV)+list(c7.MV)+list(c8.MV)+list(c9.MV)+list(c10.MV)+list(c11.MV)+list(c12.MV)+list(c13.MV)+list(c14.MV)+list(c15.MV)+list(c16.MV)#+list(c17.MV)
	#mb=mb+list(c1.MB)+list(c2.MB)+list(c3.MB)+list(c4.MB)+list(c5.MB)+list(c6.MB)+list(c7.MB)+list(c8.MB)+list(c9.MB)+list(c10.MB)+list(c11.MB)+list(c12.MB)+list(c13.MB)+list(c14.MB)+list(c15.MB)+list(c16.MB)#+list(c17.MV)

	#mbv=N.array(mb,'f')-N.array(mv,'f')
	bv=bv+list(c1.BVbestzclust)+list(c2.BVbestzclust)+list(c3.BVbestzclust)+list(c4.BVbestzclust)+list(c5.BVbestzclust)+list(c6.BVbestzclust)+list(c7.BVbestzclust)+list(c8.BVbestzclust)+list(c9.BVbestzclust)+list(c10.BVbestzclust)+list(c11.BVbestzclust)+list(c12.BVbestzclust)+list(c13.BVbestzclust)+list(c14.BVbestzclust)+list(c15.BVbestzclust)+list(c16.BVbestzclust)#+list(c17.MV)

	ub=ub+list(c1.UBbestzclust)+list(c2.UBbestzclust)+list(c3.UBbestzclust)+list(c4.UBbestzclust)+list(c5.UBbestzclust)+list(c6.UBbestzclust)+list(c7.UBbestzclust)+list(c8.UBbestzclust)+list(c9.UBbestzclust)+list(c10.UBbestzclust)+list(c11.UBbestzclust)+list(c12.UBbestzclust)+list(c13.UBbestzclust)+list(c14.UBbestzclust)+list(c15.UBbestzclust)+list(c16.UBbestzclust)#+list(c17.MV)

	lirflag=lirflag+list(c1.matchflag24)+list(c2.matchflag24)+list(c3.matchflag24)+list(c4.matchflag24)+list(c5.matchflag24)+list(c6.matchflag24)+list(c7.matchflag24)+list(c8.matchflag24)+list(c9.matchflag24)+list(c10.matchflag24)+list(c11.matchflag24)+list(c12.matchflag24)+list(c13.matchflag24)+list(c14.matchflag24)+list(c15.matchflag24)+list(c16.matchflag24)#+list(c17.matchflag24)

	lir=lir+list(c1.Lir)+list(c2.Lir)+list(c3.Lir)+list(c4.Lir)+list(c5.Lir)+list(c6.Lir)+list(c7.Lir)+list(c8.Lir)+list(c9.Lir)+list(c10.Lir)+list(c11.Lir)+list(c12.Lir)+list(c13.Lir)+list(c14.Lir)+list(c15.Lir)+list(c16.Lir)#+list(c17.Lir)

	photflag=photflag+list(c1.defmembflag)+list(c2.defmembflag)+list(c3.defmembflag)+list(c4.defmembflag)+list(c5.defmembflag)+list(c6.defmembflag)+list(c7.defmembflag)+list(c8.defmembflag)+list(c9.defmembflag)+list(c10.defmembflag)+list(c11.defmembflag)+list(c12.defmembflag)+list(c13.defmembflag)+list(c14.defmembflag)+list(c15.defmembflag)+list(c16.defmembflag)#+list(c17.defmembflag)

	specflag=specflag+list(c1.newspecmatchflag)+list(c2.newspecmatchflag)+list(c3.newspecmatchflag)+list(c4.newspecmatchflag)+list(c5.newspecmatchflag)+list(c6.newspecmatchflag)+list(c7.newspecmatchflag)+list(c8.newspecmatchflag)+list(c9.newspecmatchflag)+list(c10.newspecmatchflag)+list(c11.newspecmatchflag)+list(c12.newspecmatchflag)+list(c13.newspecmatchflag)+list(c14.newspecmatchflag)+list(c15.newspecmatchflag)+list(c16.newspecmatchflag)#+list(c17.newspecmatchflag)


	specflag=N.array(specflag,'f')
	photflag=N.array(photflag,'f')
	matchflag24=N.array(lirflag,'f')
	dr=N.array(dr,'f')
	lir=N.array(lir,'d')
	bv=N.array(bv,'f')
	ub=N.array(ub,'f')
	mv=N.array(mv,'f')

	drall=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(mv < mvcutall),dr)
	bvall=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(mv < mvcutall),bv)
	uball=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(mv < mvcutall),ub)



	drsall=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(specflag > 0.1)&(mv < mvcutall),dr)
	bvsall=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(specflag > 0.1)&(mv < mvcutall),bv)
	ubsall=N.compress((photflag > 0.1)&(matchflag24 > 0.1)&(lir > lirmin)&(specflag > 0.1)&(mv < mvcutall),ub)


	fieldbv=[]
	fieldub=[]
	fieldlir=[]
	fieldmatchflag24=[]
	fieldmv=[]
	fieldz=[]

	fieldbv=fieldbv+list(c1.fieldBV)+list(c2.fieldBV)+list(c3.fieldBV)+list(c4.fieldBV)+list(c5.fieldBV)+list(c6.fieldBV)+list(c7.fieldBV)+list(c8.fieldBV)+list(c9.fieldBV)+list(c10.fieldBV)+list(c11.fieldBV)+list(c12.fieldBV)+list(c13.fieldBV)+list(c14.fieldBV)+list(c15.fieldBV)+list(c16.fieldBV)#+list(c17.MV)

	fieldub=fieldub+list(c1.fieldUB)+list(c2.fieldUB)+list(c3.fieldUB)+list(c4.fieldUB)+list(c5.fieldUB)+list(c6.fieldUB)+list(c7.fieldUB)+list(c8.fieldUB)+list(c9.fieldUB)+list(c10.fieldUB)+list(c11.fieldUB)+list(c12.fieldUB)+list(c13.fieldUB)+list(c14.fieldUB)+list(c15.fieldUB)+list(c16.fieldUB)#+list(c17.MV)

	fieldmv=fieldmv+list(c1.fieldMv)+list(c2.fieldMv)+list(c3.fieldMv)+list(c4.fieldMv)+list(c5.fieldMv)+list(c6.fieldMv)+list(c7.fieldMv)+list(c8.fieldMv)+list(c9.fieldMv)+list(c10.fieldMv)+list(c11.fieldMv)+list(c12.fieldMv)+list(c13.fieldMv)+list(c14.fieldMv)+list(c15.fieldMv)+list(c16.fieldMv)#+list(c17.MV)

	fieldz=fieldz+list(c1.fieldz)+list(c2.fieldz)+list(c3.fieldz)+list(c4.fieldz)+list(c5.fieldz)+list(c6.fieldz)+list(c7.fieldz)+list(c8.fieldz)+list(c9.fieldz)+list(c10.fieldz)+list(c11.fieldz)+list(c12.fieldz)+list(c13.fieldz)+list(c14.fieldz)+list(c15.fieldz)+list(c16.fieldz)#+list(c17.MV)

	fieldlir=fieldlir+list(c1.fieldLir)+list(c2.fieldLir)+list(c3.fieldLir)+list(c4.fieldLir)+list(c5.fieldLir)+list(c6.fieldLir)+list(c7.fieldLir)+list(c8.fieldLir)+list(c9.fieldLir)+list(c10.fieldLir)+list(c11.fieldLir)+list(c12.fieldLir)+list(c13.fieldLir)+list(c14.fieldLir)+list(c15.fieldLir)+list(c16.fieldLir)#+list(c17.MV)

	fieldmatchflag24=fieldmatchflag24+list(c1.fieldmatchflag24)+list(c2.fieldmatchflag24)+list(c3.fieldmatchflag24)+list(c4.fieldmatchflag24)+list(c5.fieldmatchflag24)+list(c6.fieldmatchflag24)+list(c7.fieldmatchflag24)+list(c8.fieldmatchflag24)+list(c9.fieldmatchflag24)+list(c10.fieldmatchflag24)+list(c11.fieldmatchflag24)+list(c12.fieldmatchflag24)+list(c13.fieldmatchflag24)+list(c14.fieldmatchflag24)+list(c15.fieldmatchflag24)+list(c16.fieldmatchflag24)#+list(c17.MV)

	fieldbv=N.array(fieldbv,'f')
	fieldub=N.array(fieldub,'f')
	fieldmv=N.array(fieldmv,'f')
	fieldz=N.array(fieldz,'f')
	fieldMvcut=calcMvcut(fieldz)
	fieldlir=N.array(fieldlir,'d')
	fieldmatchflag24=N.array(fieldmatchflag24,'d')

	bv=N.compress((fieldmatchflag24 > 0.1)& (fieldlir > lirmin)&(fieldmv < fieldMvcut),fieldbv)
	print "Number of ediscs field galaxies for color comparison = ",len(bv)

	fbvave=pylab.average(bv)
	fbvstd=pylab.std(bv)/pylab.sqrt(1.*len(bv))

	bv=N.compress((fieldmatchflag24 > 0.1)& (fieldlir > lirmin)&(fieldmv < fieldMvcut),fieldub)

	fubave=pylab.average(bv)
	fubstd=pylab.std(bv)/pylab.sqrt(1.*len(bv))


	pylab.clf()
	pylab.cla()
	
	(xbin,ybin,ybinerr)=my.biniterr(drall,bvall,nbin)
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k')
	pylab.plot(xbin,ybin,'wo',markersize=8)


	(xbin,ybin,ybinerr)=my.biniterr(drsall,bvsall,nbin)
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k')
	pylab.plot(xbin,ybin,'ko',markersize=8)

	pylab.axhline(y=gems.bvave,color='0.5',ls='-')
	pylab.axhline(y=gems.bvave+gems.bvstd,color='0.5',ls='--')
	pylab.axhline(y=gems.bvave-gems.bvstd,color='0.5',ls='--')


	pylab.axhline(y=fbvave,color='b',ls='-')
	pylab.axhline(y=fbvave+fbvstd,color='b',ls='--')
	pylab.axhline(y=fbvave-fbvstd,color='b',ls='--')
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)
	pylab.ylabel(r'$\rm (B-V)_{rest}$',fontsize=32)

	pylab.savefig('bvdr.eps')


	pylab.clf()
	pylab.cla()
	
	(xbin,ybin,ybinerr)=my.biniterr(drall,uball,nbin)
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k')
	pylab.plot(xbin,ybin,'wo',markersize=8)


	(xbin,ybin,ybinerr)=my.biniterr(drsall,ubsall,nbin)
	pylab.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='k')
	pylab.plot(xbin,ybin,'ko',markersize=8)

	pylab.axhline(y=gems.ubave,color='0.5',ls='-')
	pylab.axhline(y=gems.ubave+gems.ubstd,color='0.5',ls='--')
	pylab.axhline(y=gems.ubave-gems.ubstd,color='0.5',ls='--')


	pylab.axhline(y=fubave,color='b',ls='-')
	pylab.axhline(y=fubave+fubstd,color='b',ls='--')
	pylab.axhline(y=fubave-fubstd,color='b',ls='--')
	pylab.xlabel(r'$\rm r/R_{200}$',fontsize=32)
	pylab.ylabel(r'$\rm (U-B)_{rest}$',fontsize=32)

	pylab.savefig('ubdr.eps')


def SFRHaIR():#plot SFR from Halpha vs SFR IR
	pylab.cla()
	pylab.clf()
	c=cl1216
	symbol='ko'
	symbol2='k2'
	sfrhairsub(c,symbol,symbol2,c.prefix)

	c=cl105412
	symbol='bo'
	symbol2='b2'
	sfrhairsub(c,symbol,symbol2,c.prefix)

	c=cl1040
	symbol='ro'
	symbol2='r2'
	sfrhairsub(c,symbol,symbol2,c.prefix)
	x=N.arange(0.,100.)
	y=x
	pylab.plot(x,y,'k--')

	ax=pylab.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	pylab.xlabel(r'$\rm SFR(H\alpha) \ (M_\odot/yr)$',fontsize=32)
	pylab.ylabel(r'$\rm SFR(IR) \ (M_\odot/yr)$',fontsize=32)
	pylab.savefig('SFRHaIR.eps')

def sfrhairsub(c,symbol,symbol2,label):
	sfrha=[]
	sfrhaerr=[]
	sfrir=[]
	sfrirerr=[]

	sfrha2=[]
	sfrhaerr2=[]
	sfrir2=[]
	sfrirerr2=[]
	upper=[]

	for i in range(len(c.ediscsID)):
		if c.matchflagha[i] > -.1:
			if (c.matchflag24[i]) > .1:
				upper.append(0)
				sfrha.append(c.SFR[i])
				sfrhaerr.append(c.SFRerr[i])
				sfrir.append(c.SFRir[i])
				sfrirerr.append(c.SFRirerr[i])


			else:
				upper.append(1)
				sfrha2.append(c.SFR[i])
				sfrhaerr2.append(c.SFRerr[i])
				sfrir2.append(c.SFRir[i])
				sfrirerr2.append(c.SFRirerr[i])

	pylab.errorbar(sfrha,sfrir,xerr=sfrhaerr,yerr=sfrirerr,fmt=symbol,ms=10)#,uplims=upper)
	#pylab.errorbar(sfrha2,sfrir2,xerr=sfrhaerr2,fmt=symbol,ms=8)#,uplims=True)
	pylab.errorbar(sfrha2,sfrir2,fmt=symbol,ms=8)#,uplims=True)
	#pylab.errorbar(sfrha,sfrir,xerr=sfrhaerr,fmt=symbol)#,uplims=upper)
	#pylab.errorbar(sfrha2,sfrir2,xerr=sfrhaerr2,fmt=symbol2)#,uplims=True)


def plotlirstellmass():
	lirmin=10.95
	lirminloz=10.75
	xmin=8.e8
	xmax=5.e12
	ymin=2.e10
	ymax=3.e12

	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)

	pylab.subplot(2,3,1)
	pylab.plot(gems.hzstellmass[pylab.where(gems.hzwhmips<1.1)],gems.hzLir[pylab.where(gems.hzwhmips<1.1)],'.',color='0.5')
	ghzx=gems.hzstellmass[pylab.where(gems.hzwhmips<1.1)]

	pylab.subplot(2,3,2)
	x=pylab.array([pylab.compress(gems.hzredflag,gems.hzstellmass)],'d')
	y=pylab.array([pylab.compress(gems.hzredflag,gems.hzLir)],'d')
	z=pylab.array([pylab.compress(gems.hzredflag,gems.hzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	#printing out data for Pascale
	gmipsflag=gems.hzwhmips < 1.1
	plir=gems.hzLir[gmipsflag]
	pstellmass=gems.hzstellmass[gmipsflag]
	predflag=gems.hzredflag[gmipsflag]
	outfile99=open('ForPascale.GemsHiz.dat','w')
	for i in range(len(plir)):
		s='%6.4e %6.4e %i'%(pstellmass[i],plir[i],predflag[i])
		outfile99.write(s)
	outfile99.close()

	ghzrx=x[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	ghzry=y[((z<1.1)& (y>10.**lirmin) & (x > minmass))]

	ghzrx=x[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	ghzry=y[((z<1.1)& (y>10.**lirmin) & (x > minmass))]

	pylab.subplot(2,3,3)
	x=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzstellmass)],'d')
	y=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzLir)],'d')
	z=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')


	#ghzbx=x[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	#ghzby=y[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	ghzbx=x[((z<1.1)& (y>10.**lirmin) & (x > minmasshizblue))]
	ghzby=y[((z<1.1)& (y>10.**lirmin) & (x > minmasshizblue))]


	hzxall=[]
	hzyall=[]
	hzxrall=[]
	hzyrall=[]
	hzxball=[]
	hzyball=[]
	for i in range(6):
		if i == 0:
			cl=c1
		elif i== 1:
			cl=c2
		elif i == 2:
			cl=c3
		elif i==3:
			cl=c4
		elif i == 4:
			cl=c5
		elif i == 5:
			cl=c6
		pylab.subplot(2,3,1)
		(x,y)=cl.lirstellmass('all',lirmin)
		pylab.subplot(2,3,2)
		(xr,yr)=cl.lirstellmass('red',lirmin)
		pylab.subplot(2,3,3)
		(xb,yb)=cl.lirstellmass('blue',lirmin)
		hzxrall=hzxrall+xr
		hzyrall=hzyrall+yr
		hzxball=hzxball+xb
		hzyball=hzyball+yb
		hzxall=hzxall+x
		hzyall=hzyall+y

	outfile99=open('ForPascale.Ediscs.Hiz.red.dat','w')
	for i in range(len(hzxrall)):
		s='%6.4e %6.4e \n'%(hzxrall[i],hzyrall[i])
		outfile99.write(s)
	outfile99.close()
	outfile99=open('ForPascale.Ediscs.Hiz.blue.dat','w')
	for i in range(len(hzxball)):
		s='%6.4e %6.4e \n'%(hzxball[i],hzyball[i])
		outfile99.write(s)
	outfile99.close()
	

	(a,b,c)=my.ratioerror(1.*len(hzxrall),1.*len(hzxall))
	print 'Fraction of hz lirgs that are red in B-V = %5.2f - %5.2f +%5.2f'%(a,b,c)
	pylab.subplot(2,3,4)
	pylab.plot(gems.lzstellmass[pylab.where(gems.lzwhmips<1.1)],gems.lzLir[pylab.where(gems.lzwhmips<1.1)],'.',color='0.5')
	glzx=gems.lzstellmass[pylab.where(gems.lzwhmips<1.1)]

	pylab.subplot(2,3,5)
	x=pylab.array([pylab.compress(gems.lzredflag,gems.lzstellmass)],'d')
	y=pylab.array([pylab.compress(gems.lzredflag,gems.lzLir)],'d')
	z=pylab.array([pylab.compress(gems.lzredflag,gems.lzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')


	#printing out data for Pascale
	gmipsflag=gems.lzwhmips < 1.1
	plir=gems.lzLir[gmipsflag]
	pstellmass=gems.lzstellmass[gmipsflag]
	predflag=gems.lzredflag[gmipsflag]
	outfile99=open('ForPascale.GemsLoz.dat','w')
	for i in range(len(plir)):
		s='%6.4e %6.4e %i'%(pstellmass[i],plir[i],predflag[i])
		outfile99.write(s)
	outfile99.close()


	glzrx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	glzry=y[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]


	pylab.subplot(2,3,6)
	x=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzstellmass)],'d')
	y=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzLir)],'d')
	z=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	#glzbx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	#glzby=y[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	glzbx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmasslozblue))]
	glzby=y[((z<1.1)& (y>10.**lirminloz) & (x > minmasslozblue))]


	lzxall=[]
	lzyall=[]
	lzxrall=[]
	lzyrall=[]
	lzxball=[]
	lzyball=[]
	for i in range(7,17):

		if i == 7:
			cl=c7
		if i == 8:
			cl=c8
		if i == 9:
			cl=c9
		if i == 10:
			cl=c10
		if i == 11:
			cl=c11
		if i == 12:
			cl=c12
		if i == 13:
			cl=c13
		if i == 14:
			cl=c14
		if i == 15:
			cl=c15
		if i == 16:
			cl=c16
		pylab.subplot(2,3,4)
		(x,y)=cl.lirstellmass('all',lirminloz)
		pylab.subplot(2,3,5)
		(xr,yr)=cl.lirstellmass('red',lirminloz)
		pylab.subplot(2,3,6)
		(xb,yb)=cl.lirstellmass('blue',lirminloz)

		lzxrall=lzxrall+xr
		lzyrall=lzyrall+yr
		lzxball=lzxball+xb
		lzyball=lzyball+yb

		lzxall=lzxall+x
		lzyall=lzyall+y

	outfile99=open('ForPascale.Ediscs.Liz.red.dat','w')
	for i in range(len(lzxrall)):
		s='%6.4e %6.4e \n'%(lzxrall[i],lzyrall[i])
		outfile99.write(s)
	outfile99.close()
	outfile99=open('ForPascale.Ediscs.Liz.blue.dat','w')
	for i in range(len(lzxball)):
		s='%6.4e %6.4e \n'%(lzxball[i],lzyball[i])
		outfile99.write(s)
	outfile99.close()


	(a,b,c)=my.ratioerror(len(lzxrall),len(lzxall))
	print 'Fraction of lz lirgs that are red in B-V (mass cut) = %5.2f - %5.2f +%5.2f'%(a,b,c)	
	for i in range(1,7):
		pylab.subplot(2,3,i)
		ax=pylab.gca()

		if (i < 4):
			pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=minmass,color='k',ls='--',label='_nolegend_',lw=1)
		else:
			pylab.axhline(y=10.**lirminloz,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=minmassloz,color='k',ls='--',label='_nolegend_',lw=1)

#		if (i == 3) or (i == 6):
#			pylab.axvline(x=minmassblue,color='b',ls=':',label='_nolegend_')

		ax.set_yscale('log')
		ax.set_xscale('log')
		pylab.axis([xmin,xmax,ymin,ymax])
		if (i < 4):
			ax.set_xticklabels(([]))
		if (i == 3) or (i == 6):
			pylab.text(.8,0.8,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			yaxissfr(ax)
		if (i == 2) or (i == 5):
			pylab.text(.8,0.8,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			if i < 3:
				matplotlib.rc('legend',numpoints=1,fontsize=16,markerscale=1)
				pylab.legend((r'$\rm GEMS$',r'$\rm EDisCS$'), loc='upper left',numpoints=1)
		if (i == 1) or (i == 4):
			pylab.text(.8,0.8,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			if i > 2:
				pylab.text(.2,0.8,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			else:
				pylab.text(.2,0.8,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)


	pylab.text(-.5,-0.3,r'$\rm M_* \ (M_\odot)$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm L_{IR}\ (L_\odot)$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.text(1.2,1.,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,verticalalignment='center',rotation=270,transform=ax.transAxes)

	pylab.savefig('lirstellmass.eps')


	#stats
	x1=ghzrx
	y1=pylab.log10(ghzry)
	x2=pylab.array(hzxrall,'d')
	y2=pylab.log10(pylab.array(hzyrall,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))

	print "Lir vs Stellar Mass: Red Hi z galaxies ks2d2s"
	print "GEMS: Average Stellar Mass and Lir of red galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Stellar Mass and Lir of red galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))
	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	#(a,b)=stats.ks_2samp(x1,x2)
	#print '2-sample KS test comparing hiz Red M* D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing hiz Red Lir (mass cut) D-value = %6.4f p-value=%6.4f'%(a,b)


	x1=(ghzbx)
	y1=pylab.log10(ghzby)
	x2=pylab.array(hzxball,'d')
	y2=pylab.log10(pylab.array(hzyball,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))
	print "Lir vs Stellar Mass: Blue Hi z galaxies ks2d2s"
	print "GEMS: Average Stellar Mass and Lir of blue galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Stellar Mass and Lir of blue galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))

	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob


	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing hiz Blue M* D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing hiz Blue Lir (mass cut) D-value = %6.4f p-value=%6.4f'%(a,b)


	x1=glzrx
	y1=pylab.log10(glzry)
	x2=pylab.array(lzxrall,'d')
	y2=pylab.log10(pylab.array(lzyrall,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))
	print "Lir vs Stellar Mass: Red Lo z galaxies ks2d2s"
	print "GEMS: Average Stellar Mass and Lir of red galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Stellar Mass and Lir of red galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))
	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing loz Red M* D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing loz Red Lir (mass cut) D-value = %6.4f p-value=%6.4f'%(a,b)


	x1=glzbx
	y1=pylab.log10(glzby)
	x2=pylab.array(lzxball,'d')
	y2=pylab.log10(pylab.array(lzyball,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))
	print "Lir vs Stellar Mass: Blue Lo z galaxies ks2d2s"
	print "GEMS: Average Stellar Mass and Lir of blue galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Stellar Mass and Lir of blue galaxies = M*= %3.2e (%3.2e), Lir= %3.2e(%3.2e)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))

	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing loz Blue M* D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing loz Blue Lir (mass cut) D-value = %6.4f p-value=%6.4f'%(a,b)




##	#plot cumulative distributions
##	pylab.clf()
##	pylab.cla()
##	#matplotlib.rc('figure',figsize='10, 5')
##	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.25,hspace=0.001)
##	ymin=-.02
##	ymax=1.05
##	Nbin=15
##	pylab.subplot(2,2,1)
##	pylab.hist(pylab.log10(hzxrall),Nbin,color='r',ls='solid',cumulative=True)
##	pylab.hist(pylab.log10(ghzrx),Nbin,color='r',ls='dashed',cumulative=True)
##	pylab.hist(pylab.log10(hzxball),Nbin,color='b',ls='solid',cumulative=True)
##	pylab.hist(pylab.log10(ghzbx),Nbin,color='b',ls='dashed',cumulative=True)
##	#old way to do it before coma upgrade
##	#my.plotcumulative(pylab.log10(hzxrall),Nbin,'r-',magplot=0)
##	#my.plotcumulative(pylab.log10(ghzrx),Nbin,'r:',magplot=0)
##	#my.plotcumulative(pylab.log10(hzxball),Nbin,'b-',magplot=0)
##	#my.plotcumulative(pylab.log10(ghzbx),Nbin,'b:',magplot=0)

##	print 'Number of galaxies that make M* cut'
##	print 'EDisCS hiz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(hzxall),len(hzxrall),len(hzxball))
##	print 'EDisCS loz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(lzxall),len(lzxrall),len(lzxball))
##	print 'GEMS hiz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(ghzx),len(ghzrx),len(ghzbx))
##	print 'GEMS loz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(glzx),len(glzrx),len(glzbx))
##	ax=pylab.gca()

##	#locs,labels=pylab.xticks()
##	#pylab.xticks([10.5,11.0,11.5])
##	ax.set_xticklabels(([]))
##	pylab.axis([10.5,11.5,ymin,ymax])
##	#pylab.axis([-20,-24.5,ymin,ymax])

##	pylab.subplot(2,2,3)
##	my.plotcumulative(pylab.log10(lzxrall),Nbin,'r-',magplot=0)
##	my.plotcumulative(pylab.log10(glzrx),Nbin,'r:',magplot=0)
##	my.plotcumulative(pylab.log10(lzxball),Nbin,'b-',magplot=0)
##	my.plotcumulative(pylab.log10(glzbx),Nbin,'b:',magplot=0)
##	#pylab.xticks([10.5,11.0,11.5])
##	pylab.axis([10.5,11.5,ymin,ymax])

##	#pylab.xticks(pylab.arange(min(locs),max(locs),1))
##	#pylab.axis([-20,-24.5,ymin,ymax])

##	pylab.subplot(2,2,2)
##	my.plotcumulative(pylab.log10(hzyrall),Nbin,'r-')
##	my.plotcumulative(pylab.log10(ghzry),Nbin,'r:')
##	my.plotcumulative(pylab.log10(hzyball),Nbin,'b-')
##	my.plotcumulative(pylab.log10(ghzby),Nbin,'b:')

##	ax=pylab.gca()
##	pylab.axis([10.5,12,ymin,ymax])
##	pylab.xticks(pylab.arange(10.5,12.5,.5))
##	ax.set_xticklabels(([]))

##	pylab.subplot(2,2,4)
##	my.plotcumulative(pylab.log10(lzyrall),Nbin,'r-')
##	my.plotcumulative(pylab.log10(glzry),Nbin,'r:')
##	my.plotcumulative(pylab.log10(lzyball),Nbin,'b-')
##	my.plotcumulative(pylab.log10(glzby),Nbin,'b:')
##	pylab.axis([10.5,12,ymin,ymax])
##	pylab.xticks(pylab.arange(10.5,12.5,.5))

##	matplotlib.rc('legend',numpoints=5,fontsize=12,markerscale=.5)
##	pylab.legend((r'$\rm EDisCS \ Red$',r'$\rm GEMS \ Red$',r'$\rm EDisCS \ Blue$', r'$\rm GEMS \ Blue$'), loc='lower right',numpoints=5)

##	pylab.subplot(2,2,1)
##	ax=pylab.gca()
##	pylab.text(.95,.45,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

##	pylab.subplot(2,2,2)
##	ax=pylab.gca()
##	pylab.text(.95,.45,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

##	pylab.subplot(2,2,3)
##	ax=pylab.gca()
##	pylab.text(.95,.45,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

##	pylab.subplot(2,2,4)
##	ax=pylab.gca()
##	pylab.text(.95,.45,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

##	pylab.subplot(2,2,3)
##	ax=pylab.gca()
##	pylab.text(.5,-0.25,r'$\rm log_{10}(M_* \ (M_\odot)) $',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
##	pylab.text(-.32,1.,r'$\rm Cumulative \ Distribution $',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
##	pylab.subplot(2,2,4)
##	ax=pylab.gca()
##	pylab.text(.5,-.25,r'$\rm log_{10}(L_{IR}\ (L_\odot))$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)

##	pylab.savefig('lirstellmasscum.eps')



##	#plot differential distributions
##	pylab.clf()
##	pylab.cla()
##	#matplotlib.rc('figure',figsize='10, 5')
##	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.25,hspace=0.001)
##	ymin=-.02
##	ymax=6
##	Nbin=15
##	Nbinred=10
##	pylab.subplot(2,2,1)
##	#my.plotnormhist(pylab.log10(hzxrall),Nbinred,'r-',magplot=0)
##	#my.plotnormhist(pylab.log10(ghzrx),Nbinred,'r:',magplot=0)
##	my.plotnormhist(pylab.log10(hzxball),Nbin,'b-',magplot=0)
##	my.plotnormhist(pylab.log10(ghzbx),Nbin,'b:',magplot=0)
	
##	ax=pylab.gca()

##	#locs,labels=pylab.xticks()
##	#pylab.xticks([10.5,11.0,11.5])
##	ax.set_xticklabels(([]))
##	pylab.axis([10.4,11.5,ymin,ymax])
##	#pylab.axis([-20,-24.5,ymin,ymax])

##	pylab.subplot(2,2,3)
##	#my.plotnormhist(pylab.log10(lzxrall),Nbinred,'r-',magplot=0)
##	#my.plotnormhist(pylab.log10(glzrx),Nbinred,'r:',magplot=0)
##	my.plotnormhist(pylab.log10(lzxball),Nbin,'b-',magplot=0)
##	my.plotnormhist(pylab.log10(glzbx),Nbin,'b:',magplot=0)
##	#pylab.xticks([10.5,11.0,11.5])
##	pylab.axis([10.4,11.5,ymin,ymax])

##	#pylab.xticks(pylab.arange(min(locs),max(locs),1))
##	#pylab.axis([-20,-24.5,ymin,ymax])

##	pylab.subplot(2,2,2)
##	#my.plotnormhist(pylab.log10(hzyrall),Nbinred,'r-')
##	#my.plotnormhist(pylab.log10(ghzry),Nbinred,'r:')
##	my.plotnormhist(pylab.log10(hzyball),Nbin,'b-')
##	my.plotnormhist(pylab.log10(ghzby),Nbin,'b:')

##	ax=pylab.gca()
##	pylab.axis([10.5,12,ymin,ymax])
##	pylab.xticks(pylab.arange(10.5,12.5,.5))
##	ax.set_xticklabels(([]))

##	pylab.subplot(2,2,4)
##	#my.plotnormhist(pylab.log10(lzyrall),Nbinred,'r-')
##	#my.plotnormhist(pylab.log10(glzry),Nbinred,'r:')
##	my.plotnormhist(pylab.log10(lzyball),Nbin,'b-')
##	my.plotnormhist(pylab.log10(glzby),Nbin,'b:')
##	pylab.axis([10.5,12,ymin,ymax])
##	pylab.xticks(pylab.arange(10.5,12.5,.5))

##	matplotlib.rc('legend',numpoints=5,fontsize=12,markerscale=.5)
##	pylab.legend((r'$\rm EDisCS \ Red$',r'$\rm GEMS \ Red$',r'$\rm EDisCS \ Blue$', r'$\rm GEMS \ Blue$'), loc='upper right',numpoints=5)

##	pylab.subplot(2,2,1)
##	ax=pylab.gca()
##	pylab.text(.05,.95,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

##	pylab.subplot(2,2,2)
##	ax=pylab.gca()
##	pylab.text(.05,.95,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

##	pylab.subplot(2,2,3)
##	ax=pylab.gca()
##	pylab.text(.05,.95,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

##	pylab.subplot(2,2,4)
##	ax=pylab.gca()
##	pylab.text(.05,.95,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

##	pylab.subplot(2,2,3)
##	ax=pylab.gca()
##	pylab.text(.5,-0.25,r'$\rm log_{10}(M_* \ (M_\odot)) $',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
##	pylab.text(-.32,1.,r'$\rm Probability \ Distribution $',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
##	pylab.subplot(2,2,4)
##	ax=pylab.gca()
##	pylab.text(.5,-.25,r'$\rm log_{10}(L_{IR}\ (L_\odot))$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)

##	pylab.savefig('lirstellmassdifferential.eps')

def plotssfrstellmass():
	lirmin=10.95
	lirminloz=10.75
	xmin=5.e8
	xmax=5.e12
	ymin=1.e-11
	ymax=5.e-9

	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)

	pylab.subplot(2,3,1)
	x=gems.hzstellmass[pylab.where(gems.hzwhmips<1.1)]
	y=gems.hzSFRir[pylab.where(gems.hzwhmips<1.1)]
	pylab.plot(x,y/x,'.',color='0.5')

	pylab.subplot(2,3,2)
	x=pylab.array([pylab.compress(gems.hzredflag,gems.hzstellmass)],'d')
	y=pylab.array([pylab.compress(gems.hzredflag,gems.hzSFRir)],'d')
	z=pylab.array([pylab.compress(gems.hzredflag,gems.hzwhmips)],'d')
	xp=x[pylab.where(z<1.1)]
	yp=y[pylab.where(z<1.1)]
	pylab.plot(xp,yp/xp,'.',color='0.5')

	ghzrx=x[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	ghzry=y[((z<1.1)& (y>10.**lirmin) & (x > minmass))]

	pylab.subplot(2,3,3)
	x=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzstellmass)],'d')
	y=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzSFRir)],'d')
	z=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzwhmips)],'d')
	xp=x[pylab.where(z<1.1)]
	yp=y[pylab.where(z<1.1)]
	pylab.plot(xp,yp/xp,'.',color='0.5')


	#ghzbx=x[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	#ghzby=y[((z<1.1)& (y>10.**lirmin) & (x > minmass))]
	ghzbx=x[((z<1.1)& (y>10.**lirmin) & (x > minmasshizblue))]
	ghzby=y[((z<1.1)& (y>10.**lirmin) & (x > minmasshizblue))]


	hzxall=[]
	hzyall=[]
	hzxrall=[]
	hzyrall=[]
	hzxball=[]
	hzyball=[]
	for i in range(6):
		if i == 0:
			cl=c1
		elif i== 1:
			cl=c2
		elif i == 2:
			cl=c3
		elif i==3:
			cl=c4
		elif i == 4:
			cl=c5
		elif i == 5:
			cl=c6
		pylab.subplot(2,3,1)
		(x,y)=cl.ssfrstellmass('all',lirmin)
		pylab.subplot(2,3,2)
		(xr,yr)=cl.ssfrstellmass('red',lirmin)
		pylab.subplot(2,3,3)
		(xb,yb)=cl.ssfrstellmass('blue',lirmin)
		hzxrall=hzxrall+xr
		hzyrall=hzyrall+yr
		hzxball=hzxball+xb
		hzyball=hzyball+yb
		hzxall=hzxall+x
		hzyall=hzyall+y



	(a,b,c)=my.ratioerror(1.*len(hzxrall),1.*len(hzxall))
	print 'Fraction of hz lirgs that are red in B-V = %5.2f - %5.2f +%5.2f'%(a,b,c)
	pylab.subplot(2,3,4)
	xp=gems.lzstellmass[pylab.where(gems.lzwhmips<1.1)]
	yp=gems.lzSFRir[pylab.where(gems.lzwhmips<1.1)]
	pylab.plot(xp,yp/xp,'.',color='0.5')

	pylab.subplot(2,3,5)
	x=pylab.array([pylab.compress(gems.lzredflag,gems.lzstellmass)],'d')
	y=pylab.array([pylab.compress(gems.lzredflag,gems.lzSFRir)],'d')
	z=pylab.array([pylab.compress(gems.lzredflag,gems.lzwhmips)],'d')
	xp=x[pylab.where(z<1.1)]
	yp=y[pylab.where(z<1.1)]
	pylab.plot(xp,yp/xp,'.',color='0.5')



	glzrx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	glzry=y[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]


	pylab.subplot(2,3,6)
	x=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzstellmass)],'d')
	y=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzSFRir)],'d')
	z=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzwhmips)],'d')
	xp=x[pylab.where(z<1.1)]
	yp=y[pylab.where(z<1.1)]
	pylab.plot(xp,yp/xp,'.',color='0.5')

	#glzbx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	#glzby=y[((z<1.1)& (y>10.**lirminloz) & (x > minmassloz))]
	glzbx=x[((z<1.1)& (y>10.**lirminloz) & (x > minmasslozblue))]
	glzby=y[((z<1.1)& (y>10.**lirminloz) & (x > minmasslozblue))]


	lzxall=[]
	lzyall=[]
	lzxrall=[]
	lzyrall=[]
	lzxball=[]
	lzyball=[]
	for i in range(7,17):

		if i == 7:
			cl=c7
		if i == 8:
			cl=c8
		if i == 9:
			cl=c9
		if i == 10:
			cl=c10
		if i == 11:
			cl=c11
		if i == 12:
			cl=c12
		if i == 13:
			cl=c13
		if i == 14:
			cl=c14
		if i == 15:
			cl=c15
		if i == 16:
			cl=c16
		pylab.subplot(2,3,4)
		(x,y)=cl.ssfrstellmass('all',lirminloz)
		pylab.subplot(2,3,5)
		(xr,yr)=cl.ssfrstellmass('red',lirminloz)
		pylab.subplot(2,3,6)
		(xb,yb)=cl.ssfrstellmass('blue',lirminloz)

		lzxrall=lzxrall+xr
		lzyrall=lzyrall+yr
		lzxball=lzxball+xb
		lzyball=lzyball+yb

		lzxall=lzxall+x
		lzyall=lzyall+y


	(a,b,c)=my.ratioerror(len(lzxrall),len(lzxall))
	print 'Fraction of lz lirgs that are red in B-V (mass cut) = %5.2f - %5.2f +%5.2f'%(a,b,c)	
	for i in range(1,7):
		pylab.subplot(2,3,i)
		ax=pylab.gca()
		xl=pylab.arange(8,13.,.5)
		xl=10.**xl

		if (i < 4):
			yl=(10.**lirmin)*bellconv*Lsol/xl
			pylab.plot(xl,yl,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=minmass,color='k',ls='--',label='_nolegend_',lw=1)
		else:
			yl=10.**lirminloz*bellconv*Lsol/xl
			pylab.plot(xl,yl,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=minmassloz,color='k',ls='--',label='_nolegend_',lw=1)

		if (i == 3):
			pylab.axvline(x=log10(minmasshizblue),color='b',ls=':',label='_nolegend_',lw=1)
		if (i == 6):
			pylab.axvline(x=log10(minmasslozblue),color='b',ls=':',label='_nolegend_',lw=1)

#		if (i == 3) or (i == 6):
#			pylab.axvline(x=minmassblue,color='b',ls=':',label='_nolegend_')

		ax.set_yscale('log')
		ax.set_xscale('log')
		pylab.axis([xmin,xmax,ymin,ymax])
		if (i < 4):
			ax.set_xticklabels(([]))
		if (i == 3) or (i == 6):
			pylab.text(.8,0.8,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			#yaxissfr(ax)
		if (i == 2) or (i == 5):
			pylab.text(.8,0.8,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			if i < 3:
				matplotlib.rc('legend',numpoints=1,fontsize=16,markerscale=1)
				pylab.legend((r'$\rm GEMS$',r'$\rm EDisCS$'), loc='upper left',numpoints=1)
		if (i == 1) or (i == 4):
			pylab.text(.8,0.8,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			if i > 2:
				pylab.text(.2,0.8,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			else:
				pylab.text(.2,0.8,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)


	pylab.text(-.5,-0.3,r'$\rm M_* \ (M_\odot)$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm SFR \ (M_\odot/yr)/ M_{*}\ (M_\odot)$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('ssfrstellmass.eps')



def plotfraclirstellmass():
	lirmin=10.95
	lirminloz=10.75
	xmin=8.e8
	xmax=5.e12
	ymin=-0.02
	ymax=1.05
	xmin=8.75#mass will be in log10
	xmax=12.75

	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)

	dbin=.5
	bins=pylab.arange(9.,12.75,dbin)#mass bins, in log


	pylab.subplot(2,3,1)
	#pylab.plot(gems.hzstellmass[pylab.where(gems.hzwhmips<1.1)],gems.hzLir[pylab.where(gems.hzwhmips<1.1)],'.',color='0.5')

	flag=(gems.hzwhmips < 1.1) & (gems.hzLir > 10.**lirmin)
	x1=pylab.log10(gems.hzstellmass)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,2)
	x=pylab.array(pylab.compress(gems.hzredflag,gems.hzstellmass),'d')
	y=pylab.array(pylab.compress(gems.hzredflag,gems.hzLir),'d')
	z=pylab.array(pylab.compress(gems.hzredflag,gems.hzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=pylab.log10(x)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,3)
	x=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzstellmass),'d')
	y=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzLir),'d')
	z=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=pylab.log10(x)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,4)
	flag=(gems.lzwhmips < 1.1) & (gems.lzLir > 10.**lirmin)
	x1=pylab.log10(gems.lzstellmass)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,5)
	x=pylab.array(pylab.compress(gems.lzredflag,gems.lzstellmass),'d')
	y=pylab.array(pylab.compress(gems.lzredflag,gems.lzLir),'d')
	z=pylab.array(pylab.compress(gems.lzredflag,gems.lzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=pylab.log10(x)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,6)
	x=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzstellmass),'d')
	y=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzLir),'d')
	z=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=pylab.log10(x)
	plotfracsub(bins,x1,flag,'-','0.5')

	hzxall=[]
	hzyall=[]
	hzxrall=[]
	hzyrall=[]
	hzxball=[]
	hzyball=[]
	for i in range(6):
		if i == 0:
			cl=c1
		elif i== 1:
			cl=c2
		elif i == 2:
			cl=c3
		elif i==3:
			cl=c4
		elif i == 4:
			cl=c5
		elif i == 5:
			cl=c6

		(x,y)=cl.fraclirstellmass('all',lirmin)

		(xr,yr)=cl.fraclirstellmass('red',lirmin)

		(xb,yb)=cl.fraclirstellmass('blue',lirmin)
		hzxrall=hzxrall+xr
		hzyrall=hzyrall+yr
		hzxball=hzxball+xb
		hzyball=hzyball+yb
		hzxall=hzxall+x
		hzyall=hzyall+y

	hzxrall=pylab.array(hzxrall,'d')
	hzyrall=pylab.array(hzyrall,'d')
	hzxball=pylab.array(hzxball,'d')
	hzyball=pylab.array(hzyball,'d')
	hzxall=pylab.array(hzxall,'d')
	hzyall=pylab.array(hzyall,'d')

	print 'got to 1'

	pylab.subplot(2,3,1)
	flag=(hzyall)
	x1=pylab.log10(hzxall)
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 2'

	pylab.subplot(2,3,2)
	flag=(hzyrall)
	x1=pylab.log10(hzxrall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 3'

	pylab.subplot(2,3,3)
	flag=(hzyball)
	x1=pylab.log10(hzxball)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')


	lzxall=[]
	lzyall=[]
	lzxrall=[]
	lzyrall=[]
	lzxball=[]
	lzyball=[]
	for i in range(7,17):

		if i == 7:
			cl=c7
		if i == 8:
			cl=c8
		if i == 9:
			cl=c9
		if i == 10:
			cl=c10
		if i == 11:
			cl=c11
		if i == 12:
			cl=c12
		if i == 13:
			cl=c13
		if i == 14:
			cl=c14
		if i == 15:
			cl=c15
		if i == 16:
			cl=c16
		(x,y)=cl.fraclirstellmass('all',lirminloz)
		(xr,yr)=cl.fraclirstellmass('red',lirminloz)
		(xb,yb)=cl.fraclirstellmass('blue',lirminloz)

		lzxrall=lzxrall+xr
		lzyrall=lzyrall+yr
		lzxball=lzxball+xb
		lzyball=lzyball+yb
		lzxall=lzxall+x
		lzyall=lzyall+y


	lzxrall=pylab.array(lzxrall,'d')
	lzyrall=pylab.array(lzyrall,'d')
	lzxball=pylab.array(lzxball,'d')
	lzyball=pylab.array(lzyball,'d')
	lzxall=pylab.array(lzxall,'d')
	lzyall=pylab.array(lzyall,'d')

	print 'got to 4'

	pylab.subplot(2,3,4)
	flag=(lzyall)
	x1=pylab.log10(lzxall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 5'

	pylab.subplot(2,3,5)
	flag=(lzyrall)
	x1=pylab.log10(lzxrall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 6'

	pylab.subplot(2,3,6)
	flag=(lzyball )
	x1=pylab.log10(lzxball)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')


	for i in range(1,7):
		pylab.subplot(2,3,i)
		ax=pylab.gca()

		if (i < 4):
			pylab.axvline(x=log10(minmass),color='k',ls='--',label='_nolegend_',lw=1)
		else:
			pylab.axvline(x=log10(minmassloz),color='k',ls='--',label='_nolegend_',lw=1)

		if (i == 3):
			pylab.axvline(x=log10(minmasshizblue),color='b',ls=':',label='_nolegend_',lw=1)
		if (i == 6):
			pylab.axvline(x=log10(minmasslozblue),color='b',ls=':',label='_nolegend_',lw=1)

		pylab.axis([xmin,xmax,ymin,ymax])
		pylab.xticks(pylab.arange(9.,xmax,.5),[r'$9.0$','',r'$10.0$','',r'$11.0$','',r'$12.0$'])
		if (i < 4):
			ax.set_xticklabels(([]))
		if (i == 3) or (i == 6):
			pylab.text(.2,0.75,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
		if (i == 2) or (i == 5):
			pylab.text(.2,0.75,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			if i < 3:
				ax=pylab.gca()
				pylab.text(.52,0.9,r'$\rm \bf{EDisCS}$',color='k',fontsize=12,horizontalalignment='left',transform=ax.transAxes)
				pylab.text(.52,0.85,r'$\rm \bf{GEMS}$',color='0.5',fontsize=12,horizontalalignment='left',transform=ax.transAxes)
		if (i == 1) or (i == 4):
			pylab.text(.2,0.75,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			if i > 2:
				pylab.text(.2,0.88,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			else:
				pylab.text(.2,0.88,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)


	pylab.text(-.5,-0.3,r'$\rm log_{10}(M_* \ (M_\odot))$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm f(L_{IR} > L_{min})$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('fraclirstellmass.eps')

def plotfracsub(bins,x1,flag,style,lcolor):
	all=pylab.ones(len(x1),'f')
	ybin=my.binyfixedbins(bins,x1,flag)
	ytot=my.binyfixedbins(bins,x1,all)
	(a,b,c)=my.ratioerror(ybin,ytot)
	yerrall=pylab.zeros([2,len(ybin)],'f')
	#print len(z),len(a),len(b),len(c)
	yerrall[0]=b
	yerrall[1]=c
	pylab.errorbar(bins,a,yerr=yerrall,fmt=None,ecolor=lcolor,lw=1)
	my.drawhistpylab2(bins-0.5*(bins[1]-bins[0]),a,style,lcolor) #draw histogram of binned data, x=left side of bins, y=number per bin


def plotfraclirMv():
	lirmin=10.95
	lirminloz=10.75
	xmin=8.e8
	xmax=5.e12
	ymin=-0.02
	ymax=1.05
	xmin=-19.5#mass will be in log10
	xmax=-24.5

	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)

	dbin=.5
	bins=pylab.arange(-25.,-19.,dbin)#mass bins, in log


	pylab.subplot(2,3,1)
	#pylab.plot(gems.hzstellmass[pylab.where(gems.hzwhmips<1.1)],gems.hzLir[pylab.where(gems.hzwhmips<1.1)],'.',color='0.5')

	flag=(gems.hzwhmips < 1.1) & (gems.hzLir > 10.**lirmin)
	x1=(gems.hzgMv)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,2)
	x=pylab.array(pylab.compress(gems.hzredflag,gems.hzgMv),'d')
	y=pylab.array(pylab.compress(gems.hzredflag,gems.hzLir),'d')
	z=pylab.array(pylab.compress(gems.hzredflag,gems.hzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=x
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,3)
	x=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzgMv),'d')
	y=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzLir),'d')
	z=pylab.array(pylab.compress((gems.hzredflag < 0.1),gems.hzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=x
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,4)
	flag=(gems.lzwhmips < 1.1) & (gems.lzLir > 10.**lirmin)
	x1=(gems.lzgMv)
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,5)
	x=pylab.array(pylab.compress(gems.lzredflag,gems.lzgMv),'d')
	y=pylab.array(pylab.compress(gems.lzredflag,gems.lzLir),'d')
	z=pylab.array(pylab.compress(gems.lzredflag,gems.lzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=x
	plotfracsub(bins,x1,flag,'-','0.5')

	pylab.subplot(2,3,6)
	x=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzgMv),'d')
	y=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzLir),'d')
	z=pylab.array(pylab.compress((gems.lzredflag < 0.1),gems.lzwhmips),'d')
	flag=((z<1.1)& (y>10.**lirmin))
	x1=x
	plotfracsub(bins,x1,flag,'-','0.5')

	hzxall=[]
	hzyall=[]
	hzxrall=[]
	hzyrall=[]
	hzxball=[]
	hzyball=[]
	for i in range(6):
		if i == 0:
			cl=c1
		elif i== 1:
			cl=c2
		elif i == 2:
			cl=c3
		elif i==3:
			cl=c4
		elif i == 4:
			cl=c5
		elif i == 5:
			cl=c6

		(x,y)=cl.fraclirMv('all',lirmin)

		(xr,yr)=cl.fraclirMv('red',lirmin)

		(xb,yb)=cl.fraclirMv('blue',lirmin)
		hzxrall=hzxrall+xr
		hzyrall=hzyrall+yr
		hzxball=hzxball+xb
		hzyball=hzyball+yb
		hzxall=hzxall+x
		hzyall=hzyall+y

	hzxrall=pylab.array(hzxrall,'d')
	hzyrall=pylab.array(hzyrall,'d')
	hzxball=pylab.array(hzxball,'d')
	hzyball=pylab.array(hzyball,'d')
	hzxall=pylab.array(hzxall,'d')
	hzyall=pylab.array(hzyall,'d')

	print 'got to 1'

	pylab.subplot(2,3,1)
	flag=(hzyall)
	x1=(hzxall)
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 2'

	pylab.subplot(2,3,2)
	flag=(hzyrall)
	x1=(hzxrall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 3'

	pylab.subplot(2,3,3)
	flag=(hzyball)
	x1=(hzxball)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')


	lzxall=[]
	lzyall=[]
	lzxrall=[]
	lzyrall=[]
	lzxball=[]
	lzyball=[]
	for i in range(7,17):

		if i == 7:
			cl=c7
		if i == 8:
			cl=c8
		if i == 9:
			cl=c9
		if i == 10:
			cl=c10
		if i == 11:
			cl=c11
		if i == 12:
			cl=c12
		if i == 13:
			cl=c13
		if i == 14:
			cl=c14
		if i == 15:
			cl=c15
		if i == 16:
			cl=c16
		(x,y)=cl.fraclirMv('all',lirminloz)
		(xr,yr)=cl.fraclirMv('red',lirminloz)
		(xb,yb)=cl.fraclirMv('blue',lirminloz)

		lzxrall=lzxrall+xr
		lzyrall=lzyrall+yr
		lzxball=lzxball+xb
		lzyball=lzyball+yb
		lzxall=lzxall+x
		lzyall=lzyall+y


	lzxrall=pylab.array(lzxrall,'d')
	lzyrall=pylab.array(lzyrall,'d')
	lzxball=pylab.array(lzxball,'d')
	lzyball=pylab.array(lzyball,'d')
	lzxall=pylab.array(lzxall,'d')
	lzyall=pylab.array(lzyall,'d')

	print 'got to 4'

	pylab.subplot(2,3,4)
	flag=(lzyall)
	x1=(lzxall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 5'

	pylab.subplot(2,3,5)
	flag=(lzyrall)
	x1=(lzxrall)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')

	print 'got to 6'

	pylab.subplot(2,3,6)
	flag=(lzyball )
	x1=(lzxball)
	x1=pylab.array(x1,'f')
	plotfracsub(bins,x1,flag,'-','k')


	for i in range(1,7):
		pylab.subplot(2,3,i)
		ax=pylab.gca()

		if (i < 4):
			pylab.axvline(x=Mvcuthz,color='k',ls='--',label='_nolegend_',lw=1)
		else:
			pylab.axvline(x=Mvcut,color='k',ls='--',label='_nolegend_',lw=1)

#		if (i == 3):
#			pylab.axvline(x=log10(minmasshizblue),color='b',ls=':',label='_nolegend_',lw=1)
#		if (i == 6):
#			pylab.axvline(x=log10(minmasslozblue),color='b',ls=':',label='_nolegend_',lw=1)

		pylab.axis([xmin,xmax,ymin,ymax])
		#pylab.xticks(pylab.arange(9.,xmax,.5),[r'$9.0$','',r'$10.0$','',r'$11.0$','',r'$12.0$'])
		if (i < 4):
			ax.set_xticklabels(([]))
		if (i == 3) or (i == 6):
			pylab.text(.2,0.75,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
		if (i == 2) or (i == 5):
			pylab.text(.2,0.75,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			if i < 3:
				ax=pylab.gca()
				pylab.text(.52,0.9,r'$\rm \bf{EDisCS}$',color='k',fontsize=12,horizontalalignment='left',transform=ax.transAxes)
				pylab.text(.52,0.85,r'$\rm \bf{GEMS}$',color='0.5',fontsize=12,horizontalalignment='left',transform=ax.transAxes)
		if (i == 1) or (i == 4):
			pylab.text(.2,0.75,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			if i > 2:
				pylab.text(.2,0.88,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			else:
				pylab.text(.2,0.88,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)


	pylab.text(-.5,-0.3,r'$\rm M_V$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm f(L_{IR} > L_{min})$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)

	pylab.savefig('fraclirMv.eps')




def plotlirMvredblue():
	lirmin=10.95
	lirminloz=10.75
	xmin=8.e8
	xmax=5.e12
	xmin=-17.8
	xmax=-26.2


	ymin=2.e10
	ymax=3.e12

	dm=0.4

	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)

	pylab.subplot(2,3,1)
	gems.hzgMv=pylab.array(gems.hzgMv,'f')
	pylab.plot(gems.hzgMv[pylab.where(gems.hzwhmips<1.1)],gems.hzLir[pylab.where(gems.hzwhmips<1.1)],'.',color='0.5')

	pylab.subplot(2,3,2)
	x=pylab.array([pylab.compress(gems.hzredflag,gems.hzgMv)],'d')
	y=pylab.array([pylab.compress(gems.hzredflag,gems.hzLir)],'d')
	z=pylab.array([pylab.compress(gems.hzredflag,gems.hzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	ghzrx=x[((z<1.1)& (y>10.**lirmin) & (x < Mvcuthz))]
	ghzry=y[((z<1.1)& (y>10.**lirmin) & (x < Mvcuthz))]


	pylab.subplot(2,3,3)
	x=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzgMv)],'d')
	y=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzLir)],'d')
	z=pylab.array([pylab.compress((gems.hzredflag < 0.1),gems.hzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	ghzbx=x[((z<1.1)& (y>10.**lirmin) & (x < Mvcuthz))]
	ghzby=y[((z<1.1)& (y>10.**lirmin) & (x < Mvcuthz))]



	hzxall=[]
	hzyall=[]

	hzxrall=[]
	hzyrall=[]
	hzxball=[]
	hzyball=[]


	for i in range(6):
		if i == 0:
			cl=c1
		elif i== 1:
			cl=c2
		elif i == 2:
			cl=c3
		elif i==3:
			cl=c4
		elif i == 4:
			cl=c5
		elif i == 5:
			cl=c6
		pylab.subplot(2,3,1)
		(x,y)=cl.lirstellmass2('all',lirmin)
		pylab.subplot(2,3,2)
		(xr,yr)=cl.lirstellmass2('red',lirmin)
		pylab.subplot(2,3,3)
		(xb,yb)=cl.lirstellmass2('blue',lirmin)

		hzxall=hzxall+x
		hzyall=hzyall+y


		hzxrall=hzxrall+xr
		hzyrall=hzyrall+yr
		hzxball=hzxball+xb
		hzyball=hzyball+yb



	(a,b,c)=my.ratioerror(len(hzxrall),len(hzxall))
	print 'Fraction of hz lirgs that are red in B-V (Mv cut) = %5.2f - %5.2f +%5.2f'%(a,b,c)	

	pylab.subplot(2,3,4)
	gems.lzgMv=pylab.array(gems.lzgMv,'f')
	pylab.plot(gems.lzgMv[pylab.where(gems.lzwhmips<1.1)],gems.lzLir[pylab.where(gems.lzwhmips<1.1)],'.',color='0.5')

	pylab.subplot(2,3,5)
	x=pylab.array([pylab.compress(gems.lzredflag,gems.lzgMv)],'d')
	y=pylab.array([pylab.compress(gems.lzredflag,gems.lzLir)],'d')
	z=pylab.array([pylab.compress(gems.lzredflag,gems.lzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	glzrx=x[((z<1.1)& (y>10.**lirminloz) & (x < Mvcut))]
	glzry=y[((z<1.1)& (y>10.**lirminloz) & (x < Mvcut))]



	pylab.subplot(2,3,6)
	x=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzgMv)],'d')
	y=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzLir)],'d')
	z=pylab.array([pylab.compress((gems.lzredflag < 0.1),gems.lzwhmips)],'d')
	pylab.plot(x[pylab.where(z<1.1)],y[pylab.where(z<1.1)],'.',color='0.5')

	glzbx=x[((z<1.1)& (y>10.**lirminloz) & (x < Mvcut))]
	glzby=y[((z<1.1)& (y>10.**lirminloz) & (x < Mvcut))]


	lzxall=[]
	lzyall=[]
	lzxrall=[]
	lzyrall=[]
	lzxball=[]
	lzyball=[]

	for i in range(7,17):

		if i == 7:
			cl=c7
		if i == 8:
			cl=c8
		if i == 9:
			cl=c9
		if i == 10:
			cl=c10
		if i == 11:
			cl=c11
		if i == 12:
			cl=c12
		if i == 13:
			cl=c13
		if i == 14:
			cl=c14
		if i == 15:
			cl=c15
		if i == 16:
			cl=c16
		pylab.subplot(2,3,4)
		(x,y)=cl.lirstellmass2('all',lirminloz)
		pylab.subplot(2,3,5)
		(xr,yr)=cl.lirstellmass2('red',lirminloz)
		pylab.subplot(2,3,6)
		(xb,yb)=cl.lirstellmass2('blue',lirminloz)

		lzxall=lzxall+x
		lzyall=lzyall+y

		lzxrall=lzxrall+xr
		lzyrall=lzyrall+yr
		lzxball=lzxball+xb
		lzyball=lzyball+yb

	(a,b,c)=my.ratioerror(len(lzxrall),len(lzxall))
	print 'Fraction of lz lirgs that are red in B-V (Mv cut) = %5.2f - %5.2f +%5.2f'%(a,b,c)	


	for i in range(1,7):
		pylab.subplot(2,3,i)
		ax=pylab.gca()


		if (i < 4):
			pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=Mvcut-dm,color='k',ls='--',label='_nolegend_',lw=1)
		else:
			pylab.axhline(y=10.**lirminloz,color='k',ls='--',label='_nolegend_',lw=1)
			pylab.axvline(x=Mvcut,color='k',ls='--',label='_nolegend_',lw=1)



		ax.set_yscale('log')
		#ax.set_xscale('log')
		pylab.axis([xmin,xmax,ymin,ymax])
		if (i < 4):
			ax.set_xticklabels(([]))
		else:
			ind=pylab.arange(-18,-27,1)
			#ax.set_xticklabels(('', r'$-19$', '', r'$-21$', '', r'$-23$','',r'$-25$'))
			ax.set_xticklabels(('','', r'$-25$', '', r'$-23$', '', r'$-21$','',r'$-19$'))
		if (i == 3) or (i == 6):
			pylab.text(.8,0.8,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			yaxissfr(ax)
		if (i == 2) or (i == 5):
			pylab.text(.8,0.8,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			ax.set_yticklabels(([]))
			if i < 3:
				matplotlib.rc('legend',numpoints=1,fontsize=16,markerscale=1)
				pylab.legend((r'$\rm GEMS$',r'$\rm EDisCS$'), loc='upper left',numpoints=1)
		if (i == 1) or (i == 4):
			pylab.text(.8,0.8,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			if i > 2:
				pylab.text(.2,0.8,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
			else:
				pylab.text(.2,0.8,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)


	pylab.text(-.5,-0.3,r'$\rm M_V $',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm L_{IR}\ (L_\odot)$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.text(1.2,1.,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,verticalalignment='center',rotation=270,transform=ax.transAxes)

	pylab.savefig('lirMvredblue.eps')


	#stats
	print 'Number of galaxies that make Mv cut'
	print 'EDisCS hiz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(hzxall),len(hzxrall),len(hzxball))
	print 'EDisCS loz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(lzxall),len(lzxrall),len(lzxball))
	print 'GEMS hiz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(ghzrx)+len(ghzbx),len(ghzrx),len(ghzbx))
	print 'GEMS loz: Ntot, Nred, Nblue = %4i, %4i, %4i'%(len(glzrx)+len(glzbx),len(glzrx),len(glzbx))



	pylab.subplot(2,2,1)
	x1=ghzrx
	y1=pylab.log10(ghzry)
	x2=pylab.array(hzxrall,'d')
	y2=pylab.log10(pylab.array(hzyrall,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))
	print "Lir vs Mv: Red Hi z galaxies ks2d2s"
	print "GEMS: Average Mv and Lir of red galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Mv and Lir of red galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))
	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing hiz Red Mv D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing hiz Red Lir (Mv cut) D-value = %6.4f p-value=%6.4f'%(a,b)


	x1=ghzbx
	y1=pylab.log10(ghzby)
	x2=pylab.array(hzxball,'d')
	y2=pylab.log10(pylab.array(hzyball,'d'))
	(d,prob)=nr.ks2d2s(x1,y1,len(x1),x2,y2,len(x2))
	print "Lir vs MV: Blue Hi z galaxies ks2d2s"
	print "GEMS: Average Mv and Lir of blue galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Mv and Lir of blue galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))

	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob



	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing hiz Blue Mv D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing hiz Blue Lir (Mv cut) D-value = %6.4f p-value=%6.4f'%(a,b)



	x1=glzrx
	y1=pylab.log10(glzry)
	x2=pylab.array(lzxrall,'d')
	y2=pylab.log10(pylab.array(lzyrall,'d'))
	#print x2,y2
	print "Lir vs Mv: Red Lo z galaxies ks2d2s"
	print "GEMS: Average Mv and Lir of red galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Mv and Lir of red galaxies = Mv= %3.1f (%3.1f), Lir= %3.2f(%3.2f)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))
	(d,prob)=nr.ks2d2s(-1*x1,y1,len(x1),-1*x2,y2,len(x2))
	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing loz Red Mv D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing loz Red Lir (Mv cut) D-value = %6.4f p-value=%6.4f'%(a,b)


	x1=glzbx
	y1=glzby
	x2=pylab.array(lzxball,'d')
	y2=pylab.array(lzyball,'d')
	#print x2,y2
	(d,prob)=nr.ks2d2s(-1*x1,y1,len(x1),-1*x2[y2<1.e12],y2[y2<1.e12],len(x2[y2<1.e12]))
	print "Lir vs Mv: Blue Lo z galaxies ks2d2s"
	print "GEMS: Average Mv and Lir of blue galaxies = Mv= %3.1f (%3.1f), Lir= %3.2e(%3.2e)"%(pylab.average(x1),pylab.std(x1),pylab.average(y1),pylab.std(y1))
	print "EDisCS: Average Mv and Lir of blue galaxies = MV= %3.1f (%3.1f), Lir= %3.2e(%3.2e)"%(pylab.average(x2),pylab.std(x2),pylab.average(y2),pylab.std(y2))

	print 'D = ',d
	print 'Prob that two pops are drawn from same parent = ',prob

	(a,b)=stats.ks_2samp(x1,x2)
	print '2-sample KS test comparing loz Blue Mv D-value = %6.4f p-value=%6.4f'%(a,b)

	(a,b)=stats.ks_2samp(y1,y2)
	print '2-sample KS test comparing loz Blue Lir (Mv cut) D-value = %6.4f p-value=%6.4f'%(a,b)



	#plot cumulative distributions
	pylab.clf()
	pylab.cla()
	#matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.25,hspace=0.001)
	ymin=-.02
	ymax=1.05
	Nbin=15
	pylab.subplot(2,2,1)
	my.plotcumulative(hzxrall,Nbin,'r-',magplot=1)
	my.plotcumulative(ghzrx,Nbin,'r:',magplot=1)
	my.plotcumulative(hzxball,Nbin,'b-',magplot=1)
	my.plotcumulative(ghzbx,Nbin,'b:',magplot=1)
	
	ax=pylab.gca()
	ax.set_xticklabels(([]))
	pylab.axis([-20,-24.5,ymin,ymax])
	#locs,labels=pylab.xticks()
	pylab.xticks([-20.0,-21.0,-22.0,-23.0,-24.0])
	pylab.axis([-20,-24.5,ymin,ymax])
	#pylab.axis([-20,-24.5,ymin,ymax])

	pylab.subplot(2,2,3)
	my.plotcumulative(lzxrall,Nbin,'r-',magplot=1)
	my.plotcumulative(glzrx,Nbin,'r:',magplot=1)
	my.plotcumulative(lzxball,Nbin,'b-',magplot=1)
	my.plotcumulative(glzbx,Nbin,'b:',magplot=1)
	pylab.axis([-20,-24.5,ymin,ymax])
	pylab.xticks([-20.0,-21.0,-22.0,-23.0,-24.0])
	pylab.axis([-20,-24.5,ymin,ymax])

	#pylab.xticks(pylab.arange(min(locs),max(locs),1))
	#pylab.axis([-20,-24.5,ymin,ymax])

	pylab.subplot(2,2,2)
	my.plotcumulative(pylab.log10(hzyrall),Nbin,'r-')
	my.plotcumulative(pylab.log10(ghzry),Nbin,'r:')
	my.plotcumulative(pylab.log10(hzyball),Nbin,'b-')
	my.plotcumulative(pylab.log10(ghzby),Nbin,'b:')

	ax=pylab.gca()
	pylab.axis([10.5,12,ymin,ymax])
	pylab.xticks(pylab.arange(10.5,12.5,.5))
	ax.set_xticklabels(([]))

	pylab.subplot(2,2,4)
	my.plotcumulative(pylab.log10(lzyrall),Nbin,'r-')
	my.plotcumulative(pylab.log10(glzry),Nbin,'r:')
	my.plotcumulative(pylab.log10(lzyball),Nbin,'b-')
	my.plotcumulative(pylab.log10(glzby),Nbin,'b:')
	pylab.axis([10.5,12,ymin,ymax])
	pylab.xticks(pylab.arange(10.5,12.5,.5))

	matplotlib.rc('legend',numpoints=5,fontsize=12,markerscale=.5)
	pylab.legend((r'$\rm EDisCS \ Red$',r'$\rm GEMS \ Red$',r'$\rm EDisCS \ Blue$', r'$\rm GEMS \ Blue$'), loc='lower right',numpoints=5)

	pylab.subplot(2,2,1)
	ax=pylab.gca()
	pylab.text(.95,.45,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

	pylab.subplot(2,2,2)
	ax=pylab.gca()
	pylab.text(.95,.45,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

	pylab.subplot(2,2,3)
	ax=pylab.gca()
	pylab.text(.95,.45,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

	pylab.subplot(2,2,4)
	ax=pylab.gca()
	pylab.text(.95,.45,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='right',transform=ax.transAxes)

	pylab.subplot(2,2,3)
	ax=pylab.gca()
	pylab.text(.5,-0.25,r'$\rm M_V $',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-.32,1.,r'$\rm Cumulative \ Distribution $',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.subplot(2,2,4)
	ax=pylab.gca()
	pylab.text(.5,-.25,r'$\rm log_{10}(L_{IR}\ (L_\odot))$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)

	pylab.savefig('lirMvredbluecum.eps')


	#plot differential distributions
	pylab.clf()
	pylab.cla()
	#matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.25,hspace=0.001)
	ymin=-.02
	ymax=6
	Nbin=15
	Nbinred=10
	pylab.subplot(2,2,1)
	#my.plotnormhist(pylab.log10(hzxrall),Nbinred,'r-',magplot=0)
	#my.plotnormhist(pylab.log10(ghzrx),Nbinred,'r:',magplot=0)
	my.plotnormhist(pylab.log10(hzxball),Nbin,'b-',magplot=0)
	my.plotnormhist(pylab.log10(ghzbx),Nbin,'b:',magplot=0)
	
	ax=pylab.gca()

	#locs,labels=pylab.xticks()
	#pylab.xticks([10.5,11.0,11.5])
	ax.set_xticklabels(([]))
	pylab.axis([10.4,11.5,ymin,ymax])
	#pylab.axis([-20,-24.5,ymin,ymax])

	pylab.subplot(2,2,3)
	#my.plotnormhist(pylab.log10(lzxrall),Nbinred,'r-',magplot=0)
	#my.plotnormhist(pylab.log10(glzrx),Nbinred,'r:',magplot=0)
	my.plotnormhist(pylab.log10(lzxball),Nbin,'b-',magplot=0)
	my.plotnormhist(pylab.log10(glzbx),Nbin,'b:',magplot=0)
	#pylab.xticks([10.5,11.0,11.5])
	pylab.axis([10.4,11.5,ymin,ymax])

	#pylab.xticks(pylab.arange(min(locs),max(locs),1))
	#pylab.axis([-20,-24.5,ymin,ymax])

	pylab.subplot(2,2,2)
	#my.plotnormhist(pylab.log10(hzyrall),Nbinred,'r-')
	#my.plotnormhist(pylab.log10(ghzry),Nbinred,'r:')
	my.plotnormhist(pylab.log10(hzyball),Nbin,'b-')
	my.plotnormhist(pylab.log10(ghzby),Nbin,'b:')

	ax=pylab.gca()
	pylab.axis([10.5,12,ymin,ymax])
	pylab.xticks(pylab.arange(10.5,12.5,.5))
	ax.set_xticklabels(([]))

	pylab.subplot(2,2,4)
	#my.plotnormhist(pylab.log10(lzyrall),Nbinred,'r-')
	#my.plotnormhist(pylab.log10(glzry),Nbinred,'r:')
	my.plotnormhist(pylab.log10(lzyball),Nbin,'b-')
	my.plotnormhist(pylab.log10(glzby),Nbin,'b:')
	pylab.axis([10.5,12,ymin,ymax])
	pylab.xticks(pylab.arange(10.5,12.5,.5))

	matplotlib.rc('legend',numpoints=5,fontsize=12,markerscale=.5)
	pylab.legend((r'$\rm EDisCS \ Red$',r'$\rm GEMS \ Red$',r'$\rm EDisCS \ Blue$', r'$\rm GEMS \ Blue$'), loc='upper right',numpoints=5)

	pylab.subplot(2,2,1)
	ax=pylab.gca()
	pylab.text(.05,.95,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

	pylab.subplot(2,2,2)
	ax=pylab.gca()
	pylab.text(.05,.95,r'$\rm 0.6 < z < 0.8 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

	pylab.subplot(2,2,3)
	ax=pylab.gca()
	pylab.text(.05,.95,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

	pylab.subplot(2,2,4)
	ax=pylab.gca()
	pylab.text(.05,.95,r'$\rm 0.42 < z < 0.6 $',fontsize=14,horizontalalignment='left',transform=ax.transAxes)

	pylab.subplot(2,2,3)
	ax=pylab.gca()
	pylab.text(.5,-0.25,r'$\rm log_{10}(M_* \ (M_\odot)) $',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-.32,1.,r'$\rm Probability \ Distribution $',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.subplot(2,2,4)
	ax=pylab.gca()
	pylab.text(.5,-.25,r'$\rm log_{10}(L_{IR}\ (L_\odot))$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)

	pylab.savefig('lirMvdifferential.eps')




def plotlirMvredblueold():
	lirmin=10.95
	lirminloz=10.75
	xmin=-17.8
	xmax=-26.2
	ymin=4.e10
	ymax=4.e12
	dm=0.4
	#plot for all hi z
	pylab.clf()
	pylab.cla()
	matplotlib.rc('figure',figsize='10, 5')
	pylab.subplots_adjust(left=0.2, right=.95,bottom=.2,top=0.95,wspace=0.001,hspace=0.001)
	pylab.subplot(2,3,1)
	pylab.plot(gems.hzgMv,gems.hzLir,'.',color='0.5')
	c1.lirstellmass2('all')
	c2.lirstellmass2('all')
	c3.lirstellmass2('all')
	c4.lirstellmass2('all')
	c5.lirstellmass2('all')
	c6.lirstellmass2('all')
	pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut-dm,color='k',ls='--',label='_nolegend_',lw=1)
	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')
	pylab.axis([xmin,xmax,ymin,ymax])

	#yaxissfr(ax)
	ax=pylab.gca()
	ax.set_xticklabels(([]))

	pylab.text(.8,0.8,r'$\rm All$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(.2,0.8,r'$\rm High z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)



	#plot for all hi z red
	pylab.subplot(2,3,2)
	x=pylab.compress(gems.hzredflag,gems.hzgMv)
	y=pylab.compress(gems.hzredflag,gems.hzLir)
	pylab.plot(x,y,'.',color='0.5')
	c1.lirstellmass2('red')
	c2.lirstellmass2('red')
	c3.lirstellmass2('red')
	c4.lirstellmass2('red')
	c5.lirstellmass2('red')
	c6.lirstellmass2('red')
	pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut-dm,color='k',ls='--',label='_nolegend_',lw=1)
	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')


	pylab.axis([xmin,xmax,ymin,ymax])
	pylab.text(.8,0.8,r'$\rm Red$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
	#yaxissfr(ax)

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax.set_yticklabels(([]))



	#plot for all hi z blue
	pylab.subplot(2,3,3)
	x=pylab.compress((gems.hzredflag < 0.1),gems.hzgMv)
	y=pylab.compress((gems.hzredflag < 0.1),gems.hzLir)
	pylab.plot(x,y,'.',color='0.5')

	c1.lirstellmass2('blue')
	c2.lirstellmass2('blue')
	c3.lirstellmass2('blue')
	c4.lirstellmass2('blue')
	c5.lirstellmass2('blue')
	c6.lirstellmass2('blue')
	pylab.axhline(y=10.**lirmin,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut-dm,color='k',ls='--',label='_nolegend_',lw=1)

	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')

	pylab.axis([xmin,xmax,ymin,ymax])

	ax=pylab.gca()
	ax.set_xticklabels(([]))
	ax.set_yticklabels(([]))

	pylab.text(.8,0.8,r'$\rm Blue$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
	yaxissfr(ax)



	#pylab.savefig('hzlirstellmass.eps')

	#plot for all low z
	#pylab.clf()
	#pylab.cla()
	pylab.subplot(2,3,4)


	pylab.plot(gems.lzgMv,gems.lzLir,'.',color='0.5')
	c7.lirstellmass2('all')
	c8.lirstellmass2('all')
	c9.lirstellmass2('all')
	c10.lirstellmass2('all')
	c11.lirstellmass2('all')
	c12.lirstellmass2('all')
	c13.lirstellmass2('all')
	c14.lirstellmass2('all')
	c15.lirstellmass2('all')
	c16.lirstellmass2('all')
	pylab.axhline(y=10.**lirminloz,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut,color='k',ls='--',label='_nolegend_',lw=1)
	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')
	pylab.axis([xmin,xmax,ymin,ymax])
	#yaxissfr(ax)

	pylab.text(.2,0.8,r'$\rm Low z$',fontsize=18,horizontalalignment='center',transform=ax.transAxes)
	#plot for all lo z red
	pylab.subplot(2,3,5)

	x=pylab.compress(gems.lzredflag,gems.lzgMv)
	y=pylab.compress(gems.lzredflag,gems.lzLir)
	pylab.plot(x,y,'.',color='0.5')


	c7.lirstellmass2('red')
	c8.lirstellmass2('red')
	c9.lirstellmass2('red')
	c10.lirstellmass2('red')
	c11.lirstellmass2('red')
	c12.lirstellmass2('red')
	c13.lirstellmass2('red')
	c14.lirstellmass2('red')
	c15.lirstellmass2('red')
	c16.lirstellmass2('red')
	pylab.axhline(y=10.**lirminloz,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut,color='k',ls='--',label='_nolegend_',lw=1)
	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')

	pylab.axis([xmin,xmax,ymin,ymax])
	#yaxissfr(ax)

	ax=pylab.gca()
	#ax.set_xticklabels(([]))
	ax.set_yticklabels(([]))
	matplotlib.rc('legend',numpoints=1,fontsize=18,markerscale=1)
	pylab.legend((r'$\rm GEMS$',r'$\rm EDisCS$'), loc='upper center',numpoints=1)

	matplotlib.rc('legend',numpoints=3,fontsize=20,markerscale=1)
	#plot for all lo z blue
	pylab.subplot(2,3,6)

	x=pylab.compress((gems.lzredflag < 0.1),gems.lzgMv)
	y=pylab.compress((gems.lzredflag < 0.1),gems.lzLir)
	pylab.plot(x,y,'.',color='0.5')


	c7.lirstellmass2('blue')
	c8.lirstellmass2('blue')
	c9.lirstellmass2('blue')
	c10.lirstellmass2('blue')
	c11.lirstellmass2('blue')
	c12.lirstellmass2('blue')
	c13.lirstellmass2('blue')
	c14.lirstellmass2('blue')
	c15.lirstellmass2('blue')
	c16.lirstellmass2('blue')
	pylab.axhline(y=10.**lirminloz,color='k',ls='--',label='_nolegend_',lw=1)
	pylab.axvline(x=Mvcut,color='k',ls='--',label='_nolegend_',lw=1)
	ax=pylab.gca()
	ax.set_yscale('log')
	#ax.set_xscale('log')


	pylab.axis([xmin,xmax,ymin,ymax])
	ax.set_yticklabels(([]))
	yaxissfr(ax)
	

	#ax.set_xticklabels(([]))


	pylab.text(-.5,-0.3,r'$\rm M_V$',fontsize=30,horizontalalignment='center',transform=ax.transAxes)
	pylab.text(-2.5,1.,r'$\rm L_{IR} \ (L_\odot)$',fontsize=30,verticalalignment='center',rotation='vertical',transform=ax.transAxes)
	pylab.text(1.2,1.,r'$\rm SFR \ (M_\odot/yr)$',fontsize=30,verticalalignment='center',rotation=270,transform=ax.transAxes)

	pylab.savefig('lirMvredblue.eps')



def writeapj():
	apjtab=open('/home/rfinn/research/clusters/spitzer/papers/paper2/tab1.tex','w')
	#apjtab.write("\\begin{deluxetable}{lcccccrrrrrrr} \n")
	apjtab.write("\\begin{deluxetable}{lccccccc} \n")
	apjtab.write("\\tablecaption{Summary of 24\micron \ Detections \label{detect}}\n")
	#apjtab.write("\\tablehead{  & & & & & & \multicolumn{3}{c}{All} & &\multicolumn{3}{c}{dr $< R_{200}$} \\\\ \colhead{Cluster} &\colhead{z} &\colhead{$\sigma$} &\colhead{f80} & \colhead{$\\rm L_{IR}(80)$} & \colhead{SFR$_{80}$}&\colhead{N$_{spec}$} & \colhead{N$_{phot}$} & \colhead{N$_{total}$} &$R_{200}$ & \colhead{N$_{spec}$} & \colhead{N$_{phot}$} & \colhead{N$_{total}$} \\\\ &  & (km/s) & ($\mu$Jy) & log$\\rm _{10}(L_{IR}/L_\odot)$& (\\smy) &   & & & (\\arcmin)} \n")
	apjtab.write("\\tablehead{  \colhead{Cluster} &\colhead{z} &\colhead{$\sigma$} & \colhead{$\\rm f_{min}^a$} & \colhead{$\\rm f_{max}^b$}&\colhead{$\\rm f_{80}^c$} & \\colhead{$\\rm L_{IR}(80)^d$} & \colhead{SFR$_{80}^e$}\\\\ &  & (km/s) & ($\mu$Jy) & ($\mu$Jy) & ($\mu$Jy) & log$\\rm _{10}(L_{IR}/L_\odot)$& (\\smy)} \n")
	apjtab.write('\\startdata \n')
	nspec=0
	nphot=0
	nspec2=0
	nphot2=0
	(n1,n2,n3,n4)=c1.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c2.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c3.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c4.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c5.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c6.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c7.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c8.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c9.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c10.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c11.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c12.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c13.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c14.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c15.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	(n1,n2,n3,n4)=c16.writeapjtbl(apjtab)
	nspec += n1
	nphot += n2
	nspec2 += n3
	nphot2 += n4
	#(n1,n2,n3,n4)=c17.writeapjtbl(apjtab)
	#nspec += n1
	#nphot += n2
	#nspec2 += n3
	#nphot2 += n4
	#apjtab.write('\\\\ \n')
	s='{\\bf Total} & &  & & & & %i & %i &%i & &%i & %i & %i\\\\ \n'%(nspec,nphot,(nspec+nphot),nspec2,nphot2,(nspec2+nphot2))
	#apjtab.write(s)
	apjtab.write('\enddata \n')
	apjtab.write('\\tablenotetext{a}{Minimum 24\micron \ flux detected from sources with SNR$>$2.5} \n')
	apjtab.write('\\tablenotetext{b}{Maximum 24\micron \ flux detected from sources with SNR$>$2.5.}\n')
	apjtab.write('\\tablenotetext{c}{24\micron \ flux corresponding to 80\% completeness limit.}\n')
	apjtab.write('\\tablenotetext{d}{L$\\rm_{IR}$ corresponding to 80\% completeness limit.  Relative error is the same as for f$_{80}$.}\n')
	apjtab.write('\\tablenotetext{e}{SFR corresponding to 80\% completeness limit.  Relative error is the same as for f$_{80}$.} \n')
	s='\end{deluxetable}\n'
	apjtab.write(s)

def compareBV():

	mb=mb+list(c1.MB)+list(c2.MB)+list(c3.MB)+list(c4.MB)+list(c5.MB)+list(c6.MB)+list(c7.MB)+list(c8.MB)+list(c9.MB)+list(c10.MB)+list(c11.MB)+list(c12.MB)+list(c13.MB)+list(c14.MB)+list(c15.MB)+list(c16.MB)#+list(c17.MB)

	mv=mv+list(c1.MV)+list(c2.MV)+list(c3.MV)+list(c4.MV)+list(c5.MV)+list(c6.MV)+list(c7.MV)+list(c8.MV)+list(c9.MV)+list(c10.MV)+list(c11.MV)+list(c12.MV)+list(c13.MV)+list(c14.MV)+list(c15.MV)+list(c16.MV)#+list(c17.MV)



def o2stats():
	nOIIno24all=N.array([c1.nOIIno24,c2.nOIIno24,c3.nOIIno24,c4.nOIIno24,c5.nOIIno24,c6.nOIIno24,c7.nOIIno24,c8.nOIIno24,c9.nOIIno24,c10.nOIIno24,c11.nOIIno24,c12.nOIIno24,c13.nOIIno24,c14.nOIIno24,c15.nOIIno24,c16.nOIIno24],'f')
	
	nOIIall=N.array([c1.nOII,c2.nOII,c3.nOII,c4.nOII,c5.nOII,c6.nOII,c7.nOII,c8.nOII,c9.nOII,c10.nOII,c11.nOII,c12.nOII,c13.nOII,c14.nOII,c15.nOII,c16.nOII],'f')
	
	n24specnoOIIall=N.array([c1.n24specnoOII,c2.n24specnoOII,c3.n24specnoOII,c4.n24specnoOII,c5.n24specnoOII,c6.n24specnoOII,c7.n24specnoOII,c8.n24specnoOII,c9.n24specnoOII,c10.n24specnoOII,c11.n24specnoOII,c12.n24specnoOII,c13.n24specnoOII,c14.n24specnoOII,c15.n24specnoOII,c16.n24specnoOII],'f')
	
	n24specall=N.array([c1.n24spec,c2.n24spec,c3.n24spec,c4.n24spec,c5.n24spec,c6.n24spec,c7.n24spec,c8.n24spec,c9.n24spec,c10.n24spec,c11.n24spec,c12.n24spec,c13.n24spec,c14.n24spec,c15.n24spec,c16.n24spec],'f')
	
	print "OII/24 statistics for entire sample"
	
	n=N.sum(nOIIno24all)
	nt=N.sum(nOIIall)
	
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb OII w/no 24 = %i, Numb OII = %i, frac = %6.3f - %6.3f +%6.3f"%(n,nt,a,b,c)
	

	n=N.sum(n24specnoOIIall)
	nt=N.sum(n24specall)
	
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb 24-detected spec memb w/ OII = %i, Numb OII = %i, frac = %6.3f - %6.3f +%6.3f"%(n,nt,a,b,c)

	nOIIno24=nOIIno24all
	nOII=nOIIall
	n24noOII=n24specnoOIIall
	n24=n24specall


	nspecall=N.array([c1.nspec,c2.nspec,c3.nspec,c4.nspec,c5.nspec,c6.nspec,c7.nspec,c8.nspec,c9.nspec,c10.nspec,c11.nspec,c12.nspec,c13.nspec,c14.nspec,c15.nspec,c16.nspec],'f')


	(fraco2,errfraco2,errfraco2hi)=my.ratioerror(nOII,nspecall)
	(frac24,errfrac24,errfrac24hi)=my.ratioerror(n24,nspecall)
	print "[OII] frac vs IR frac"
	x=frac24
	y=fraco2
	my.dospear(x,y)

	pylab.cla()
	pylab.clf()
	errfraco2=zip(errfraco2,errfraco2hi)
	errfrac24=zip(errfrac24,errfrac24hi)
	pylab.errorbar(frac24,fraco2,yerr=errfraco2,xerr=errfrac24,fmt='ko')
	pylab.xlabel(r'IR Fraction',fontsize=32)
	pylab.ylabel(r'[OII] Fraction',fontsize=32)
	pylab.savefig('o2fracvsIRfrac.eps')






##################  Main   #######################

(Lid,Lmorph,Lvistype)=readLaiMorph()

gems=GemsGalaxies()
vandanaout=open('IDs24noOII.dat','w')
gocl105412=0.
gocl105411=0.
gocl1037=0.
gocl1040=0.
gocl1216=0.    
gocl1227=1.
gocl1232=1.
gocl1354=0.
gocl1018=0.
gocl1059=0.
gocl1103=0.
gocl1138=0.
gocl1202=0.
gocl1301=0.
gocl1353=0.
gocl1411=0.
gocl1420=0.

#use a dictionary to store and access average ra and dec offsets to apply to 24um coordinates
mydra24={'cl1216':.22,'cl1354': 0.22,'cl105412': 0.22,'cl1040':0.37,'cl105411': 0.26,'cl1227': 0.44,'cl1353':-0.20,'cl1037': 0.45,'cl1232': 0.42,'cl1411': 0.14,'cl1420': 0.12,'cl1301': 0.51,'cl1138': 0.15,'cl1018': 0.30,'cl1059': 0.54,'cl1202': 0.11,'cl1103':0.0}
myddec24={'cl1216':.19,'cl1354': 0.04,'cl105412': -0.57,'cl1040':0.39,'cl105411': 0.73,'cl1227': 0.99,'cl1353':0.62,'cl1037': 0.07,'cl1232': 0.69,'cl1411': 0.58,'cl1420': 0.17,'cl1301': -0.27,'cl1138': 0.03,'cl1018': -0.43,'cl1059': -0.27,'cl1202': 0.06,'cl1103':0.0}
EdiscsName={'cl1216':'CL1216.8$-$1201','cl1354': 'CL1354.2$-$1230','cl105412': 'CL1054.7$-$1245','cl1040':'CL1040.7$-$1155','cl105411': 'CL1054.4$-$1146','cl1227': 'CL1227.9$-$1138','cl1353':'CL1353.0$-$1137','cl1037':'CL1037.9$-$1243','cl1232': 'CL1232.5$-$1250','cl1411': 'CL1411.1$-$1148','cl1420': 'CL1420.3$-$1236','cl1301': 'CL1301.7$-$1139','cl1138': 'CL1138.2$-$1133','cl1018': 'CL1018.8$-$1211','cl1059': 'CL1059.2$-$1253','cl1202': 'CL1202.7$-$1224','cl1103':'CL1103.7$-$1245'}

ErrorF80={'cl1216':'5','cl1354': '3','cl105412': '3','cl1040':'3','cl105411': '3','cl1227': '3','cl1353':'3','cl1037': '3','cl1232': '3','cl1411': '3','cl1420': '5','cl1301': '3','cl1138': '3','cl1018': '3','cl1059': '3','cl1202': '3','cl1103':'3'}

allclusters=1.
if allclusters > .1:
	gocl105412=1.
	gocl105411=1.
	
	
	gocl1037=1.
	gocl1040=1.
	gocl1216=1.    
	gocl1227=1.
	gocl1232=1.
	gocl1354=1.
	
	gocl1018=1.
	gocl1059=1.
	gocl1103=1.
	gocl1138=1.
	gocl1202=1.
	gocl1301=1.
	gocl1353=1.
	gocl1411=1.
	gocl1420=1.

if (gocl1216 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.8#apcor for
	apcor3=2.0
	apcor2=2.56
	f80=75.#90.#flux in uJy at 80% completeness
	z=0.7943
	sigma=1018.
	errsigmap=73.
	errsigmam=77.
	prefix='cl1216'
	fullprefix='cl1216-1201'
	tablename='CL1216.'
	idra='1216'
	tabra='1216.'
	iddec='1201'
	rvalue='r13810688'
	print "CL1216"
	mipsimage='/home/rfinn/research/clusters/spitzer/cl1216/mips/24/r13810688/ch1/bcd/pbcd/Combine/mosaic.fits'
	#mipsimage='/home/rfinn/research/clusters/spitzer/cl1216/hinz/mosaic.fits'#use Joannah's image
	rac='184.1879'
	decc='-12.0214'
	cl1216 = ediscs24()
	cl1216.gotoitHa(catalogpath,'cl1216','1216','1201','cl1216.v1.spe','cl1216.id',mipsimage)
	#cl1216.readoptIRphot(catalogpath+'photometry/cl1216-1201vrijk_optIRphot.v23.cat')
	##cl1216.read24('/home/rfinn/research/clusters/spitzer/MasterTables/cl1216mosaic_extract.tbl',mipsimage)
	#cl1216.read24('/home/rfinn/research/clusters/spitzer/cl1216/hinz/output_apex_step2/mosaic_extract.tbl',mipsimage)#use Joannah's image
	#print "reading irac"
	#cl1216.readirac('cl1216')
	##cl1216.read24('/home/rfinn/research/clusters/spitzer/cl1216/mips/24/r13810688/ch1/pbcd/output/mosaic_extract.tbl')#SSC pbcd version
	#print "reading photoz catalog"
	#cl1216.readphotoz(catalogpath+'photoz/zcat.final.cl1216-1201vrijk.V2.4.1.dat')
	#print "reading H-alpha file"
	#cl1216.readhafile(catalogpath+'hafiles/cl1216sfrtableapj.dat')
        ## ##cl1216.x2=cl1216.x2+40.
        ## ##cl1216.y2=cl1216.y2+38.
        ## ##cl1216.idra=cl1216.ra+0.8/3600.*15.
        ## ##cl1216.raha=cl1216.raha+0.8/3600.*15.
	#print "reading gim2d morphology"
	#cl1216.readgimmorphology(catalogpath+'morphology/HST/GIM2D/cl1216.id',catalogpath+'morphology/HST/GIM2D/cl1216-1201_drz_sci.struct2.cat')
	#print "reading visible morphology"
	#cl1216.readvismorphology(catalogpath+'morphology/HST/Visual/cl1216_1201.f.cat')
	#print "reading spectroscopy"
	#cl1216.readspectroscopy(catalogpath+'spectroscopy/cl1216.v1.spe')
	#print "running memb()"
	#cl1216.memb()
	#print "magnitude to flux conversion"
	#cl1216.magnitudetofluxconversion('cl1216',1.)
	#print "writing master table"
	#cl1216.mastertable('cl1216')
	##cl1216.readwide24('/home/rfinn/research/clusters/spitzer/desai-mosaic_extract.tbl','cl1216wide-SSC.eps')
	##cl1216.readwide24('/home/rfinn/research/clusters/spitzer/desai-mosaic_extract.tbl','cl1216wide.eps')
if (gocl1040 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.87#apcor for
	apcor3=2.01
	apcor2=2.57
	f80=80.#85.#flux in uJy at 80% completeness

	print "CL1040"
	z=0.7043
	sigma=418.
	errsigmap=55.
	errsigmam=46.

	prefix='cl1040'
	fullprefix='cl1040-1155'
	idra='1040'
	iddec='1155'
	mipsimage='/home/rfinn/research/clusters/spitzer/cl1040/mips/24/r13810432/ch1/bcd/pbcd/Combine/mosaic.fits'
	rac='160.1733'
	decc='-11.9308'
	cl1040 = ediscs24()
	cl1040.gotoitHa(catalogpath,'cl1040','1040','1155','cl1040.v1.spe','cl1040.id',mipsimage)
if (gocl105412 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.82#apcor for
	apcor3=1.99
	apcor2=2.58
	f80=75.#100.#flux in uJy at 80% completeness

	print "CL105412"
	z=.7498
	sigma=504.
	errsigmap=113.
	errsigmam=65.

	prefix='cl105412'
	fullprefix='cl1054-1245'
	idra='1054'
	iddec='1245'
	mipsimage='/home/rfinn/research/clusters/spitzer/cl105412/mips/24/r13812480/ch1/bcd/pbcd/Combine/mosaic.fits'
	rac='163.6813'
	decc='-12.7639'
	cl105412 = ediscs24()
	cl105412.gotoitHa(catalogpath,'cl105412','1054','1245','cl105412.v1.spe','cl1054-12.id',mipsimage)


if (gocl105411 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.86#apcor for
	apcor3=2.00
	apcor2=2.58
	f80=72.5#90.#flux in uJy at 80% completeness

	print "CL105411"
	z=.6972
	sigma=589.
	errsigmap=78.
	errsigmam=70.

	prefix='cl105411'
	fullprefix='cl1054-1146'
	idra='1054'
	iddec='1146'
	rvalue='r13810176'
	rac='163.6008'
	decc='-11.7717'

	cl105411 = ediscs24()
	cl105411.gotoit()

if (gocl1354 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.82#apcor for
	apcor3=1.97
	apcor2=2.55
	f80=82.5#90.#flux in uJy at 80% completeness

	print "CL1354"
	z=.7620
	sigma=648.
	errsigmap=105.
	errsigmam=110.

	prefix='cl1354'
	fullprefix='cl1354-1230'
	idra='1354'
	iddec='1230'
	rvalue='r13810944'
	rac='208.5396'
	decc='-12.5164'

	cl1354 = ediscs24()
	#(self,catalogpath,prefix,ra,dec,speccat,morphidcat,mipsimage):
	cl1354.gotoit()

if (gocl1037 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.81#apcor for
	apcor3=1.97
	apcor2=2.53
	f80=82.5#140.#flux in uJy at 80% completeness

	print "CL1037"
	z=.5783
	sigma=319.
	errsigmap=53.
	errsigmam=52.

	prefix='cl1037'
	fullprefix='cl1037-1243'
	idra='1037'
	iddec='1243'
	rvalue='r17053184'
	rac='159.46875'
	decc='-12.72916666667'
	cl1037 = ediscs24()
	cl1037.gotoit()

if (gocl1232 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.83#apcor for
	apcor3=1.99
	apcor2=2.57
	f80=97.5#90.#flux in uJy at 80% completeness

	print "CL1232"
	z=.5414
	sigma=1080.
	errsigmap=119.
	errsigmam=89.

	prefix='cl1232'
	fullprefix='cl1232-1250'
	idra='1232'
	iddec='1250'
	rvalue='r17053696'
	rac='188.125833333'
	decc='-12.8433333'

	cl1232 = ediscs24()
	cl1232.gotoit()



if (gocl1018 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.90#apcor for
	apcor3=2.01
	apcor2=2.54
	f80=80.#90.#flux in uJy at 80% completeness

	print "CL1018"
	z=.4734
	sigma=486.
	errsigmap=59.
	errsigmam=63.

	prefix='cl1018'
	fullprefix='cl1018-1211'
	idra='1018'
	iddec='1211'
	rvalue='r17783552'
	rac='154.69458333'
	decc='-12.19805555'

	cl1018 = ediscs24()
	cl1018.gotoit()


if (gocl1059 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.86#apcor for
	apcor3=2.00
	apcor2=2.55
	f80=82.5#90.#flux in uJy at 80% completeness

	print "CL1059"
	z=.4564
	sigma=510.
	errsigmap=52.
	errsigmam=56.

	prefix='cl1059'
	fullprefix='cl1059-1253'
	idra='1059'
	iddec='1253'
	rvalue='r17783808'
	rac='164.7895833333'
	decc='-12.88888888'
	cl1059 = ediscs24()
	cl1059.gotoit()
if (gocl1103 > .1):
	acsflag=1#cluster has acs data
	apcor4=1.91#apcor for
	apcor3=2.03
	apcor2=2.60
	f80=82.5#80.#flux in uJy at 80% completeness

	print "CL1103"
	z=.6261
	sigma=336.
	errsigmap=36.
	errsigmam=40.

	prefix='cl1103'
	fullprefix='cl1103-1245'
	idra='1103'
	iddec='1245'
	rvalue='r17785344'
	rac='165.92916666667'
	decc='-12.7602777777'
	cl1103 = ediscs24()
	cl1103.gotoit()
if (gocl1138 > .1):
	acsflag=1#cluster has acs data
	apcor4=1.82#apcor for
	apcor3=1.96
	apcor2=2.51
	f80=85.#90.#flux in uJy at 80% completeness

	print "CL1138"
	z=.4796
	sigma=732.
	errsigmap=72.
	errsigmam=76.

	prefix='cl1138'
	fullprefix='cl1138-1133'
	idra='1138'
	iddec='1133'
	rvalue='r17785600'
	rac='174.5416666666'
	decc='-11.560555555'
	cl1138 = ediscs24()
	cl1138.gotoit()
if (gocl1202 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.76#apcor for
	apcor3=1.93
	apcor2=2.52
	f80=97.#90.#flux in uJy at 80% completeness

	print "CL1202"
	z=.4240
	sigma=518.
	errsigmap=92.
	errsigmam=104.

	prefix='cl1202'
	fullprefix='cl1202-1224'
	idra='1202'
	iddec='1224'
	rvalue='r17784064'
	rac='180.679583333'
	decc='-12.4083333333'
	cl1202 = ediscs24()
	cl1202.gotoit()
if (gocl1301 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.81#apcor for
	apcor3=1.95
	apcor2=2.56
	f80=85.#100.#flux in uJy at 80% completeness

	print "CL1301"
	z=.4828
	sigma=687.
	errsigmap=81.
	errsigmam=86.

	prefix='cl1301'
	fullprefix='cl1301-1139'
	idra='1301'
	iddec='1139'
	rvalue='r17784320'
	rac='195.417083333'
	decc='-11.656388888'
	cl1301 = ediscs24()
	cl1301.gotoit()
if (gocl1353 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.87#apcor for
	apcor3=2.01
	apcor2=2.58
	f80=100.#flux in uJy at 80% completeness

	print "CL1353"
	z=.5882
	sigma=666.
	errsigmap=136.
	errsigmam=139.

	prefix='cl1353'
	fullprefix='cl1353-1137'
	idra='1353'
	iddec='1137'
	rvalue='r17784576'
	rac='208.255833333'
	decc='-11.624444444'
	cl1353 = ediscs24()
	cl1353.gotoit()
if (gocl1411 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.86#apcor for
	apcor3=1.97
	apcor2=2.49
	f80=97.5#100.#flux in uJy at 80% completeness

	print "CL1411"
	z=.5195
	sigma=710.
	errsigmap=125.
	errsigmam=133.

	prefix='cl1411'
	fullprefix='cl1411-1148'
	idra='1411'
	iddec='1148'
	rvalue='r17784832'
	rac='212.77000000'
	decc='-11.808055555'
	cl1411 = ediscs24()
	cl1411.gotoit()
if (gocl1420 > .1):
	acsflag=0#cluster has acs data
	apcor4=1.89#apcor for
	apcor3=1.99
	apcor2=2.55
	f80=95.#110.#flux in uJy at 80% completeness

	print "CL1420"
	z=.4962
	sigma=218.
	errsigmap=43.
	errsigmam=50.

	fullprefix='cl1420-1236'
	prefix='cl1420'
	idra='1420'
	iddec='1236'
	rvalue='r17785088'
	rac='215.0829166666'
	decc='-12.608333333'
	cl1420 = ediscs24()
	cl1420.gotoit()
if (gocl1227 > 0.1):
	acsflag=1#cluster has acs data
	apcor4=1.89#apcor for
	apcor3=2.01
	apcor2=2.54
	f80=80.#80.#flux in uJy at 80% completeness

	print "CL1227"
	z=.6357
	sigma=574.
	errsigmap=72.
	errsigmam=75.

	prefix='cl1227'
	fullprefix='cl1227-1138'
	idra='1227'
	iddec='1138'
	rvalue='r17053952'
	rac='186.9746'#12 27 58.9
	decc='-11.6389'#-11 35 13.5
	cl1227 = ediscs24()
	cl1227.gotoit()

morphout.close()

vandanaout.close()




#print "CL1216 - newspecmatchflag"
#print cl1216.newspecmatchflag
#for i in range(len(cl1216.newspecmatchflag)):
#	if abs(cl1216.newspecmatchflag[i]-0.) < 0.1:
#		print i,cl1216.newspecmatchflag[i]


c1=cl1216
c2=cl1354
c3=cl105412
c4=cl1040
c5=cl105411
c6=cl1227
c7=cl1353
c8=cl1037
c9=cl1232
c10=cl1411
c11=cl1420
c12=cl1301
c13=cl1138
c14=cl1018
c15=cl1059
c16=cl1202

calcmatchstats()


IRoptdist()

calcIRfraco2()

fieldlir=[]
fieldz=[]
fieldmatchflag24=[]
fieldMv=[]
fieldlir=fieldlir+list(c1.fieldLir)+list(c2.fieldLir)+list(c3.fieldLir)+list(c4.fieldLir)+list(c5.fieldLir)+list(c6.fieldLir)+list(c7.fieldLir)+list(c8.fieldLir)+list(c9.fieldLir)+list(c10.fieldLir)+list(c11.fieldLir)+list(c12.fieldLir)+list(c13.fieldLir)+list(c14.fieldLir)+list(c15.fieldLir)+list(c16.fieldLir)#+list(c17.fieldLir)
fieldz=fieldz+list(c1.fieldz)+list(c2.fieldz)+list(c3.fieldz)+list(c4.fieldz)+list(c5.fieldz)+list(c6.fieldz)+list(c7.fieldz)+list(c8.fieldz)+list(c9.fieldz)+list(c10.fieldz)+list(c11.fieldz)+list(c12.fieldz)+list(c13.fieldz)+list(c14.fieldz)+list(c15.fieldz)+list(c16.fieldz)#+list(c17.fieldz)
fieldmatchflag24=fieldmatchflag24+list(c1.fieldmatchflag24)+list(c2.fieldmatchflag24)+list(c3.fieldmatchflag24)+list(c4.fieldmatchflag24)+list(c5.fieldmatchflag24)+list(c6.fieldmatchflag24)+list(c7.fieldmatchflag24)+list(c8.fieldmatchflag24)+list(c9.fieldmatchflag24)+list(c10.fieldmatchflag24)+list(c11.fieldmatchflag24)+list(c12.fieldmatchflag24)+list(c13.fieldmatchflag24)+list(c14.fieldmatchflag24)+list(c15.fieldmatchflag24)+list(c16.fieldmatchflag24)#+list(c17.fieldmatchflag24)

fieldMv=fieldMv+list(c1.fieldMv)+list(c2.fieldMv)+list(c3.fieldMv)+list(c4.fieldMv)+list(c5.fieldMv)+list(c6.fieldMv)+list(c7.fieldMv)+list(c8.fieldMv)+list(c9.fieldMv)+list(c10.fieldMv)+list(c11.fieldMv)+list(c12.fieldMv)+list(c13.fieldMv)+list(c14.fieldMv)+list(c15.fieldMv)+list(c16.fieldMv)#+list(c17.fieldMv)

fieldlir=N.array(fieldlir,'d')
fieldz=N.array(fieldz,'d')
fieldMv=N.array(fieldMv,'d')
fieldmatchflag24=N.array(fieldmatchflag24,'d')
field=N.average(N.compress(fieldmatchflag24>.1,fieldlir))
print "average Lir of EDisCS field galaxies = ",field

avefieldz=N.average(N.compress(fieldmatchflag24> .1,fieldz))
print "average z of field galaxies = ",avefieldz,"+/-",pylab.std(fieldz)

fieldMvcut=calcMvcut(fieldz)

fieldlirloz=N.compress((fieldz < 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldlir)
fieldlirhiz=N.compress((fieldz > 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldlir)
fieldzloz=N.compress((fieldz < 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldz)
fieldzhiz=N.compress((fieldz > 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldz)

nfieldlirloz=len(N.compress((fieldz < 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldlir))
nfieldlirhiz=len(N.compress((fieldz > 0.6) & (fieldmatchflag24 > .1) & (fieldMv < fieldMvcut),fieldlir))

ntotfieldloz=len(N.compress((fieldz < 0.6) & (fieldMv < fieldMvcut),fieldlir))
ntotfieldhiz=len(N.compress((fieldz > 0.6) & (fieldMv < fieldMvcut),fieldlir))


fieldlirlozave=N.average(fieldlirloz)
fieldlirlozstd=pylab.std(fieldlirloz)/N.sqrt(len(fieldlirloz))

fieldlirhizave=N.average(fieldlirhiz)
fieldlirhizstd=pylab.std(fieldlirhiz)/N.sqrt(len(fieldlirhiz))


print "Number of field galaxies w/z < 0.6 = %i, n/ntot = %5.2f, ave z = %5.2f, ave lir = %5.2e "%(len(fieldlirloz),1.*len(fieldlirloz)/(1.*len(fieldlir)),N.average(fieldzloz),N.average(fieldlirloz))
print "Number of field galaxies w/z > 0.6 = %i, n/ntot = %5.2f, ave z = %5.2f, ave lir = %5.2e "%(len(fieldlirhiz),1.*len(fieldlirhiz)/(1.*len(fieldlir)),N.average(fieldzhiz),N.average(fieldlirhiz))
	
print "total number of field galaxies (w/ and w/out 24um) = ",len(fieldz)


def sffzvcolor():#plot IR frac vs z, separating blue and red galaxies
	blueIRfrac=[]
	redIRfrac=[]
	z=[]

	for i in range(ncl):
		if i == 0:
			cl=c1
		if i == 1:
			cl=c2
		if i == 2:
			cl=c3
		if i == 3:
			cl=c4
		if i == 4:
			cl=c5
		if i == 5:
			cl=c6
		if i == 6:
			cl=c7
		if i == 7:
			cl=c8
		if i == 8:
			cl=c9
		if i == 9:
			cl=c10
		if i == 10:
			cl=c11
		if i == 11:
			cl=c12
		if i == 12:
			cl=c13
		if i == 13:
			cl=c14
		if i == 14:
			cl=c15
		if i == 15:
			cl=c16
		if i == 16:
			cl=c17


writeapj()

getgroups()

makeplots=0
plotlirMrall()
plotcolormagallUV()
plotssfrstellmass()
plotfraclirMv()
plotfraclirstellmass()
plotlirstellmass()
#skipping for now to make catalogs for Pascale
#plotlirMvredblue()


#plotcolormagsalt()
#plotcolormagallBV()

if makeplots > 0.1:

##########    plots   ############
#SFRHaIR()
	plotbvdr()
	plotlirdistall()
	plotlirdistcomb()

	plotlirdistz()
	plotlirpositionsall()
	plotsffdrcomb()
	plotcolormagall()
	plotsigmaz()	

	
#  no longer using these next plots
#plotcolormagall()
#plotLIRLI()
	plotcolormagmorphall()

	plotveldistall()



# calculate fraction of 24 um sources with multiple optical matches
nmultimatch=c1.nmultimatch+c2.nmultimatch+c3.nmultimatch+c4.nmultimatch+c5.nmultimatch+c6.nmultimatch+c7.nmultimatch+c8.nmultimatch+c9.nmultimatch+c10.nmultimatch+c11.nmultimatch+c12.nmultimatch+c13.nmultimatch+c14.nmultimatch+c15.nmultimatch+c16.nmultimatch#+c17.nmultimatch

n24match=c1.n24match+c2.n24match+c3.n24match+c4.n24match+c5.n24match+c6.n24match+c7.n24match+c8.n24match+c9.n24match+c10.n24match+c11.n24match+c12.n24match+c13.n24match+c14.n24match+c15.n24match+c16.n24match#+c17.n24match

(a,b,c)=my.ratioerror(nmultimatch,n24match)
print "Fraction of 24um sources w/>1 optical match = %5.4f - %5.4f + %5.4f (%i/%i)"%(a,b,c,nmultimatch,n24match)




#ediscsID,ra,dec,xcorr,ycorr,starflag,EWha,EWhaerr,SFR,SFRerr,matchflaghaediscs,SF_flag, hstgimtype,matchflagmorphgimediscs,hstvisnumtype,matchflagmorphvisediscs,matchflag24,f24,errf24,nmatchediscs24,misoV,misoeVapsim,misoR,misoeRapsim,misoI,misoeIapsim,misoJ,misoeJapsim,misoK,misoeKapsim, magV,mageVapsim,magR,mageRapsim,magI,mageIapsim,magJ,mageJapsim,magK,mageKapsim, membflag,newspecmatchflag,defmembflag,specz,spectype,specEWOII,matchflagspecediscs,specEWOIIflag,bestz, lowz, highz, wmin, Pclust,LUlowzclust,LUbestzclust,LUhighzclust,LBlowzclust,LBbestzclust,LBhighzclust,LVlowzclust,LVbestzclust,LVhighzclust,LRlowzclust,LRbestzclust,LRhighzclust,LIlowzclust,LIbestzclust,LIhighzclust,LJlowzclust,LJbestzclust,LJhighzclust,LKlowzclust,LKbestzclust,LKhighzclust,fluxK,UBlowzclust,UBbestzclust,UBhighzclust,BVlowzclust,BVbestzclust,BVhighzclust,UVlowzclust,UVbestzclust,UVhighzclust,matchflagediscsirac,iracf1,iracf2,iracf3,iracf4,erriracf1,erriracf2,erriracf3,erriracf4,iracsexflag0,iracsexflag1,iracwch1,iracwch2,iracwch3,iracwch4,iracwmin,nmatchediscsirac,l24,errl24,lha,lhaerr,snr24,imagex24,imagey24,fap1,fap2,fap3,fap4,fap5,fap6,fap7,fap8,fap9,fap10,errfap1,errfap2,errfap3,errfap4,errfap5,errfap6,errfap7,errfap8,errfap9,errfap10

#o2stats()
try:
	print 'Average contam of galaxies w/in 24um aperture = ',N.average(blendall)
	print 'Average contam of 24um galaxies w/in 24um aperture = ',N.average(blendall24)
except:
	print "didn't measure contamination this time"

