#!/usr/bin/env python
""" program to look at H-alpha properties as a function of local environment 

"""
#import sys
#sys.path.append('/Users/rfinn/bin')

#import z08readfits
from pylab import *
import os
import pyfits
import atpy
from mpl_toolkits.mplot3d import Axes3D
import mystuff as my
from z08common import *

h100=.7
zphotmin=0.7
zphotmax=0.9
zmin=.78
zmax=0.82

def findnearest(x1,y1,x2,y2,delta):#use where command
	matchflag=1
	nmatch=0
	d=sqrt((x1-x2)**2 + (y1-y2)**2)#x2 and y2 are arrays
	t=d[d<delta]
	matches=t
	if len(matches) > 0:
		nmatch=len(matches)
		if nmatch > 1:
			try:
				z=[]
				for im in matches:
					z.append(d[int(im)])
				z=array(z,'f')
				imatch=matches[(z == z.min())]
			except IndexError:
				print len(d),matches,array(matches,'i'),d[matches[0]],d[matches[1]],d[matches],d[array(matches,'i')]
		else:
			imatch=matches[0]
			
	else:
		imatch = 0
		matchflag = 0


	return imatch, matchflag,nmatch


def biniterr(x,y,n):#bin arrays x, y into n equally-populated bins, returning xbin,ybin
    nx=len(x)
    y=take(y,argsort(x))
    x=take(x,argsort(x))
    xbin=zeros(n,'f')
    xbin=zeros(n,'f')
    ybin=zeros(n,'f')
    ybinerr=zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        #xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        #ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        xbin[i]=average(x[nmin:nmax])
        ybin[i]=average(y[nmin:nmax])
        ybinerr[i]=std(y[nmin:nmax])/sqrt(1.*(nmax-nmin))
    return xbin, ybin, ybinerr

class cluster:
    def __init__(self,infile):
	self.name=name
	in1=open(infile,'r')
	ngal=0
	for line in in1:
	    if line.find('#') > -1:
		continue
	    ngal += 1
	in1.close()

	self.J=zeros(ngal,'f')
	self.EW=zeros(ngal,'f')
	self.EWerr=zeros(ngal,'f')
	self.sfflag=zeros(ngal,'f')

	in1=open(infile,'r')
	i=0
	for line in in1:
	    if line.find('#') > -1:
		continue

	    t=line.split()
	    self.J[i]=float(t[7])
	    self.EW[i]=float(t[12])
	    self.EWerr[i]=float(t[13])
	    self.sfflag[i]=int(t[16])
	    i += 1
	in1.close()

	self.Jm=self.J[(self.sfflag>0.1)]
	self.EWm=self.EW[(self.sfflag>0.1)]
	self.EWerrm=self.EWerr[(self.sfflag>0.1)]

class newfirmMedBandSurvey:
    def __init__(self):
	    print 'initializing nmbs data'
	    self.readphotcat()
	    self.readzcat()
    def readphotcat(self):
# id  x  y  ra  dec Kaper eKaper CH4 eCH4 wCH4 CH3 eCH3 wCH3 CH2 eCH2 wCH2 CH1 eCH1 wCH1 Ks eKs wKs K eK wK H eH wH H2 eH2 wH2 H1 eH1 wH1 J eJ wJ J3 eJ3 wJ3 J2 eJ2 wJ2 J1 eJ1 wJ1 Zp eZp wZp z ez wz Ip eIp wIp I eI wI Rp eRp wRp R eR wR G eG wG V eV wV B eB wB U eU wU NUV eNUV wNUV FUV eFUV wFUV IA427 eIA427 IA464 eIA464 IA484 eIA484 IA505 eIA505 IA527 eIA527 IA574 eIA574 IA624 eIA624 IA679 eIA679 IA709 eIA709 IA738 eIA738 IA767 eIA767 IA827 eIA827 ftot24um_uJy f24um_uJy e24um_uJy w24um wmin wmin_irac z_spec star_flag ap_col ap_tot totcor K_ellip K_theta_J2000 K_R50 K_class_star K_flags UH2_flags Near_Star CH1_contam CH2_contam CH3_contam CH4_contam NUV_contam FUV_contam contam_flag nchild id_parent use
	data=mlab.load('/Users/rfinn/research/FieldSurvey/nmbs/cosmos-1.deblend.v5.1.cat')
	self.id=array(data[:,0],'f')
	self.ra=array(data[:,3],'f')
	self.dec=array(data[:,4],'f')
    def readzcat(self):
	data=mlab.load('/Users/rfinn/research/FieldSurvey/nmbs/cosmos-1.deblend.redshifts/cosmos-1.deblend.v5.1.zout')
	self.zid=array(data[:,0],'f')
	self.zspec=array(data[:,1],'f')
	self.za=array(data[:,2],'f')
	self.zm1=array(data[:,3],'f')
	self.chia=array(data[:,4],'f')
	self.zp=array(data[:,5],'f')
	self.chip=array(data[:,6],'f')
	self.zm2=array(data[:,7],'f')
	self.odds=array(data[:,8],'f')
	self.l68=array(data[:,9],'f')
	self.u68=array(data[:,10],'f')
	self.l95=array(data[:,11],'f')
	self.u95=array(data[:,12],'f')
	self.l99=array(data[:,13],'f')
	self.u99=array(data[:,14],'f')
	self.nfilt=array(data[:,15],'f')
	self.qz=array(data[:,16],'f')
	self.zpeak=array(data[:,17],'f')
	self.peakprob=array(data[:,18],'f')
	self.zmc=array(data[:,19],'f')
	self.fieldflag=(self.zp > zmin) & (self.zp < zmax)

class primus:
    def __init__(self):
	zcat=homedir+'research/PRIMUS/PRIMUS_2013_zcat_v1_topcat.fits'
	self.p=atpy.Table(zcat,type='fits')
	ramin_xmm=33.
	ramax_xmm=36.4
	decmin_xmm=-6.3
	decmax_xmm=-4.4
	self.zflag=(self.p.Z > zmin) & (self.p.Z < zmax)
	self.xmmflag=(self.p.RA > ramin_xmm) & (self.p.RA < ramax_xmm) & (self.p.DEC > decmin_xmm) & (self.p.DEC < decmax_xmm) & self.zflag
	ramin_cos=149.2
	ramax_cos=151.
	decmin_cos=1.6
	decmax_cos=3.
	self.cosmosflag=(self.p.RA > ramin_cos) & (self.p.RA < ramax_cos) & (self.p.DEC > decmin_cos) & (self.p.DEC < decmax_cos) & self.zflag

class field:
    def __init__(self,infile,infile2):
	self.name=name
	self.decmin=decmin
	self.decmax=decmax
	self.ramin=ramin
	self.ramax=ramax
	""" 
	['ID',
	'DET',
	'SE_X',
	'SE_Y',
	'SE_RA',
	'SE_DEC',
	'NB1190_APER',
	'NB1190_APER_E',
	'J_APER',
	'J_APER_E',
	'J-NB1190_APER',
	'J-NB1190_APER_E',
	'S/N',
	'EW',
	'EW_ERR',
	'F_LINE',
	'F_LINE_ERR',
	'NB1190_SE_FLAG',
	'J_SE_FLAG',
	'NB1190_BEST',
	'NB1190_BEST_E',
	'J_BEST',
	'J_BEST_E',
	'NB1190_EXPTIME',
	'J_EXPTIME',
	'ZPHOT-ID',
	'ZPHOT',
	'RADIFF1',
	'DECDIFF1',
	'B',
	'V',
	'R',
	'I',
	'Z',
	'ZSPEC-ID',
	'ZSPEC',
	'RADIFF2',
	'DECDIFF2',
	'QUAL',
	'MAG',
	'ZCOMMENT',
	'NB1190_IM',
	'J_IM']
	"""
	data, hdr = pyfits.getdata(infile,1,header=True)
	self.ID=data.field('ID')
	self.DET=data.field('DET')
	self.x=data.field('SE_X')
	self.y=data.field('SE_Y')
	self.ra=data.field('SE_RA')
	self.dec=data.field('SE_DEC')
	self.NB1190APER=data.field('NB1190_APER')
	self.NB1190APER_E=data.field('NB1190_APER_E')
	self.JAPER=data.field('J_APER')
	self.JAPERerr=data.field('J_APER_E')
	self.JNB1190APER=data.field('J_NB1190_APER')
	self.JNB1190APERE=data.field('J_NB1190_APER_E')
	self.SNR=data.field('SNR')
	self.EW=data.field('EW')
	self.EWerr=data.field('EW_ERR')
	self.FLINE=data.field('F_LINE')
	self.FLINEerr=data.field('F_LINE_ERR')
	self.NB1190SEFLAG=data.field('NB1190_SE_FLAG')
	self.JSEFLAG=data.field('J_SE_FLAG')
	self.NB1190BEST=data.field('NB1190_BEST')
	self.NB1190BEST_E=data.field('NB1190_BEST_E')
	self.JBEST=data.field('J_BEST')
	self.JBESTE=data.field('J_BEST_E')
	self.nbexptime=data.field('NB1190_EXPTIME')
	self.jexptime=data.field('J_EXPTIME')
	self.zphotid=data.field('ZPHOT_ID')
	self.zphot=data.field('ZPHOT')
	self.radiff1=data.field('RADIFF1')
	self.decdiff1=data.field('DECDIFF1')
	self.B=data.field('B')
	self.V=data.field('V')
	self.R=data.field('R')
	self.I=data.field('I')
	self.Z=data.field('Z')
	self.zspecid=data.field('ZSPEC_ID')
	self.zspec=data.field('ZSPEC')
	#self.zspec=zeros(len(temp),'f')
	self.zspecCode=[]
	#for i in range(len(temp)):
	#	self.zspec[i]=temp[i][0:-1]
	#	self.zspecCode.append(temp[i][6])
	
	self.radiff2=data.field('RADIFF2')
	self.decdiff2=data.field('DECDIFF2')
	self.qual=data.field('QUAL')
	self.mag=data.field('MAG')
	self.zcomment=data.field('ZCOMMENT')
	self.nb1190im=data.field('NB1190_IM')
	self.Jim=data.field('J_IM')

	membflag=zeros(len(self.zspec),'i')
	photmembflag=zeros(len(self.zspec),'i')
	supermembflag=zeros(len(self.zspec),'i')#uses spec and phot-z information

	i=0
	z=[]
	for zs in self.zspec:
	    t=float(zs[0:len(zs)-1])
	    z.append(t)
	    if (t > zmin) & (t < zmax):
		membflag[i]=1
	    if t > 9.:#no spectroscopy
		membflag[i]=-1	    
	    t=float(self.zphot[i])
	    if (t > zphotmin) & (t < zphotmax):
		photmembflag[i]=1
	    i += 1
	self.z=array(z)

	self.membflag=array(membflag)
	self.photmembflag=array(photmembflag)


	for i in range(len(membflag)):
	    if membflag[i] > -0.1:
		    if membflag[i] > .1:
			    supermembflag[i] = 1
		    continue
	    else:
		supermembflag[i]=photmembflag[i]


	#green flag
	vi=self.V-self.I
	i=self.I
	dm=3.3
	ymin=-.09*(i-20)+dm-1.6
	
	dm=3.6
	ymax=-.09*(i-20)+dm-1.6
	self.greenflag=(vi > ymin) &  (vi < ymax)
	

	self.ra=array(self.ra)
	self.dec=array(self.dec)
	self.x=array(self.x)
	self.y=array(self.y)
	self.EW=array(self.EW)
	self.EWerr=array(self.EWerr)
	self.JBEST=array(self.JBEST)
	self.JBESTE=array(self.JBESTE)
	self.R=array(self.R)
	#cull members
	self.ram=array(self.ra[membflag>0])
	self.decm=array(self.dec[membflag>0])
	self.xm=array(self.x[membflag>0])
	self.ym=array(self.y[membflag>0])
	self.zm=array(self.z[membflag>0])
	self.EWm=array(self.EW[membflag>0])/(1+zfield)
	self.EWerrm=array(self.EWerr[membflag>0])/(1+zfield)
	self.Jm=array(self.JBEST[membflag>0])
	self.Jerrm=array(self.JBESTE[membflag>0])
	#print len(self.R),len(membflag)
	self.Rm=array(self.R[membflag>0])

	self.Vm=array(self.V[membflag>0])
	self.Im=array(self.I[membflag>0])
	self.greenflagm=array(self.greenflag[membflag>0])
	#cull zphot members w/Halpha emission
	self.ramp=array(self.ra[supermembflag>0])
	self.decmp=array(self.dec[supermembflag>0])
	self.xmp=array(self.x[supermembflag>0])
	self.ymp=array(self.y[supermembflag>0])
	self.zmp=array(self.z[supermembflag>0])
	self.EWmp=array(self.EW[supermembflag>0])/(1+zfield)
	self.EWerrmp=array(self.EWerr[supermembflag>0])/(1+zfield)
	self.Jmp=array(self.JBEST[supermembflag>0])
	self.Jerrmp=array(self.JBESTE[supermembflag>0])

	self.Rmp=array(self.R[supermembflag>0])

	self.Vmp=array(self.V[supermembflag>0])
	self.Imp=array(self.I[supermembflag>0])
	self.greenflagmp=array(self.greenflag[supermembflag>0])

	fname='/Users/rfinn/research/FieldSurvey/LocalDensity/'+self.name+'coordsz.txt'

	outfile=open(fname,'w')
	for i in range(len(self.ram)):
		s="%12.8e %12.8e  %12.8e %i \n"%(self.ram[i],self.decm[i],self.zm[i],self.greenflagm[i])
		#print i,s
		outfile.write(s)
	outfile.close()



	print 'in init', infile2
	self.getphotozcat(infile2)
	self.getspeczcat()
	#self.getnearest()
	self.getnearestall()
	self.getnearestallp()
	self.getnearestp()
	
    def plot3d(self):
	    fig = plt.figure()
	    ax = fig.add_subplot(111, projection='3d')
	    flag=self.membflag == 1
	    sp=ax.scatter(self.ra[flag], self.dec[flag], self.z[flag], c=log10(self.EW[flag]),vmin=log10(40),vmax=log10(500))
	    cb=colorbar(sp)
	    ax.set_xlabel('RA')
	    ax.set_ylabel('Dec')
	    ax.set_zlabel('Redshift')
	    title(self.name)


    def plotpositionszslice(self):
	    dz=.0025
	    zmin=.790
	    zmax=.82
	    zrange=arange(zmin,zmax,dz)
	    colors=['#000099','#110088','#220077','#330066','#440055','#550044','#660033','#770022','#880011','#990000']
	    colors=colors+colors
	    figure()
	    for i in range(len(zrange)):
		    zmin=zrange[i]
		    zmax=zmin+dz
		    flag=(self.zm<zmax) & (self.zm > zmin)
		    #print flag
		    plot(self.ram[flag],self.decm[flag],ls='None',marker='o',color=colors[i])
	    

    def getphotozcat(self,infile):

	#print 'in getphotozcat',infile
	in1=open(infile,'r')
	ngal=0
	for line in in1:
	    if line.find('#') > -1:
		continue
	    ngal += 1
	in1.close()
	#all variables start with zphot
	self.zphotz=zeros(ngal,'f')
	self.zphotRA=zeros(ngal,'f')
	self.zphotDec=zeros(ngal,'f')
	self.zphotB=zeros(ngal,'f')
	self.zphotV=zeros(ngal,'f')
	self.zphotR=zeros(ngal,'f')
	self.zphoti=zeros(ngal,'f')
	self.zphotzmag=zeros(ngal,'f')

	#print 'opening for second time ', infile
	in1=open(infile,'r')
	ngal=0
	for line in in1:
	    if line.find('#') > -1:
		continue
	    #print line,infile
	    t=line.split()
	    #print ngal,t
	    self.zphotz[ngal]=float(t[1])
	    self.zphotRA[ngal]=float(t[2])
	    self.zphotDec[ngal]=float(t[3])
	    self.zphotB[ngal]=float(t[4])
	    self.zphotV[ngal]=float(t[5])
	    self.zphotR[ngal]=float(t[6])
	    self.zphoti[ngal]=float(t[7])
	    self.zphotzmag[ngal]=float(t[8])
	    ngal += 1
	in1.close()


	zmin=.75
	zmax=0.85
	i=0
	z=[]
	zphotmembflag=zeros(len(self.zphotz),'i')
	i=0
	for z in self.zphotz:
	    if (z > zmin) & (z < zmax):
		zphotmembflag[i]=1
	    i += 1
	self.zphotRAm=self.zphotRA[zphotmembflag > .1]
	self.zphotDecm=self.zphotDec[zphotmembflag > .1]
	self.zphotBm=self.zphotB[zphotmembflag > .1]
	self.zphotVm=self.zphotV[zphotmembflag > .1]
	self.zphotRm=self.zphotR[zphotmembflag > .1]
	self.zphotim=self.zphoti[zphotmembflag > .1]
	self.zphotzmagm=self.zphotzmag[zphotmembflag > .1]

	#find photoz members that are not emission line objects

	Halphamatchflag=zeros(len(self.zphotRAm))
	Halphamatchindex=zeros(len(self.zphotRAm))
	delta=.5/3600.#offset for matching spec and Halpha catalog
	for i in range(len(self.zphotRAm)):
	    
	    imatch, matchflag,nmatch=findnearest(self.zphotRAm[i],self.zphotDecm[i],self.ram,self.decm,delta)	    
	    if matchflag > .1:
		Halphamatchflag[i]=matchflag
		Halphamatchindex[i]=imatch
		#print 'Number of matches = ',nmatch
	    #else:
		#print "spec member w/no Halpha"

		


	#  !!!!! FIX THIS  !!!!!

	#photoz sources w/out Halpha emission
	self.RAnoHp=self.zphotRAm[(Halphamatchflag < .1) & (self.zphotDecm > self.decmin) & (self.zphotDecm < self.decmax)]
	self.DecnoHp=self.zphotDecm[(Halphamatchflag < .1) & (self.zphotDecm > self.decmin)  & (self.zphotDecm < self.decmax)]
	#match to photometry to get mags
	self.EWnoHp=zeros(len(self.RAnoHp),'f')
	self.RAallp=concatenate((self.ramp,self.RAnoHp))
	self.Decallp=concatenate((self.decmp,self.DecnoHp))
	self.EWallp=concatenate((self.EWmp,self.EWnoHp))

	#check to make sure sources are w/in overlap region of photoz and Halpha
	# -5.49 < Dec < -5.025
	# 34.17 < RA < 34.65
	# Dec < -5.302 AND RA < 34.197

	if self.name.find('SXDS') > -1:
		flag1=((self.Decallp > -5.49)& (self.Decallp < -5.025))
		flag2=((self.RAallp > 34.17)& (self.RAallp < 34.65))
		notflag3=((self.Decallp < -5.302)& (self.RAallp < 34.197))
		self.posflagallp=flag1 & flag2 & ~notflag3
	#ignoring this for now


	flag1=((self.Decallp > self.decmin)& (self.Decallp < self.decmax))
	flag2=((self.RAallp > self.ramin)& (self.RAallp < self.ramax))
	self.posflagallp=flag1 & flag2 

	self.RAallp=self.RAallp[self.posflagallp]
	self.Decallp=self.Decallp[self.posflagallp]
	self.EWallp=self.EWallp[self.posflagallp]
	self.photozHamatchflag=Halphamatchflag[self.posflagallp]


    def getspeczcat(self):
	
	if self.name.find('SX') > -1:
	    infile='/Users/rfinn/research/FieldSurvey/catalogs/SXDS_zspec.cat'
	if self.name.find('COS') > -1:
	    infile='/Users/rfinn/research/FieldSurvey/catalogs/COSMOS_zspec.cat'
	#print 'in getphotozcat',infile
	in1=open(infile,'r')
	ngal=0
	for line in in1:
	    if line.find('#') > -1:
		continue
	    ngal += 1
	in1.close()
	#all variables start with zphot
	self.zspecz=zeros(ngal,'f')
	self.zspecRA=zeros(ngal,'f')
	self.zspecDec=zeros(ngal,'f')
	self.zspecCode=[]


	#print 'opening for second time ', infile
	in1=open(infile,'r')
	ngal=0
	for line in in1:
	    if line.find('#') > -1:
		continue
	    #print line,infile
	    t=line.split()
	    #print ngal,t
	    z=t[4]
	    self.zspecCode.append(z[len(z)-1])
	    z=z[0:(len(z)-1)]#remove last character which is code
	    #print z
	    self.zspecz[ngal]=float(z)
	    self.zspecRA[ngal]=float(t[1])
	    self.zspecDec[ngal]=float(t[2])
	    #self.zphoti[ngal]=float(t[7])
	    ngal += 1
	in1.close()


	zmin=.78
	zmax=0.82
	i=0
	z=[]
	membflag=zeros(len(self.zspecz),'i')
	i=0
	for z in self.zspecz:
	    if (z > zmin) & (z < zmax):
		membflag[i]=1
	    i += 1
	self.zspecRAm=self.zspecRA[membflag > .1]#m = members
	self.zspecDecm=self.zspecDec[membflag > .1]
	print "number of spec members from spec catalog = ",len(self.zspecRAm)

	#match spec members with Halpha catalog
	#create list of spec members that are not in Halpha catalog
	#need to create member list that includes both sets
	Halphamatchflag=zeros(len(self.zspecRAm))
	Halphamatchindex=zeros(len(self.zspecRAm))
	delta=1./3600.#offset for matching spec and Halpha catalog
	for i in range(len(self.zspecRAm)):
	    
	    imatch, matchflag,nmatch=findnearest(self.zspecRAm[i],self.zspecDecm[i],self.ram,self.decm,delta)	    
	    if matchflag > .1:
		Halphamatchflag[i]=matchflag
		Halphamatchindex[i]=imatch
		#print 'Number of matches = ',nmatch
	    #else:
		#print "spec member w/no Halpha"
	#extras are 
	self.RAnoH=self.zspecRAm[(Halphamatchflag < .1) & (self.zspecDecm > -5.5) & (self.zspecDecm < -5.)]
	self.DecnoH=self.zspecDecm[(Halphamatchflag < .1) & (self.zspecDecm > -5.5)  & (self.zspecDecm < -5.)]
	#match to photometry to get mags
	self.EWnoH=zeros(len(self.RAnoH),'f')

	self.RAall=concatenate((self.ram,self.RAnoH))
	self.Decall=concatenate((self.decm,self.DecnoH))
	self.EWall=concatenate((self.EWm,self.EWnoH))

    def getnearestall(self):
	#calculate nearest neighbor distances from all z=0.8 spec galaxies, including non-Ha galaxies
	self.sig5=zeros(len(self.RAall),'f')
        self.sig10=zeros(len(self.RAall),'f')
        self.nearest=zeros(len(self.RAall),'f')#distance to nearest neighbor
        self.dmagnearest=zeros(len(self.RAall),'f')#mag diff b/w spec obj and nearest


        for i in range(len(self.RAall)):
	    dspec=sqrt((self.RAall[i]-self.RAall)**2+(self.Decall[i]-self.Decall)**2)#sorted array of distances in degrees

	    dspecsort=take(dspec,argsort(dspec))

	    self.sig5[i]=5./(pi)/(sum(dspecsort[1:6])*3600.*angdist/1000.)**2
	    if (len(dspecsort) > 10):
		self.sig10[i]=10./(pi)/(sum(dspecsort[1:11])*3600.*angdist/1000.)**2
	    else:
		self.sig10[i]=0.

	    self.nearest[i]=dspecsort[1]*3600.*angdist#first element of array is dist from galaxy to itself

    def getnearestallp(self):
	#calculate nearest neighbor distances from all z=0.8 photz galaxies, including non-Ha galaxies
	self.sig5allp=zeros(len(self.RAallp),'f')
        self.sig10allp=zeros(len(self.RAallp),'f')
        self.nearestallp=zeros(len(self.RAallp),'f')#distance to nearest neighbor
        self.dmagnearestallp=zeros(len(self.RAallp),'f')#mag diff b/w spec obj and nearest


        for i in range(len(self.RAallp)):
	    #dspec=sqrt((self.RAallp[i]-self.RAallp)**2+(self.Decallp[i]-self.Decallp)**2)#sorted array of distances in degrees
	    dspec=sqrt((self.RAallp[i]-self.zphotRAm)**2+(self.Decallp[i]-self.zphotDecm)**2)#sorted array of distances in degrees

	    dspecsort=take(dspec,argsort(dspec))

	    self.sig5allp[i]=5./(pi)/(sum(dspecsort[1:6])*3600.*angdist/1000.)**2
	    if (len(dspecsort) > 10):
		self.sig10allp[i]=10./(pi)/(sum(dspecsort[1:11])*3600.*angdist/1000.)**2
	    else:
		self.sig10allp[i]=0.

	    self.nearestallp[i]=dspecsort[1]*3600.*angdist#first element of array is dist from galaxy to itself

    def getnearest(self):
	self.sig5=zeros(len(self.ram),'f')
        self.sig10=zeros(len(self.ram),'f')
        self.nearest=zeros(len(self.ram),'f')#distance to nearest neighbor
        self.dmagnearest=zeros(len(self.ram),'f')#mag diff b/w spec obj and nearest


        for i in range(len(self.ram)):
	    dspec=sqrt((self.ram[i]-self.zspecRAm)**2+(self.decm[i]-self.zspecDecm)**2)#sorted array of distances in degrees
	    #Jsort=take(self.zphotim,argsort(dspec))
	    dspecsort=take(dspec,argsort(dspec))
	    #self.sig5[i]=5./(pi)/(dspecsort[5]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
	    self.sig5[i]=5./(pi)/(sum(dspecsort[1:6])*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
	    if (len(dspecsort) > 10):
		self.sig10[i]=10./(pi)/(sum(dspecsort[1:11])*3600.*angdist/1000.)**2
		#self.sig10[i]=10./(pi)/(dspecsort[10]*3600.*angdist/1000.)**2
	    else:
		self.sig10[i]=0.

	    self.nearest[i]=dspecsort[1]*3600.*angdist#first element of array is dist from galaxy to itself
	    #self.dmagnearest[i]=Jsort[1]-Jsort[0]

    def getnearestp(self):
	self.sig5p=zeros(len(self.ramp),'f')
        self.sig10p=zeros(len(self.ramp),'f')
        self.nearestp=zeros(len(self.ramp),'f')#distance to nearest neighbor
        self.dmagnearestp=zeros(len(self.ramp),'f')#mag diff b/w spec obj and nearest


        for i in range(len(self.ramp)):
	    dspec=sqrt((self.ramp[i]-self.zphotRAm)**2+(self.decmp[i]-self.zphotDecm)**2)#sorted array of distances in degrees
	    Jsort=take(self.zphotim,argsort(dspec))
	    dspecsort=take(dspec,argsort(dspec))
	    #self.sig5p[i]=5./(pi)/(dspecsort[5]*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
	    self.sig5p[i]=5./(pi)/(sum(dspecsort[1:4])*3600.*angdist/1000.)**2#convert from deg to arcsec, multiply by DA (kpc/arcsec), divide by 1000 to convert to Mpc, index 5 element b/c 0 is itself
	    if (len(dspecsort) > 10):
		self.sig10p[i]=10./(pi)/(sum(dspecsort[1:11])*3600.*angdist/1000.)**2
	    else:
		self.sig10p[i]=0.

	    self.nearestp[i]=dspecsort[1]*3600.*angdist#first element of array is dist from galaxy to itself
	    self.dmagnearestp[i]=Jsort[1]-Jsort[0]


    def plotfig1(self):
	figure(1)
	clf()
	ra=array(self.ra)
	dec=array(self.dec)
	plot(ra[self.membflag>0],dec[self.membflag>0],'ko',markersize=3)
	#volume=(15.*self.EWm/max(self.EWm))**2
	#scatter(self.ram,self.decm,s=volume,alpha=0.75)
	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)
	xmin,xmax=xlim()
	xlim(xmax,xmin)

    def plotfig2(self):
	figure(2)
	clf()
	hist(self.z[self.membflag > 0])

    def plotprimuscoverage(self):
	figure()
	clf()
	plot(self.ramp[self.greenflagmp],self.decmp[self.greenflagmp],'co',markersize=18,label='Green')
	plot(self.ramp,self.decmp,'bo',markersize=10, label='Halpha+specz')
	xmin,xmax=xlim()
	ymin,ymax=ylim()


	plot(self.ramp,self.decmp,'go',markersize=8,label='Halpha+photz')
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=5,label='Specz')
	plot(self.ram,self.decm,'ro',markersize=5,label='Specz')
	#plot(self.zphotRAm,self.zphotDecm,'ko',markersize=2,label='Photz')


	if self.name.find('SXDS') > -1:
		plot(prim.p.RA[prim.p.XMMFLAG],prim.p.DEC[prim.p.XMMFLAG],'ko',markersize=2,label='PRIMUS')
	if self.name.find('COSMOS') > -1:
		plot(prim.p.RA[prim.p.COSMOSFLAG],prim.p.DEC[prim.p.COSMOSFLAG],'ko',markersize=2,label='PRIMUS')



	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)


	#plot(self.RAnoH,self.DecnoH,'gs',markersize=5)
	#legend(loc='center left',numpoints=1)
	legend(numpoints=1, prop={'size':10},loc='lower right')
	gca().invert_xaxis()
	s=self.name+' PRIMUS Coverage'
	title(s)
	savefig(self.name+'primuscoverage.eps')

    def plotnmbscoverage(self):
	figure()
	clf()
	plot(self.ramp[self.greenflagmp],self.decmp[self.greenflagmp],'co',markersize=18,label='Green')
	plot(self.ramp,self.decmp,'bo',markersize=10, label='Halpha+specz')
	xmin,xmax=xlim()
	ymin,ymax=ylim()


	plot(self.ramp,self.decmp,'go',markersize=8,label='Halpha+photz')
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=5,label='Specz')
	plot(self.ram,self.decm,'ro',markersize=5,label='Specz')
	#plot(self.zphotRAm,self.zphotDecm,'ko',markersize=2,label='Photz')


	if self.name.find('SXDS') > -1:
		#plot(prim.p.RA[prim.p.XMMFLAG],prim.p.DEC[prim.p.XMMFLAG],'ko',markersize=2,label='PRIMUS')
		print 'no nmbs coverage'
	if self.name.find('COSMOS') > -1:
		plot(nmbs.ra[nmbs.fieldflag],nmbs.dec[nmbs.fieldflag],'ko',markersize=2,label='NMBS')



	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)


	#plot(self.RAnoH,self.DecnoH,'gs',markersize=5)
	#legend(loc='center left',numpoints=1)
	legend(numpoints=1, prop={'size':10},loc='lower right')
	gca().invert_xaxis()
	s=self.name+' NMBS Coverage'
	title(s)
	savefig(self.name+'nmbscoverage.eps')

    def plotpositions(self):
	figure()
	clf()

	plot(self.ramp[self.greenflagmp],self.decmp[self.greenflagmp],'co',markersize=18,label='Green')
	plot(self.ramp,self.decmp,'bo',markersize=10, label='Halpha+specz')
	xmin,xmax=xlim()
	ymin,ymax=ylim()


	plot(self.ramp,self.decmp,'go',markersize=8,label='Halpha+photz')
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=5,label='Specz')
	plot(self.ram,self.decm,'ro',markersize=5,label='Specz')
	#plot(self.zphotRAm,self.zphotDecm,'ko',markersize=2,label='Photz')

	if self.name.find('SXDS') > -1:
		plot(prim.p.RA[prim.xmmflag],prim.p.DEC[prim.xmmflag],'ko',markersize=2,label='PRIMUS')
	if self.name.find('COSMOS') > -1:
		plot(prim.p.RA[prim.cosmosflag],prim.p.DEC[prim.cosmosflag],'ko',markersize=2,label='PRIMUS')




	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)


	#plot(self.RAnoH,self.DecnoH,'gs',markersize=5)
	#legend(loc='center left',numpoints=1)
	legend(numpoints=1,prop={'size':10})

	#plot(self.zphotRA,self.zphotDec,'k.')
	axhline(y=-5.49,color='c')
	axhline(y=-5.025,color='c')
	axvline(x=34.65,color='c')
	axvline(x=34.17,color='c')
	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)




	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)

#	figure()
#	plot(self.zphotRAm,self.zphotDecm,'k.')
#	axhline(y=-5.49,color='c')
#	axhline(y=-5.025,color='c')
#	axvline(x=34.65,color='c')
#	axvline(x=34.17,color='c')
#	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
#	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)
	dx=0.08
	dy=dx

	if self.name.find('SXDS') >-1:
		ramin=34.
		ramax=35.
		decmin=-6.
		decmax=-5.
		primflag=prim.xmmflag

	if self.name.find('COSMOS') >-1:
		ramin=150.
		ramax=150.7
		decmin=1.7
		decmax=2.4
		primflag=prim.cosmosflag

	ramin=self.ramin
	ramax=self.ramax
	decmin=self.decmin
	decmax=self.decmax

	x=arange(ramin,(ramax+dx),dx)
	y=arange(decmin,(decmax+dx),dy)

	backgroundRA=self.zphotRAm
	backgroundDEC=self.zphotDecm
	


	z=zeros([len(y),len(x)],'f')
	ylength=len(y)
	xlength=len(x)
	for i in range(len(backgroundRA)):
	    xbin=int(round((backgroundRA[i]-x[0])/dx))
	    ybin=int(round((backgroundDEC[i]-y[0])/dy))
	    if (xbin > -1) & (xbin < xlength):
		if (ybin > -1) & (ybin < ylength):
		    z[ybin,xbin] += 1
		  
	#print z
	ncon=10
	contour(x,y,z,ncon)
	contourf(x,y,z,ncon,cmap=cm.Greys)

	xlim(xmax,xmin)
	ylim(ymin,ymax)
	title(self.name)

	xlim(xmax,xmin)
	ylim(ymin,ymax)

	xlim(ramax,ramin)
	ylim(decmin,decmax)
	fname=self.name+'positions.eps'
	savefig(fname)

    def plotpositionsprimus(self):
	figure()
	clf()
	plot(self.ramp[self.greenflagmp],self.decmp[self.greenflagmp],'co',markersize=18,label='Green')
	plot(self.ramp,self.decmp,'bo',markersize=10, label='Halpha+specz')
	xmin,xmax=xlim()
	ymin,ymax=ylim()


	plot(self.ramp,self.decmp,'go',markersize=8,label='Halpha+photz')
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=5,label='Specz')
	plot(self.ram,self.decm,'ro',markersize=5,label='Specz')
	#plot(self.zphotRAm,self.zphotDecm,'ko',markersize=2,label='Photz')


	if self.name.find('SXDS') > -1:
		plot(prim.p.RA[prim.xmmflag],prim.p.DEC[prim.xmmflag],'ko',markersize=2,label='PRIMUS')
	if self.name.find('COSMOS') > -1:
		plot(prim.p.RA[prim.cosmosflag],prim.p.DEC[prim.cosmosflag],'ko',markersize=2,label='PRIMUS')



	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)


	#plot(self.RAnoH,self.DecnoH,'gs',markersize=5)
	#legend(loc='center left',numpoints=1)
	legend(numpoints=1,prop={'size':10})

	#plot(self.zphotRA,self.zphotDec,'k.')
	axhline(y=-5.49,color='c')
	axhline(y=-5.025,color='c')
	axvline(x=34.65,color='c')
	axvline(x=34.17,color='c')
	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)




	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)

#	figure()
#	plot(self.zphotRAm,self.zphotDecm,'k.')
#	axhline(y=-5.49,color='c')
#	axhline(y=-5.025,color='c')
#	axvline(x=34.65,color='c')
#	axvline(x=34.17,color='c')
#	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
#	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)
	dx=0.08
	dy=dx

	if self.name.find('SXDS') >-1:
		ramin=34.
		ramax=35.
		decmin=-6.
		decmax=-5.
		primflag=prim.xmmflag

	if self.name.find('COSMOS') >-1:
		ramin=150.
		ramax=150.7
		decmin=1.7
		decmax=2.4
		primflag=prim.cosmosflag

	ramin=self.ramin
	ramax=self.ramax
	decmin=self.decmin
	decmax=self.decmax

	x=arange(ramin,(ramax+dx),dx)
	y=arange(decmin,(decmax+dx),dy)

	backgroundRA=prim.p.RA[primflag]
	backgroundDEC=prim.p.DEC[primflag]
	


	z=zeros([len(y),len(x)],'f')
	ylength=len(y)
	xlength=len(x)
	for i in range(len(backgroundRA)):
	    xbin=int(round((backgroundRA[i]-x[0])/dx))
	    ybin=int(round((backgroundDEC[i]-y[0])/dy))
	    if (xbin > -1) & (xbin < xlength):
		if (ybin > -1) & (ybin < ylength):
		    z[ybin,xbin] += 1
		  
	#print z
	ncon=10
	contour(x,y,z,ncon)
	contourf(x,y,z,ncon,cmap=cm.Greys)

	xlim(xmax,xmin)
	ylim(ymin,ymax)
	title(self.name)

	xlim(xmax,xmin)
	ylim(ymin,ymax)

	xlim(ramax,ramin)
	ylim(decmin,decmax)
	fname=self.name+'positionsprimus.eps'
	savefig(fname)

    def plotpositionsprimushex(self):
	figure()
	clf()
	plot(self.ramp[self.greenflagmp],self.decmp[self.greenflagmp],'co',markersize=18)
	plot(self.ramp,self.decmp,'bo',markersize=10, label='Halpha+specz')
	xmin,xmax=xlim()
	ymin,ymax=ylim()


	plot(self.ramp,self.decmp,'go',markersize=8,label='Halpha+photz')
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=5,label='Specz')
	plot(self.ram,self.decm,'ro',markersize=5,label='Specz')
	#plot(self.zphotRAm,self.zphotDecm,'ko',markersize=2,label='Photz')





	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)

	#plot(self.RAnoH,self.DecnoH,'gs',markersize=5)
	legend(numpoints=1,prop={'size':10})

	#plot(self.zphotRA,self.zphotDec,'k.')
	axhline(y=-5.49,color='c')
	axhline(y=-5.025,color='c')
	axvline(x=34.65,color='c')
	axvline(x=34.17,color='c')
	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)




	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)

#	figure()
#	plot(self.zphotRAm,self.zphotDecm,'k.')
#	axhline(y=-5.49,color='c')
#	axhline(y=-5.025,color='c')
#	axvline(x=34.65,color='c')
#	axvline(x=34.17,color='c')
#	axvline(x=34.197,ymin=0,ymax=.4,color='r',ls='--',lw=2)
#	axhline(y=-5.302,xmin=0.9,xmax=1,color='r',ls='--',lw=2)
	dx=0.08
	dy=dx
	if self.name.find('SXDS') >-1:
		ramin=34.
		ramax=35.
		decmin=-6.
		decmax=-5.
		primflag=prim.xmmflag

	if self.name.find('COSMOS') >-1:
		ramin=150.
		ramax=150.7
		decmin=1.7
		decmax=2.4
		primflag=prim.cosmosflag

	ramin=self.ramin
	ramax=self.ramax
	decmin=self.decmin
	decmax=self.decmax
	plot(prim.p.RA[primflag],prim.p.DEC[primflag],'k.')
	ngridy=int(round((decmax-decmin)/.02))
	ngridx=int(round((ramax-ramin)/.02))
	hexbin(prim.p.RA[primflag],prim.p.DEC[primflag],cmap=cm.gray_r,gridsize=[ngridx,ngridy],alpha=0.5,vmin=0,vmax=10,edgecolors='None')
	axis('equal')
	#x=arange(ramin,(ramax+dx),dx)
	#y=arange(decmin,(decmax+dx),dy)


	#z=zeros([len(y),len(x)],'f')
	#ylength=len(y)
	#xlength=len(x)
	#for i in range(len(self.zphotRAm)):
	#    xbin=int(round((self.zphotRAm[i]-x[0])/dx))
	#    ybin=int(round((self.zphotDecm[i]-y[0])/dy))
	#    if (xbin > -1) & (xbin < xlength):
	#	if (ybin > -1) & (ybin < ylength):
	#	    z[ybin,xbin] += 1
		  
	#print z
	#ncon=10
	#contour(x,y,z,ncon)
	#contourf(x,y,z,ncon,cmap=cm.Greys)

	
	#xlim(xmax,xmin)
	#ylim(ymin,ymax)
	title(self.name)

	#xlim(xmax,xmin)
	#ylim(ymin,ymax)

	xlim(ramax,ramin)
	ylim(decmin,decmax)
	fname=self.name+'positionsPRIMUShex.eps'
	savefig(fname)

    def plotcolormag(self):
	    figure()
	    clf()
	    #plot photoz's
	    vi=self.zphotVm-self.zphotim
	    plot(self.zphotim,vi,'ko',markersize=2)

	    vi=self.Vmp-self.Imp
	    plot(self.Imp,vi,'bo',markersize=8)

	    vi=self.Vm-self.Im
	    plot(self.Im,vi,'r.')
	    

	    x=arange(20.5,25.1,1.)
	    dm=4.403
	    y=-.09*(x-20)+dm-1.6
	    plot(x,y,'k-',label="_nolegend_")
	    y2=y-.3
	    plot(x,y2,ls='--',color='0.5',label="_nolegend_")
	    y2=y+.3
	    plot(x,y2,ls='--',color='0.5',label="_nolegend_")

	    dm=3.3
	    y=-.09*(x-20)+dm-1.6
	    plot(x,y,'r-',label="_nolegend_")

	    dm=3.6
	    y=-.09*(x-20)+dm-1.6
	    plot(x,y,'r-',label="_nolegend_")

	    axis([19.5,26.,0.,3.5])
	    title(self.name)
	    xlabel('I')
	    ylabel('V-I')
	    fname=self.name+'colormag.eps'
	    savefig(fname)


    def plotscatter(self):
	figure()
	clf()
	#ra=array(self.ra)
	#dec=array(self.dec)
	#plot(ra[self.membflag>0],dec[self.membflag>0],'ko',markersize=3)

	volume=(20.*self.EWm/350.)**2
	scatter(self.ram,self.decm,s=volume,color='b',alpha=0.5)
	plot(self.ramp,self.decmp,'ko',markersize=1)
	xmin,xmax=xlim()
	ymin,ymax=ylim()
	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)
	plot(self.zspecRAm,self.zspecDecm,'ro',markersize=2.5)
	plot(self.RAnoH,self.DecnoH,'gs',markersize=5)

	#a=imread('/Users/rfinn/research/FieldSurvey/images/SXSDSS.png')
	#imshow(a)
	#scatter(self.xm,self.ym,s=volume,alpha=0.5)
	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	#axis([34.7,34.25,-5.5,-5])
	xticks(fontsize=14)
	yticks(fontsize=14)

	xlim(xmax,xmin)
	ylim(ymin,ymax)


	file2=str(self.name)+'.spec.reg'
	output111=open(file2,'w')
	output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	for i in range(len(self.ram)):
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") #color=red\n"%(self.ram[i],self.decm[i],sqrt(volume[i])/15.*25.)
	    output111.write(string1)
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") # color=blue \n"%(self.ram[i],self.decm[i],2.66)
	    output111.write(string1)
	output111.close()

    def plotscattersig(self):#scale pt size w/local density
	figure()
	clf()
	#ra=array(self.ra)
	#dec=array(self.dec)
	#plot(ra[self.membflag>0],dec[self.membflag>0],'ko',markersize=3)

	volume=(1./min(self.sig5)*self.sig5)#**2
	scatter(self.ram,self.decm,s=volume,color='b',alpha=0.5)
	plot(self.ramp,self.decmp,'ko',markersize=1)
	xmin,xmax=xlim()
	ymin,ymax=ylim()
	#plot(self.zphotRAm,self.zphotDecm,'ro',markersize=2)
	#plot(self.zspecRAm,self.zspecDecm,'ro',markersize=2.5)

	#a=imread('/Users/rfinn/research/FieldSurvey/images/SXSDSS.png')
	#imshow(a)
	#scatter(self.xm,self.ym,s=volume,alpha=0.5)
	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	#axis([34.7,34.25,-5.5,-5])
	xticks(fontsize=14)
	yticks(fontsize=14)

	xlim(xmax,xmin)
	ylim(ymin,ymax)


	file2=str(self.name)+'.spec.reg'
	output111=open(file2,'w')
	output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	for i in range(len(self.ram)):
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") #color=red\n"%(self.ram[i],self.decm[i],sqrt(volume[i])/15.*25.)
	    output111.write(string1)
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") # color=blue \n"%(self.ram[i],self.decm[i],2.66)
	    output111.write(string1)
	output111.close()


    def plotscatterp(self):#using supermembflag information
	figure()
	clf()
	#ra=array(self.ra)
	#dec=array(self.dec)
	#plot(ra[self.membflag>0],dec[self.membflag>0],'ko',markersize=3)

	volume=(20.*self.EWmp/350.)**2
	scatter(self.ramp,self.decmp,s=volume,alpha=0.5)
	#plot(self.ramp,self.decmp,'ko',markersize=3)
	#a=imread('/Users/rfinn/research/FieldSurvey/images/SXSDSS.png')
	#imshow(a)
	#scatter(self.xm,self.ym,s=volume,alpha=0.5)
	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)
	xmin,xmax=xlim()
	xlim(xmax,xmin)
	#ylim(-5.5,-5.)
	title(self.name)

	file2=str(self.name)+'.supermemb.reg'
	output111=open(file2,'w')
	output111.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	for i in range(len(self.ramp)):
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") #color=red\n"%(self.ramp[i],self.decmp[i],sqrt(volume[i])/15.*25.)
	    output111.write(string1)
	    string1 = "circle(%12.8f, %12.8f, %5.1f\") # color=blue \n"%(self.ramp[i],self.decmp[i],2.66)
	    output111.write(string1)
	output111.close()


    def plotscatterwim(self):
	figure()
	clf()
	volume=(20.*self.EWm/max(self.EWm))**2
	#scatter(self.ram,self.decm,s=volume,alpha=0.5)
	#scatter(self.xm,self.ym,s=volume,alpha=0.5)
	xlabel('RA',fontsize=24)
	ylabel('Dec',fontsize=24)
	xticks(fontsize=14)
	yticks(fontsize=14)

	a=imread('/Users/rfinn/research/FieldSurvey/images/SXSDSS.png')
	#a.shape = max(self.x),max(self.y)
	imshow(a,aspect='equal')
	xmin,xmax=xlim()
	ymin,ymax=ylim()
	scatter(self.xm/max(self.x)*(xmax),(self.ym)/max(self.y)*(ymax),s=volume,alpha=0.5)
	#scatter(self.xm,(self.ym),s=volume,alpha=0.5)
	#xlim(xmax,xmin)



    def plotenv(self):
	figure()
	clf()
	title(self.name)
	subplot(2,2,1)
    
	plot(self.sig5[self.Jm>22],self.EWm[self.Jm>22],'bo')
	plot(self.sig5[self.Jm<22],self.EWm[self.Jm<22],'ro')
	#errorbar(self.sig5,self.EWm,yerr=self.EWerrm,fmt=None)
	xlabel(r'$\Sigma_5 \ (gal/Mpc^2)$')
	ylabel(r'$EW$')

	(xbin,ybin,ybinerr)=biniterr(self.sig5,self.EWm,5)
	#plot(xbin,ybin,'ro',markersize=8)
	#plot(xbin,ybin,'r')
	#errorbar(xbin,ybin,yerr=ybinerr,fmt=None)
	ax=gca()
	#ax.set_xscale('log')
	axhline(y=50.,ls='--',color='k')
	subplot(2,2,2)
	plot(self.sig10,self.EWm,'bo')
	xlabel(r'$\Sigma_{10} \ (gal/Mpc^2)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.sig10,self.EWm,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)


	ax=gca()
	#ax.set_xscale('log')

	subplot(2,2,3)
	plot(self.nearest,self.EWm,'bo')
	xlabel(r'$d_{Nearest} \ (kpc)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.nearest,self.EWm,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)

	ax=gca()
	#ax.set_xscale('log')
	
	subplot(2,2,4)
	plot(self.nearest,self.EWm,'bo')
	xlabel(r'$d_{Nearest} \ (kpc)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.nearest,self.EWm,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)
	xlim(0.,50.)


    def plotenvp(self):
	figure()
	clf()
	title(self.name)
	subplot(2,2,1)
    
	plot(self.sig5p[self.Jmp>22],self.EWmp[self.Jmp>22],'bo')
	plot(self.sig5p[self.Jmp<22],self.EWmp[self.Jmp<22],'ro')
	#errorbar(self.sig5,self.EWm,yerr=self.EWerrm,fmt=None)
	xlabel(r'$\Sigma_7 \ (gal/Mpc^2)$')
	ylabel(r'$EW$')

	(xbin,ybin,ybinerr)=biniterr(self.sig5p,self.EWmp,5)
	#plot(xbin,ybin,'ro',markersize=8)
	#plot(xbin,ybin,'r')
	#errorbar(xbin,ybin,yerr=ybinerr,fmt=None)
	ax=gca()
	#ax.set_xscale('log')
	axhline(y=50.,ls='--',color='k')
	subplot(2,2,2)
	plot(self.sig10p,self.EWmp,'bo')
	xlabel(r'$\Sigma_{10} \ (gal/Mpc^2)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.sig10p,self.EWmp,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)


	ax=gca()
	#ax.set_xscale('log')

	subplot(2,2,3)
	plot(self.nearestp,self.EWmp,'bo')
	xlabel(r'$d_{Nearest} \ (kpc)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.nearestp,self.EWmp,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)

	ax=gca()
	#ax.set_xscale('log')
	
	subplot(2,2,4)
	plot(self.nearestp,self.EWmp,'bo')
	xlabel(r'$d_{Nearest} \ (kpc)$')
	ylabel(r'$EW$')
	(xbin,ybin,ybinerr)=biniterr(self.nearestp,self.EWmp,5)
	plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)
	xlim(0.,50.)


    def plotlocalden(self):
	figure()
	clf()
	title(self.name)
    
	#plot(self.sig5[self.Jm>22],self.EWm[self.Jm>22],'bo')
	#plot(self.sig5[self.Jm<22],self.EWm[self.Jm<22],'ro')
	plot(self.sig5allp,self.EWallp,'bo')
	#errorbar(self.sig5,self.EWm,yerr=self.EWerrm,fmt=None)
	xlabel(r'$\Sigma_5 \ (gal/Mpc^2)$')
	ylabel(r'$EW$')

	(xbin,ybin,ybinerr)=biniterr(self.sig5allp,self.EWallp,9)
	#plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'r')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)
	ax=gca()
	ax.set_xscale('log')
	axhline(y=50.,ls='--',color='k')

    def plotsffracden(self):
	figure()
	clf()
	title(self.name)
	sfflag=self.EWall>.1
	self.sfflag=array(sfflag,'i')
	(xbin,ybin,ybinerr)=biniterr(self.sig5,self.sfflag,7)
	#plot(xbin,ybin,'ro',label='specz')
	#errorbar(xbin,ybin,yerr=ybinerr,fmt=None)

	sfflag=self.EWallp>.1
	self.sfflagallp=array(sfflag,'i')
	(xbin,ybin,ybinerr)=biniterr(self.sig5allp,self.sfflagallp,7)
	plot(xbin,ybin,'ko',label='photz')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)

	(xbin,ybin,ybinerr)=biniterr(self.sig5allp,self.photozHamatchflag,7)
	plot(xbin,ybin,'go',label='photz')
	errorbar(xbin,ybin,yerr=ybinerr,fmt=None)

	ax=gca()
	ax.set_xscale('log')
	#ax.set_yscale('log')
	legend(numpoints=1)
	xlabel(r'$\Sigma_5 \ (gal/Mpc^2)$')
	ylabel(r'$SF \ Fraction$')

	fname=self.name+'sff.eps'
	savefig(fname)
	
    def plotewmass(self):
	figure()
	clf()
	title(self.name)
	cut=median(self.sig5)
	plot(self.Jm,self.EWm,'ko')
	#errorbar(self.Jm,self.EWm,self.EWerrm,fmt=None)
	#plot(self.Jm[self.sig5>cut],self.EWm[self.sig5>cut],'ro')
	#plot(self.Jm[self.sig5<cut],self.EWm[self.sig5<cut],'bo')
	#(xbin,ybin,ybinerr)=biniterr(self.Jm[self.sig5>cut],self.EWm[self.sig5>cut],5)
	#plot(xbin,ybin,'ro',markersize=8)
	#plot(xbin,ybin,'r')
	#(xbin,ybin,ybinerr)=biniterr(self.Jm[self.sig5<cut],self.EWm[self.sig5<cut],5)
	#plot(xbin,ybin,'ro',markersize=8)
	#plot(xbin,ybin,'b')

	plot(cl1216.Jm,cl1216.EWm,'k^',mfc='w')
	plot(cl1054.Jm,cl1054.EWm,'ks',mfc='w')
	plot(cl1040.Jm,cl1040.EWm,'kd',mfc='w')

	Jm=cl1216.Jm.tolist()+cl1054.Jm.tolist()+cl1040.Jm.tolist()
	EWm=cl1216.EWm.tolist()+cl1054.EWm.tolist()+cl1040.EWm.tolist()
	
	Jm=array(Jm)
	EWm=array(EWm)
	(xbin,ybin,ybinerr)=biniterr(Jm,EWm,7)
	#plot(xbin,ybin,'ro',markersize=8)
	plot(xbin,ybin,'k')
	
	axis([19.5,24.5,0,300])
	#errorbar(self.sig5,self.EWm,yerr=self.EWerrm,fmt=None)
	xlabel(r'J')
	ylabel(r'EW')

    def savefigs(self):
	self.plotlocalden()
	s=self.name+'localden.eps'
	savefig(s)
	self.plotpositions()
	s=self.name+'positions.eps'
	savefig(s)
	self.plotsffracden()
	s=self.name+'sffracden.eps'
	savefig(s)

#sx,sxdhr=z08readfits.sxdss()
#getmembers(sx)

name='CL1216'
infile='/Users/rfinn/research/FieldSurvey/catalogs/cl1216.dat'
cl1216=cluster(infile)

name='CL1054'
infile='/Users/rfinn/research/FieldSurvey/catalogs/cl1054.dat'
cl1054=cluster(infile)

name='CL1040'
infile='/Users/rfinn/research/FieldSurvey/catalogs/cl1040.dat'
cl1040=cluster(infile)

zfield=0.805
angdist=my.DA(zfield,h100)#kpc/arcsec
#infile='/Users/rfinn/research/FieldSurvey/catalogs/SXDSS_NB1190_emitter_both.zspec.fits'
infile='/Users/rfinn/research/FieldSurvey/catalogs/SXDSS_HaNB1190_emitter_both.updated.zspec.fits'
infile2='/Users/rfinn/research/FieldSurvey/catalogs/dr1_photoz_radec_BVRiz_SXDS.dat'
decmin=-5.45
decmax=-5.02
ramin=34.17
ramax=34.65


name='SXDS-S'

#print infile2
sxs=field(infile,infile2)

infile='/Users/rfinn/research/FieldSurvey/catalogs/SXDSN_HaNB1190_emitter_both.updated.zspec.fits'
infile2='/Users/rfinn/research/FieldSurvey/catalogs/dr1_photoz_radec_BVRiz_SXDS.dat'
decmin=-5.45
decmax=-5.02
ramin=34.17
ramax=34.65


name='SXDS-N'

#print infile2
sxn=field(infile,infile2)

infile='/Users/rfinn/research/FieldSurvey/catalogs/SXDSW_HaNB1190_emitter_both.updated.zspec.fits'
infile2='/Users/rfinn/research/FieldSurvey/catalogs/dr1_photoz_radec_BVRiz_SXDS.dat'
decmin=-5.45
decmax=-5.02
ramin=34.17
ramax=34.65


name='SXDS'

#print infile2
sxw=field(infile,infile2)
#sx.plotenv()

#infile='/Users/rfinn/research/FieldSurvey/catalogs/COSMOS_NB1190_emitter_both.zspec.fits'
infile='/Users/rfinn/research/FieldSurvey/catalogs/COSMOS_HaNB1190_emitter_both.updated.zspec.fits'
infile2='/Users/rfinn/research/FieldSurvey/catalogs/dr1_photoz_radec_BVRiz_COSMOS.dat'
ramin=150.05
ramax=150.47
decmin=1.7947
decmax=2.25

name='COSMOS'
cos=field(infile,infile2)
#cos.plotenv()

prim=primus()
nmbs=newfirmMedBandSurvey()
