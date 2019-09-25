#!/usr/bin/env python
import sys, os
import Numeric as N
#import scipy
from math import *
from Scientific.Functions.Romberg import romberg
import pylab
from pyraf import iraf
from matplotlib import rc
#rc('font',family='serif', style='normal', variant='normal',weight='bold', stretch='normal', size='large')


h=0.7
Om=0.3
OmL=0.7
lamb = 23.8#micro meters
lamb = lamb*1.e-4 #cm
dlamb = 4.7#micro meters
dlamb = dlamb*1.e-4#cm
cl = 3.e10#cm/s

#rc(interactive,'False')
def func(x):
    return 1.0/(N.sqrt(Om*(1+x)**3 + OmL))

def dL(z,h):
    c=3.e5
    H0=100.*h
    s=romberg(func,(0.,z))#,tol=1.e-6)
    DL=c/H0*(1+z)*s #(Mpc/h)
    return DL
def dLcm(z,h):
    c=3.e5
    H0=100.*h
    s=romberg(func,(0.,z))#,tol=1.e-6)
    DL=c/H0*(1+z)*s #(Mpc/h)
    DL=DL*1.e6*3.08568025e18#convert to cm
    return DL


def Lum(dL,lamb,dlamb,flux):
    Lm =1.e-23*flux*4*N.pi*(dL**2.)*(cl/lamb**2.)*dlamb
    return Lm 

def SFR(l):
    star =(l)/(10.**9.8)
    return star

def flux24toSFR(flux24):
    flucs = flux24*1.e-6/.13  #Jansky's
    lum = Lum(dL,lamb,dlamb,flucs)#finds luminosity from flx, wvlngth, distance
    lum = lum/3.9e33#lum in solar lums
    Starform = SFR(lum)#SFR in solar masses per year
    return Starform


def plotallEWhavsgimmorph():

	pylab.errorbar(cl1216.gimbulgefrac,cl1216.gimEWha,cl1216.gimEWhaerr,cl1216.gimbulgefracerr, 'bo')
      	pylab.errorbar(cl1040.gimbulgefrac,cl1040.gimEWha,cl1040.gimEWhaerr,cl1040.gimbulgefracerr, 'rs')
	pylab.errorbar(cl1037.gimbulgefrac,cl1037.gimEWha,cl1037.gimEWhaerr,cl1037.gimbulgefracerr, 'go')
	pylab.xlabel('Bulgefrac')
	pylab.ylabel('Equivalent Widths')
	pylab.title('Equivalent Widths vs gimBulgefrac for Cl1216, CL1037, CL1040')
	pylab.axhline(y=10)
	pylab.axhline(y=-10)
	pylab.axhline(y=0)
	pylab.show()

def plotallEWhavsvismorph():

    pylab.errorbar(cl1216.visbulgefrac,cl1216.visEWha,cl1216.visEWhaerr,cl1216.visbulgefracerr, 'bo')
	pylab.errorbar(cl1040.visbulgefrac,cl1040.visEWha,cl1040.visEWhaerr,cl1040.visbulgefracerr, 'rs')
	pylab.errorbar(cl1037.visbulgefrac,cl1037.visEWha,cl1037.visEWhaerr,cl1037.visbulgefracerr, 'go')
	pylab.xlabel('Bulgefrac')
	pylab.ylabel('Equivalent Widths')
	pylab.title('Equivalent Widths vs visBulgefrac for Cl1216, CL1037, CL1040')
	pylab.axis([-10,15,-10,50])
       	pylab.axhline(y=10)
	pylab.axhline(y=-10)
	pylab.axhline(y=0)
	pylab.show()

class models:
    def __init__(self):
	print "model SED templates"

    def readDaleSED(self):
	file='/Users/rfinn/spitzer-teachers/programs/dale/spectra.dat'
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
	    if (float(t[0]) < 3.):
		continue
            ngal=ngal+1
        input.close()
	
	self.lam = N.zeros(ngal,'d')
	self.lamerr = N.zeros(ngal,'d')
	self.flux  = N.zeros([64,ngal],'d')
	self.alpha=N.zeros(64,'d')
	self.ircolor=N.zeros(64,'d')
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
	    if (float(t[0]) < 3.):#start models at 3 microns
		continue
	    self.lam[i] = float(t[0]) 
	    for j in range(64):
		self.flux[j][i] = float(t[int(j+1)])
	    i=i+1
	input.close()

      	input=open('/Users/rfinn/spitzer-teachers/programs/dale/alpha.dat','r')
	i=0
	for line in input:
	    t=line.split()
	    self.alpha[i]=float(t[0])
	    self.ircolor[i]=float(t[1])
	    i=i+1
	input.close()

	
	#calculate LTIR
	self.Ltir=N.zeros(len(self.alpha),'d')
	self.sfr=N.zeros(len(self.alpha),'d')
	self.logfnu24=N.zeros(len(self.alpha),'d')
	input=open('/Users/rfinn/spitzer-teachers/programs/dale/model.0.dat','r')
	(c1,c2,c3)=(1.559,0.7686,1.347)
	(nu1,nu2,nu3)=(3.e8/(23.8e-6),3.e8/(70.e-6),3.e8/(160.e-6))
	#print nu1,nu2,nu3
	i=0
	for line in input:
	    t=line.split()
	    self.logfnu24[i]=float(t[11])
	    self.Ltir[i]=(c1*nu1*10.**float(t[11])+c2*nu2*10.**float(t[12])+c3*nu3*10.**float(t[13]))*1.e-23#luminosity in erg/s
	    self.Ltir[i]=self.Ltir[i]/3.9e33#luminosity in solar luminosities
	    self.sfr[i]=self.Ltir[i]/(6.31e9)#convert to SFR
	    i += 1

	#for i in range(len(self.sfr)):
	#    if (self.alpha[i] > .999999999) & (self.alpha[i] < 2.500001):
	#	s = "%7.5f %7.5f %9.7e "% (self.alpha[i],self.logfnu24[i],self.Ltir[i])
	#	print i,s
	       

    def plotSEDtemplate(self):
	lam=N.array(self.lam,'f')	
	flux1=N.array(self.fl18,'f')
	flux2=N.array(self.fl20,'f')
	flux3=N.array(self.fl25,'f')
	flux4=N.array(self.fl30,'f')
	flux5=N.array(self.fl35,'f')

	pylab.plot(lam,flux1, 'b-')
	pylab.plot(lam,flux2, 'g-')
	pylab.plot(lam,flux3, 'r-')
	pylab.plot(lam,flux4, 'm-')
	pylab.plot(lam,flux5, 'k-')
	ax=pylab.gca()
	ax.set_xscale("log")
	pylab.xlabel(r'$\rm{Wavelength \ (\mu m)}$',fontsize=20)
	pylab.ylabel(r'$\rm{\nu F_\nu}} \ \ (Jy \ Hz)$',fontsize=20)
	pylab.title('Dust-only Models')

	pylab.show()

class ediscs24:
    def __init__(self):#individual galaxy properties
        print "Reading 24 micron data!"

    def readmastertable(self,file):
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

	self.ediscsID = []
	self.ra = N.zeros(ngal,'f')
	self.dec = N.zeros(ngal,'f')
	self.EWha = N.zeros(ngal,'f')
	self.EWhaerr = N.zeros(ngal,'f')
	self.SFR = N.zeros(ngal,'f')
	self.SFRerr = N.zeros(ngal,'f')
	self.matchflaghaediscs = N.zeros(ngal,'f')
	self.hstgimtype = N.zeros(ngal,'f')
	self.matchflagmorphgimediscs = N.zeros(ngal,'f')
	self.hstvisnumtype = N.zeros(ngal,'f')
	self.matchflagmorphvisediscs = N.zeros(ngal,'f')
	self.f24 = N.zeros(ngal,'d')
	self.errf24 = N.zeros(ngal,'d')
	self.matchflagediscs24 = N.zeros(ngal,'f')
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
	self.specmembership = N.zeros(ngal,'f')
	self.defmembership = N.zeros(ngal,'f')
	self.specz = N.zeros(ngal,'f')
	self.spectype = N.zeros(ngal,'f')
	self.specEWOII = N.zeros(ngal,'f')
	self.matchflagspecediscs = N.zeros(ngal,'f')
	self.newmatchflag = N.zeros(ngal,'f')
	self.bestz = N.zeros(ngal,'f')
	self.membflag = N.zeros(ngal,'f')
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
	self.UBlowzclust = N.zeros(ngal,'f')
	self.UBbestzclust = N.zeros(ngal,'f')
	self.UBhighzclust = N.zeros(ngal,'f')
	self.BVlowzclust = N.zeros(ngal,'f')
	self.BVbestzclust = N.zeros(ngal,'f')
	self.BVhighzclust = N.zeros(ngal,'f')
	self.UVlowzclust = N.zeros(ngal,'f')
	self.UVbestzclust = N.zeros(ngal,'f')
	self.UVhighzclust = N.zeros(ngal,'f')
	self.LJlowzclust = N.zeros(ngal,'f')
	self.LJbestzclust = N.zeros(ngal,'f')
	self.LJhighzclust = N.zeros(ngal,'f')
	self.LKlowzclust = N.zeros(ngal,'f')
	self.LKbestzclust = N.zeros(ngal,'f')
	self.LKhighzclust = N.zeros(ngal,'f')
	self.fluxK= N.zeros(ngal,'f')
	self.xcorr = N.zeros(ngal,'f')
	self.ycorr = N.zeros(ngal,'f')
	self.starflag = N.zeros(ngal,'f')
	self.SFflag = N.zeros(ngal,'f')
	self.nmatchediscs24  = N.zeros(ngal,'f')
	self.specEWOIIflag = N.zeros(ngal,'f')
	self.matchflagediscsirac  = N.zeros(ngal,'f')
	self.iracf1  = N.zeros(ngal,'d')
	self.iracf2   = N.zeros(ngal,'d')
	self.iracf3  = N.zeros(ngal,'d')
	self.iracf4  = N.zeros(ngal,'d')
	self.erriracf1  = N.zeros(ngal,'d')
	self.erriracf2  = N.zeros(ngal,'d')
	self.erriracf3  = N.zeros(ngal,'d')
	self.erriracf4  = N.zeros(ngal,'d')
	self.iracsexflag0  = N.zeros(ngal,'f')
	self.iracsexflag1  = N.zeros(ngal,'f')
	self.iracwch1  = N.zeros(ngal,'f')
	self.iracwch2  = N.zeros(ngal,'f')
	self.iracwch3  = N.zeros(ngal,'f')
	self.iracwch4  = N.zeros(ngal,'f')
	self.iracwmin  = N.zeros(ngal,'f')
	self.nmatchediscsirac  = N.zeros(ngal,'f')
	self.L24 = N.zeros(ngal,'d')
	self.errL24  = N.zeros(ngal,'d')
	self.LHa  = N.zeros(ngal,'d')
	self.errLHa  = N.zeros(ngal,'d')

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
	    index=N.arange(1,(len(t)),1)
	    for j in index:
		t[j]=float(t[j])
	    self.ediscsID.append(t[0])
	    #print "number of fields = ",len(t)
	    (self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EWha[i],self.EWhaerr[i],self.SFR[i],self.SFRerr[i],self.matchflaghaediscs[i],self.SFflag[i],self.hstgimtype[i],self.matchflagmorphgimediscs[i],self.hstvisnumtype[i],self.matchflagmorphvisediscs[i],self.matchflagediscs24[i],self.f24[i],self.errf24[i],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i], self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i],self.membflag[i],self.specmembership[i],self.defmembership[i],self.specz[i],self.spectype[i],self.specEWOII[i],self.matchflagspecediscs[i], self.specEWOIIflag[i],self.bestz[i], self.lowz[i], self.highz[i], self.wmin[i], self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i],self.nmatchediscsirac[i],self.L24[i],self.errL24[i],self.LHa[i],self.errLHa[i])=t[1:(len(t))]


	    i=i+1
							
        input.close()

	lumdist=dLcm(z,h)
	L24=self.f24*1.e-6 #micro Jansky's to Jansky's
	L24err=self.errf24*1.e-6 #micro Jansky's to Jansky's
	L24=1.e-23*L24#erg/s/cm^2/Hz
	L24err=1.e-23*L24err#erg/s/cm^2/Hz
	L24=L24*4*N.pi*((lumdist)**2.)#erg/s/Hz
	L24err=L24err*4*N.pi*((lumdist)**2.)#erg/s/Hz
	#self.nufnu24=self.L24*3.e8/23.8e-6#erg/s * nu
	#self.nufnu24err=self.L24err*3.e8/23.8e-6#erg/s *nu
	self.nufnu24=L24*3.e8/23.8e-6#erg/s * nu
	self.nufnu24err=L24err*3.e8/23.8e-6#erg/s *nu
	self.Ltir=L24*13./(3.9e33)#convert to total IR luminosity using Dale et al 2005, Fig 16
	self.sfr24=self.Ltir/(6.31e9)#convert to solar luminosities and then sfr


	self.nufnuch1=self.iracf1*1.e-23*4*N.pi*lumdist**2*3.e8/3.6e-6
	self.nufnuch2=self.iracf2*1.e-23*4*N.pi*lumdist**2*3.e8/4.5e-6
	self.nufnuch3=self.iracf3*1.e-23*4*N.pi*lumdist**2*3.e8/5.8e-6
	self.nufnuch4=self.iracf4*1.e-23*4*N.pi*lumdist**2*3.e8/8.0e-6
	self.nufnuch1err=self.erriracf1*1.e-23*4*N.pi*lumdist**2*3.e8/3.6e-6
	self.nufnuch2err=self.erriracf2*1.e-23*4*N.pi*lumdist**2*3.e8/4.5e-6
	self.nufnuch3err=self.erriracf3*1.e-23*4*N.pi*lumdist**2*3.e8/5.8e-6
	self.nufnuch4err=self.erriracf4*1.e-23*4*N.pi*lumdist**2*3.e8/8.0e-6
	#read in Greg's filter data to convert rest-frame L to nu Fnu
	input=open('filterfile.dat','r')
        i=0
	filt=[]#filtername
	lamb=[]#central wavelength
	zpconv=[]#conersions from AB to Johnson
	flamconv=[]#convert from flam to fnu
	filtfile=[]
	msol=[]
	lsol=[]#solar lum in erg/s/A
        for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
	    t=line.split()
	    lamb.append(float(t[1]))
	    flamconv.append(float(t[3]))
	    lsol.append(float(t[6]))
	lamb=N.array(lamb,'f')
	flamconv=N.array(flamconv,'f')
	lsol=N.array(lsol,'f')
	self.fnuU=self.LUbestzclust*1.e10*flamconv[0]*lsol[0]/h**2
	self.fnuB=self.LBbestzclust*1.e10*flamconv[1]*lsol[1]/h**2
	self.fnuV=self.LVbestzclust*1.e10*flamconv[2]*lsol[2]/h**2
	self.fnuR=self.LRbestzclust*1.e10*flamconv[3]*lsol[3]/h**2
	self.fnuI=self.LIbestzclust*1.e10*flamconv[4]*lsol[4]/h**2
	self.fnuJ=self.LJbestzclust*1.e10*flamconv[8]*lsol[8]/h**2
	self.fnuK=self.LKbestzclust*1.e10*flamconv[10]*lsol[10]/h**2
	self.fnuUerrp=(self.LUhighzclust-self.LUbestzclust)*1.e10*flamconv[0]*lsol[0]/h**2
	self.fnuBerrp=(self.LBhighzclust-self.LBbestzclust)*1.e10*flamconv[1]*lsol[1]/h**2
	self.fnuVerrp=(self.LVhighzclust-self.LVbestzclust)*1.e10*flamconv[2]*lsol[2]/h**2
	self.fnuRerrp=(self.LRhighzclust-self.LRbestzclust)*1.e10*flamconv[3]*lsol[3]/h**2
	self.fnuIerrp=(self.LIhighzclust-self.LIbestzclust)*1.e10*flamconv[4]*lsol[4]/h**2
	self.fnuJerrp=(self.LJhighzclust-self.LJbestzclust)*1.e10*flamconv[8]*lsol[8]/h**2
	self.fnuKerrp=(self.LKhighzclust-self.LKbestzclust)*1.e10*flamconv[10]*lsol[10]/h**2
	self.fnuUerrm=(self.LUbestzclust-self.LUlowzclust)*1.e10*flamconv[0]*lsol[0]/h**2
	self.fnuBerrm=(self.LBbestzclust-self.LBlowzclust)*1.e10*flamconv[1]*lsol[1]/h**2
	self.fnuVerrm=(self.LVbestzclust-self.LVlowzclust)*1.e10*flamconv[2]*lsol[2]/h**2
	self.fnuRerrm=(self.LRbestzclust-self.LRlowzclust)*1.e10*flamconv[3]*lsol[3]/h**2
	self.fnuIerrm=(self.LIbestzclust-self.LIlowzclust)*1.e10*flamconv[4]*lsol[4]/h**2
	self.fnuJerrm=(self.LJbestzclust-self.LJlowzclust)*1.e10*flamconv[8]*lsol[8]/h**2
	self.fnuKerrm=(self.LKbestzclust-self.LKlowzclust)*1.e10*flamconv[10]*lsol[10]/h**2
	#print "self.fnuU[5] = ",self.fnuU[5],self.LUbestzclust[5]
	self.nufnuU=self.fnuU*3.e8/(lamb[0]*1.e-10)
	self.nufnuB=self.fnuB*3.e8/(lamb[1]*1.e-10)
	self.nufnuV=self.fnuV*3.e8/(lamb[2]*1.e-10)
	self.nufnuR=self.fnuR*3.e8/(lamb[3]*1.e-10)
	self.nufnuI=self.fnuI*3.e8/(lamb[4]*1.e-10)
	self.nufnuJ=self.fnuJ*3.e8/(lamb[8]*1.e-10)
	self.nufnuK=self.fnuK*3.e8/(lamb[10]*1.e-10)
	self.nufnuUerrp=self.fnuUerrp*3.e8/(lamb[0]*1.e-10)
	self.nufnuBerrp=self.fnuBerrp*3.e8/(lamb[1]*1.e-10)
	self.nufnuVerrp=self.fnuVerrp*3.e8/(lamb[2]*1.e-10)
	self.nufnuRerrp=self.fnuRerrp*3.e8/(lamb[3]*1.e-10)
	self.nufnuIerrp=self.fnuIerrp*3.e8/(lamb[4]*1.e-10)
	self.nufnuJerrp=self.fnuJerrp*3.e8/(lamb[8]*1.e-10)
	self.nufnuKerrp=self.fnuKerrp*3.e8/(lamb[10]*1.e-10)
	self.nufnuUerrm=self.fnuUerrm*3.e8/(lamb[0]*1.e-10)
	self.nufnuBerrm=self.fnuBerrm*3.e8/(lamb[1]*1.e-10)
	self.nufnuVerrm=self.fnuVerrm*3.e8/(lamb[2]*1.e-10)
	self.nufnuRerrm=self.fnuRerrm*3.e8/(lamb[3]*1.e-10)
	self.nufnuIerrm=self.fnuIerrm*3.e8/(lamb[4]*1.e-10)
	self.nufnuJerrm=self.fnuJerrm*3.e8/(lamb[8]*1.e-10)
	self.nufnuKerrm=self.fnuKerrm*3.e8/(lamb[10]*1.e-10)
	#print "self.nufnuU[5], fnu[5], LUbestzclust[5] = ",self.nufnuU[5],self.fnuU[5],self.LUbestzclust[5]

	for i in range(len(self.newmatchflag)):
	    if (self.matchflagspecediscs[i] > 0.):
		if (self.specmembership[i] > 0.):
		    self.newmatchflag[i] = 1.
	    if (self.matchflagspecediscs[i] > 0.):
		if (self.specmembership[i] < 1.):
		    self.newmatchflag[i] = 0.
	    if (self.matchflagspecediscs[i] < 1.):
		self.newmatchflag[i] = -1.	    
		
	for i in range (len(self.newmatchflag)):
	    if self.matchflagmorphvisediscs[i] < 1.:
		self.hstvisnumtype[i] = 100.
	

    def findDale24(self):#find index of Dale model that is closest to observed wavelength of 24um band
	l24obs=23.8/(1.+z)
	min=100000.
	self.i24=-99
	for i in range(len(mod.lam)):
	    d=abs(mod.lam[i]-l24obs)
	    if d < min:
		min=d
		self.i24=i


    def readGregSED(self):
	file='Cl1037mastertable24.nuLnu.out'
	input=open(file,'r')
        #get number of galaxies
        ngal=0
	nwave=0
	for line in input:
            if line.find('#') > -1: #skip lines with '#' in them
                continue
            if line.find('\\') > -1: #skip lines with '#' in them
                continue
            if line.find('|') > -1: #skip lines with '#' in them
                continue
	    nwave += 1
	t=line.split()	
        ngal=len(t)-1#first column is wavelength, one additional column per galaxy
        input.close()
	#print "ngal, nwave = ",ngal,nwave
	self.greglam = N.zeros(nwave,'f')
	self.gregflux  = N.zeros([ngal,nwave],'f')
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
	    self.greglam[i] = float(t[0])
	    for j in range(ngal):
		self.gregflux[j][i] = float(t[int(j+1)])
	    i=i+1
	input.close()
	self.greglam=self.greglam/10000.#convert from A to um
    def plothavs24(self,name):
	    outfile=str(name)+'match24ha.dat'
	    output=open(outfile,'w')
	    flux1=[]
	    errflux1=[]
	    flux2=[]
	    errflux2=[]	
	    for i in range(len(self.EW)):
		    if self.matchflagha[i] > 0:
			    if self.matchflag24[i] > 0:
				    flux1.append(self.SFR[i])
				    flux2.append(self.flux24[i])  
				    errflux1.append(self.SFRerr[i])
				    errflux2.append(self.flux24err[i])
				    #print "found match with Halpha source!\n"

	    flux1=N.array(flux1,'f')
	    flux2=N.array(flux2,'f')
		
	    pylab.xlabel('H Alpha Flux - SFR')
	    pylab.ylabel('24 micron Flux')
	    title='24 Micron Flux Versus Ha Flux for '+str(name)
	    pylab.title(title)
	    pylab.errorbar(flux1,flux2,errflux2,errflux1, 'bo')
	    pylab.show()
	    
	    self.ha24flux1=flux1
	    self.ha24flux2=flux2
	    self.ha24errflux1=errflux1
	    self.ha24errflux2=errflux2

    def plotEWhavsgimmorph(self,name):
	    outfile=str(name)+'EWhavsgim.dat'
	    output=open(outfile,'w')
            bulgefrac=[]
	    EWha=[]
	    EWhaerr=[]
	    for i in range(len(self.EW)):
		    if self.matchflaggimtype[i] > 0:
			    if self.matchflagha[i] > 0:
				    bulgefrac.append(self.gimtype[i])
				    EWha.append(self.EW[i])
				    EWhaerr.append(self.EWerr[i])
				   
	    bulgefrac=N.array(bulgefrac,'f')
	    bulgefracerr=N.zeros(len(bulgefrac),'f')
	    EWha=N.array(EWha,'f')
	    EWhaerr=N.array(EWhaerr,'f')
	    
      	    pylab.xlabel('Morphology')
	    pylab.ylabel('Equivalent Widths')
	    title='Equivalent Widths vs Morphology for GIM2D '+str(name)
	    pylab.title(title)
	    pylab.errorbar(bulgefrac,EWha,EWhaerr,bulgefracerr, 'bo')
	    pylab.show()
	    
	    self.gimbulgefrac = bulgefrac
	    self.gimbulgefracerr = bulgefracerr
	    self.gimEWha = EWha
	    self.gimEWhaerr = EWhaerr

    def plotEWhavsvismorph(self,name):
	    outfile=str(name)+'EWhavsvis.dat'
            bulgefrac=[]
	    EWha=[]
	    EWhaerr=[]
	    for i in range(len(self.EW)):
		    if self.matchflagvistype[i] > 0:
			    if self.matchflagha[i] > 0:
				    bulgefrac.append(self.vistype[i])
				    EWha.append(self.EW[i])
				    EWhaerr.append(self.EWerr[i])
				   
	    bulgefrac=N.array(bulgefrac,'f')
	    bulgefracerr=N.zeros(len(bulgefrac),'f')
	    EWha=N.array(EWha,'f')
	    EWhaerr=N.array(EWhaerr,'f')
	    
      	    pylab.xlabel('Morphology')
	    pylab.ylabel('Equivalent Widths')
	    title='Equivalent Widths vs Morphology for Visual '+str(name)
	    pylab.title(title)
	    pylab.errorbar(bulgefrac,EWha,EWhaerr,bulgefracerr, 'gs')
	    pylab.axis([-10.,15.,-10.,50.])
	    pylab.show()

	    self.visbulgefrac = bulgefrac
	    self.visbulgefracerr = bulgefracerr
	    self.visEWha = EWha
	    self.visEWhaerr = EWhaerr

    def plotcolormag(self):
	    mi=[]
	    mv=[]
	    mr=[]
	    pylab.plot(self.magI,(self.magV-self.magI),'go')
	    pylab.xlabel('I band magnitude')
	    pylab.ylabel('V - I band magnitude')
	    pylab.title('V - I band magnitude Versus I band magnitude')
	    pylab.axis([14,30,-2,10])
	    pylab.show()


    def plotcolormagmemb(self):
	    mi=[]
	    mv=[]
	    mr=[]
	    for i in range(len(self.ediscsID)):
		    if self.membflag[i] > 0:
			    mi.append(self.magI[i])
			    mv.append(self.magV[i])
			    mr.append(self.magR[i])
	    mi=N.array(mi,'f')
	    mv=N.array(mv,'f')
	    mr=N.array(mr,'f')

	    pylab.plot(mi,(mv-mi),'go')
	    pylab.xlabel('I band magnitude')
	    pylab.ylabel('V - I band magnitude')
	    pylab.title('V - I band magnitude Versus I band magnitude')
	    pylab.axis([14,30,-2,10])
	    pylab.show()
    def plotallseds(self):
	for i in range(len(self.membflag)):
	    if self.defmembership[i] > 0:
		self.plot1sedmod(i)
		pylab.show()
		k = raw_input("press any key to continue, q to quit: ")
		if k.lower().startswith('q'): 
		    sys.exit()
		pylab.show(False)


    def plot1sed(self,i):#plot photometry w/out dust models
	pylab.clf()
	lamb=N.array([.36,.44,.55,.64,.79,1.26,23.8/(1+z)],'f')#UBVRIJK central wave
	dlamb=N.array([.15,.22,.16,.23,.19,.16,.23,4.7],'f')#UBVRIJK delta lambda
	j=0
	k=0

	if ((self.membflag[i] +self.specmembership[i]) < 1):
	    print "Warning: not a member"
	flux=[self.nufnuU[i],self.nufnuB[i],self.nufnuV[i],self.nufnuR[i],self.nufnuI[i],self.nufnuJ[i],self.nufnu24[i]]
	fluxerrp=[self.nufnuUerrp[i],self.nufnuBerrp[i],self.nufnuVerrp[i],self.nufnuRerrp[i],self.nufnuIerrp[i],self.nufnuJerrp[i],self.nufnu24err[i]]
	fluxerrm=[self.nufnuUerrm[i],self.nufnuBerrm[i],self.nufnuVerrm[i],self.nufnuRerrm[i],self.nufnuIerrm[i],self.nufnuJerrm[i],self.nufnu24err[i]]
	pylab.ylabel(r"$\rm{\nu L_\nu \ (erg \ s^{-1})}$", fontsize=16)
	pylab.xlabel(r"$\rm{Rest \ Wavelength \ (\mu m)}$", fontsize=16)
	
	#lam=mod.lam
	#modindex=N.array([15,20,25,30,39],'i')
	#model=N.zeros([len(modindex),len(mod.flux[0][:])],'d')
	#for j in range(len(modindex)):
	#    model[j][:]=mod.flux[modindex[j]][:]
	#self.i24=int(self.i24)
	#scale=N.zeros(len(modindex),'d')
	#sfr=N.zeros(len(modindex),'d')
	#Ltir=N.zeros(len(modindex),'d')

	#for j in range(len(scale)):
	#    scale[j]=N.log10(self.nufnu24[i])-model[j][self.i24]
	#    model[j][:]=10.**(model[j][:] + scale[j])
	#    sfr[j]=mod.sfr[modindex[j]]*10.**(scale[j]-3.)
	#    Ltir[j]=N.log10(mod.Ltir[modindex[j]]*10.**(scale[j]-3.)*3.9e33)
	#(line1,line2,line3,line4,line5)=pylab.plot(lam,model[0][:], 'y-',lam,model[1][:],'g-',lam,model[2][:], 'r-',lam,model[3][:], 'm-',lam,model[4][:], 'c-')
	#s=[]
	#for j in range(len(modindex)):
	    #print str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]])
	    #l=str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]]*10.**scale)
	#    l="%3.1f, %5.2f, %4.1f"%(mod.alpha[modindex[j]],Ltir[j],sfr[j])
	#    s.append(l)
	#pylab.legend((line1,line2,line3,line4,line5),(s),loc=(0.45,0.1))
	pylab.hold(True)
	#f=10**(self.gregflux[i][:])
	#pylab.plot(self.greglam,f,'b:')
	pylab.plot(lamb,flux,'ko',markersize=8)
	pylab.errorbar(lamb,flux,[fluxerrp,fluxerrm],fmt=None,ecolor='k')
	#mark positions of IRAC and 70 micron bands
	#spitzer=N.array([3.6,4.5,5.8,8.,70.],'f')/(1.+z)
	#for l in spitzer:
	#    s=str(l)
	#    pylab.axvline(x=l, color='k',linestyle='--',label=s)
	pylab.axis([.2,1000.,5.e41,1.e45])
	ax1=pylab.gca()
	ax1.set_yscale("log")
	ax1.set_xscale("log")
	s1="%2.1f" % (self.membflag[i])
	#sfr= "SFR(24) = %5.1f"%(self.sfr24[i])
	s=str(i)+", "+str(self.ediscsID[i])+", photmemb= "+str(s1)+", Ttype="+str(self.hstvisnumtype[i])
	pylab.title(s,fontsize=12)
	#pylab.text(0.495,0.04,sfr,transform=ax1.transAxes,fontsize=12)
	#print s
	#print sfr
	if self.matchflagspecediscs[i] > 0:
	    s1="%6.4f" % (self.specz[i])
	    s2="z="+str(s1)+", specmemb="+str(self.specmembership[i])
	    pylab.text(0.3,0.95,s2,horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontsize=14)
	pylab.savefig("plot1sed.eps")
    def plot1sedmod(self,i):
	if i > (len(self.specz)-1):
	    print "only have ",len(self.specz)," galaxies.  Prepare to crash!"

	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.4)
	lamb=N.array([.36,.44,.55,.64,.79,1.26,23.8/(1.+z),3.6/(1.+z),4.5/(1.+z),5.8/(1.+z),8.0/(1.+z)],'f')#UBVRIJK central wave
	dlamb=N.array([.15,.22,.16,.23,.19,.16,.23,4.7],'f')#UBVRIJK delta lambda
	j=0
	k=0

	if ((self.membflag[i] +self.specmembership[i]) < 1):
	    print "Warning: not a member"
	flux=[self.nufnuU[i],self.nufnuB[i],self.nufnuV[i],self.nufnuR[i],self.nufnuI[i],self.nufnuJ[i],self.nufnu24[i],self.nufnuch1[i],self.nufnuch2[i],self.nufnuch3[i],self.nufnuch4[i]]
	fluxerrp=[self.nufnuUerrp[i],self.nufnuBerrp[i],self.nufnuVerrp[i],self.nufnuRerrp[i],self.nufnuIerrp[i],self.nufnuJerrp[i],self.nufnu24err[i],self.nufnuch1err[i],self.nufnuch2err[i],self.nufnuch3err[i],self.nufnuch4err[i]]
	fluxerrm=[self.nufnuUerrm[i],self.nufnuBerrm[i],self.nufnuVerrm[i],self.nufnuRerrm[i],self.nufnuIerrm[i],self.nufnuJerrm[i],self.nufnu24err[i],self.nufnuch1err[i],self.nufnuch2err[i],self.nufnuch3err[i],self.nufnuch4err[i]]
	pylab.ylabel(r"$\rm{\nu L_\nu \ (erg \ s^{-1})}$", fontsize=16)
	pylab.xlabel(r"$\rm{Rest \ Wavelength \ (\mu m)}$", fontsize=16)
	
	lam=mod.lam
	modindex=N.array([15,20,25,30,39],'i')
	model=N.zeros([len(modindex),len(mod.flux[0][:])],'d')
	for j in range(len(modindex)):
	    model[j][:]=mod.flux[modindex[j]][:]
	self.i24=int(self.i24)
	scale=N.zeros(len(modindex),'d')
	sfr=N.zeros(len(modindex),'d')
	Ltir=N.zeros(len(modindex),'d')


	for j in range(len(scale)):
	    scale[j]=N.log10(self.nufnu24[i])-model[j][self.i24]
	    model[j][:]=10.**(model[j][:] + scale[j])
	    sfr[j]=mod.sfr[modindex[j]]*10.**(scale[j]-3.)
	    Ltir[j]=N.log10(mod.Ltir[modindex[j]]*10.**(scale[j]-3.)*3.9e33)
	(line1,line2,line3,line4,line5)=pylab.plot(lam,model[0][:], 'y-',lam,model[1][:],'g-',lam,model[2][:], 'r-',lam,model[3][:], 'm-',lam,model[4][:], 'c-')
	s=[]
	for j in range(len(modindex)):
	    #print str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]])
	    #l=str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]]*10.**scale)
	    l="%3.1f, %5.2f, %4.1f"%(mod.alpha[modindex[j]],Ltir[j],sfr[j])
	    s.append(l)
	pylab.legend((line1,line2,line3,line4,line5),(s),loc=(0.45,0.1))
	pylab.hold(True)
	#f=10**(self.gregflux[i][:])
	#pylab.plot(self.greglam,f,'b:')
	pylab.plot(lamb,flux,'ko',markersize=8)
	pylab.errorbar(lamb,flux,[fluxerrp,fluxerrm],fmt=None,ecolor='k')
	#mark positions of IRAC and 70 micron bands
	spitzer=N.array([3.6,4.5,5.8,8.,70.],'f')/(1.+z)
	for l in spitzer:
	    s=str(l)
	    pylab.axvline(x=l, color='k',linestyle='--',label=s)
	pylab.axis([.2,1000.,5.e41,1.e45])
	ax1=pylab.gca()
	ax1.set_yscale("log")
	ax1.set_xscale("log")
	s1="%2.1f" % (self.membflag[i])
	sfr= "SFR(24) = %5.1f"%(self.sfr24[i])
	s=str(i)+", "+str(self.ediscsID[i])+", photmemb= "+str(s1)+", Ttype="+str(self.hstvisnumtype[i])
	pylab.title(s,fontsize=12)
	pylab.text(0.495,0.04,sfr,transform=ax1.transAxes,fontsize=12)
	print s
	print sfr
	if self.matchflagspecediscs[i] > 0:
	    s1="%6.4f" % (self.specz[i])
	    s2="z="+str(s1)+", specmemb="+str(self.specmembership[i])
	    pylab.text(0.3,0.95,s2,horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontsize=14)
	pylab.savefig("plot1sedmod.eps")
  

    def plot1sedmodprop(self):
	print 'z=',z
	#rc('legend',numpoints=1,fontsize=14,markerscale=.8)
	rc('ytick',labelsize=30)
	rc('xtick',labelsize=30)
	#rc('xtick',labelsize=28)

	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.4)
	#lamb=N.array([.36,.44,.55,.64,.79,1.26,23.8/(1.+z),3.6/(1.+z),4.5/(1.+z),5.8/(1.+z),8.0/(1.+z)],'f')#UBVRIJK central wave
	lamb=N.array([23.8/(1.+z),70./(1+z)],'f')#UBVRIJK central wave

	lumdist=dLcm(z,h)
	f24=.263#MJy/sr
	area=N.pi**3*(10./3600./180.)**2
	print "area in sr = ",area
	f24=f24*1.e6*area#convert to Jy

	L24=1.e-23*f24#erg/s/cm^2/Hz
	L24=L24*4*N.pi*((lumdist)**2.)#erg/s/Hz
	#self.nufnu24=self.L24*3.e8/23.8e-6#erg/s * nu
	#self.nufnu24err=self.L24err*3.e8/23.8e-6#erg/s *nu
	nufnu24=L24*3.e8/23.8e-6#erg/s * nu

	f70=.902#MJy/sr
	f70=f70*1.e6*area#convert to Jy

	L70=1.e-23*f70#erg/s/cm^2/Hz
	L70=L70*4*N.pi*((lumdist)**2.)#erg/s/Hz
	nufnu70=L70*3.e8/70.e-6#erg/s * nu



	flux=N.array([nufnu24,nufnu70],'d')#UBVRIJK central wave
	print flux
	j=0
	k=0


	pylab.ylabel(r'$\nu$ L$_\nu$  (erg s$^{-1}$)', fontsize=36)
	pylab.xlabel(r'Rest Wavelength ($\mu$m)', fontsize=36)
	#pylab.ylabel('$L_x$',fontsize=40)
	lam=mod.lam
	modindex=N.array([15,20,25,30,39],'i')
	model=N.zeros([len(modindex),len(mod.flux[0][:])],'d')
	for j in range(len(modindex)):
	    model[j][:]=mod.flux[modindex[j]][:]
	self.i24=int(self.i24)
	scale=N.zeros(len(modindex),'d')
	sfr=N.zeros(len(modindex),'d')
	Ltir=N.zeros(len(modindex),'d')




	
	nufnu24=flux[0]

	for j in range(len(scale)):
	    scale[j]=N.log10(nufnu24)-model[j][self.i24]
	    model[j][:]=10.**(model[j][:] + scale[j])
	    sfr[j]=mod.sfr[modindex[j]]*10.**(scale[j]-3.)
	    Ltir[j]=N.log10(mod.Ltir[modindex[j]]*10.**(scale[j]-3.)*3.9e33)
	(line1,line2,line3,line4,line5)=pylab.plot(lam,model[0][:], 'k:',lam,model[1][:],'g--',lam,model[2][:], 'r-.',lam,model[3][:], 'm-',lam,model[4][:], 'c-',lw=3)
	s=[]
	for j in range(len(modindex)):
	    #print str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]])
	    #l=str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]]*10.**scale)
	    #l="%3.1f, %5.2f, %4.1f"%(mod.alpha[modindex[j]],Ltir[j],sfr[j])
	    l="log(L$_{IR})$=%5.2f, SFR=%4.2f"%(Ltir[j],sfr[j])
	    s.append(l)
	pylab.legend((line1,line2,line3,line4,line5),(s),loc=(0.15,0.03))
	pylab.hold(True)
	#f=10**(self.gregflux[i][:])
	#pylab.plot(self.greglam,f,'b:')
	pylab.plot(lamb,flux,'ko',markersize=18)
	#pylab.errorbar(lamb,flux,[fluxerrp,fluxerrm],fmt=None,ecolor='k')
	#mark positions of IRAC and 70 micron bands
	spitzer=N.array([23.8,70.],'f')/(1.+z)
	for l in spitzer:
	    s=str(l)
	    #pylab.axvline(x=l, color='k',linestyle='--',label=s)
	pylab.text(27.,4.e41,'24$\mu$m',horizontalalignment='center')
	pylab.text(70.,4.e41,'70$\mu$m',horizontalalignment='center')
	pylab.axis([3.,1000.,2.e39,1.e43])
	pylab.xticks(fontsize=32)
	pylab.yticks(fontsize=32)
	ax1=pylab.gca()
	ax1.set_yscale("log")
	ax1.set_xscale("log")
	#s1="%2.1f" % (self.membflag[i])
	#sfr= "SFR(24) = %5.1f"%(self.sfr24[i])
	#s=str(i)+", "+str(self.ediscsID[i])+", photmemb= "+str(s1)+", Ttype="+str(self.hstvisnumtype[i])
	#pylab.title(s,fontsize=12)
	#pylab.text(0.495,0.04,sfr,transform=ax1.transAxes,fontsize=12)
	#print s
	#print sfr
	#if self.matchflagspecediscs[i] > 0:
	#    s1="%6.4f" % (self.specz[i])
	#    s2="z="+str(s1)+", specmemb="+str(self.specmembership[i])
	#    pylab.text(0.3,0.95,s2,horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontsize=14)
	pylab.savefig("plot1sedmod.eps")

    def plot1sedmodDDT(self):
	print 'z=',z
	#rc('legend',numpoints=1,fontsize=14,markerscale=.8)
	rc('ytick',labelsize=30)
	rc('xtick',labelsize=30)
	#rc('xtick',labelsize=28)


	
	pylab.clf()
	pylab.subplots_adjust(left=0.1, right=.9,bottom=.1,top=0.9,wspace=.3,hspace=.4)
	#lamb=N.array([.36,.44,.55,.64,.79,1.26,23.8/(1.+z),3.6/(1.+z),4.5/(1.+z),5.8/(1.+z),8.0/(1.+z)],'f')#UBVRIJK central wave
	lamb=N.array([23.8/(1.+z),70./(1+z)],'f')#UBVRIJK central wave

	lumdist=dLcm(z,h)

	l24=3.*0.1*10.**(9.8)*3.9e33#expected luminosity at 24um from 1 Msun/yr source at z=0.42

	#nufnu24=l24*3.e8/23.8e-6#erg/s * nu
	#l24=l24*10.#10 Msun/yr source
	#l24=l24*.75
	f24=l24/(4.*N.pi*lumdist**2)#erg/s/cm^2
	#f24=f24/(3.e8*(4.7e-6)/(23.8e-6)**2)#erg/s/cm^2/Hz
	f24=f24/(3.e8/(23.8e-6))#erg/s/cm^2/Hz
	f24=f24/(1.e-23)#Jy

	print "24um flux in Jy = ",f24, f24*1.e6


	#f24=.263#MJy/sr
	area=N.pi**3*(3./3600./180.)**2#area of point source on MIPS
	#print "area in sr = ",area
	#f24=f24*1.e6*area#convert to Jy

	L24=1.e-23*f24#erg/s/cm^2/Hz
	L24=L24*4*N.pi*((lumdist)**2.)#erg/s/Hz
	nufnu24=L24*3.e8/23.8e-6#erg/s * nu
	#self.nufnu24err=self.L24err*3.e8/23.8e-6#erg/s *nu


	f70=1210.#uJy 10 cycles, 10 sec exp
	f70=f70*1.e-6
	L70=1.e-23*f70#erg/s/cm^2/Hz
	L70=L70*4*N.pi*((lumdist)**2.)#erg/s/Hz
	nufnu70=L70*3.e8/70.e-6#erg/s * nu

	firac=N.array([.716,1.25,8.05,9.39],'d')#uJy 5 repeat, 30 s
	firac=firac*1.e-6*1.e-23*4*N.pi*((lumdist)**2.)#erg/s/Hz
	wavirac=N.array([3.6,4.5,5.8,8.0],'d')
	nufnuirac=firac*3.e8/(wavirac*1.e-6)#erg/s * nu



	flux=N.array([nufnu24,nufnu70],'d')#UBVRIJK central wave
	print flux
	j=0
	k=0


	pylab.ylabel(r'$\nu$ L$_\nu$  (erg s$^{-1}$)', fontsize=36)
	pylab.xlabel(r'Rest Wavelength ($\mu$m)', fontsize=36)
	#pylab.ylabel('$L_x$',fontsize=40)
	lam=mod.lam
	modindex=N.array([15,20,25,30,39],'i')
	model=N.zeros([len(modindex),len(mod.flux[0][:])],'d')
	for j in range(len(modindex)):
	    model[j][:]=mod.flux[modindex[j]][:]
	self.i24=int(self.i24)
	scale=N.zeros(len(modindex),'d')
	sfr=N.zeros(len(modindex),'d')
	Ltir=N.zeros(len(modindex),'d')




	
	nufnu24=flux[0]

	for j in range(len(scale)):
	    scale[j]=N.log10(nufnu24)-model[j][self.i24]
	    model[j][:]=10.**(model[j][:] + scale[j])
	    sfr[j]=mod.sfr[modindex[j]]*10.**(scale[j]-3.)
	    Ltir[j]=N.log10(mod.Ltir[modindex[j]]*10.**(scale[j]-3.)*3.9e33)
	(line1,line2,line3,line4,line5)=pylab.plot(lam,model[0][:], 'k:',lam,model[1][:],'g--',lam,model[2][:], 'r-.',lam,model[3][:], 'm-',lam,model[4][:], 'c-',lw=3)
	s=[]
	for j in range(len(modindex)):
	    #print str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]])
	    #l=str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]]*10.**scale)
	    #l="%3.1f, %5.2f, %4.1f"%(mod.alpha[modindex[j]],Ltir[j],sfr[j])
	    l="log(L$_{IR})$=%5.2f, SFR=%4.2f"%(Ltir[j],sfr[j])
	    s.append(l)
	pylab.legend((line1,line2,line3,line4,line5),(s),loc=(0.15,0.03))
	pylab.hold(True)
	#f=10**(self.gregflux[i][:])
	#pylab.plot(self.greglam,f,'b:')
	pylab.plot(lamb,flux,'ko',markersize=18)
	pylab.plot(wavirac/(1.+z),nufnuirac,'ko',markersize=15)
	#pylab.errorbar(lamb,flux,[fluxerrp,fluxerrm],fmt=None,ecolor='k')
	#mark positions of IRAC and 70 micron bands
	spitzer=N.array([23.8,70.],'f')/(1.+z)
	for l in spitzer:
	    s=str(l)
	    #pylab.axvline(x=l, color='k',linestyle='--',label=s)
	pylab.text(24./(1.+z),3.e42,'24$\mu$m',horizontalalignment='center')
	pylab.text(70./(1.+z),6.e42,'70$\mu$m',horizontalalignment='center')
	pylab.text(5./(1.+z),2.e43,'IRAC',horizontalalignment='center')
	pylab.axis([2.,1000.,5.e39,1.1e44])
	pylab.xticks(fontsize=32)
	pylab.yticks(fontsize=32)
	ax1=pylab.gca()
	ax1.set_yscale("log")
	ax1.set_xscale("log")
	#s1="%2.1f" % (self.membflag[i])
	#sfr= "SFR(24) = %5.1f"%(self.sfr24[i])
	#s=str(i)+", "+str(self.ediscsID[i])+", photmemb= "+str(s1)+", Ttype="+str(self.hstvisnumtype[i])
	#pylab.title(s,fontsize=12)
	#pylab.text(0.495,0.04,sfr,transform=ax1.transAxes,fontsize=12)
	#print s
	#print sfr
	#if self.matchflagspecediscs[i] > 0:
	#    s1="%6.4f" % (self.specz[i])
	#    s2="z="+str(s1)+", specmemb="+str(self.specmembership[i])
	#    pylab.text(0.3,0.95,s2,horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes,fontsize=14)
	pylab.savefig("plot1sedmodDDT.eps")
  
    def plotsnr(self,name):
	outfile=str(name)+'snr.dat'
	output=open(outfile,'w')
	flux1=[]
	errflux1=[]
	flux2=[]
	errflux2=[]
	fluxsnr=[]
	
	for i in range(len(self.f24)):
	    if self.matchflagediscs24[i] > 0:
		fluxsnr = (self.f24[i])/(self.errf24[i])
		flux1.append(fluxsnr)
		flux2.append(self.f24[i])  
		#print "found match with Halpha source!\n"
	    i=i+1
	flux1=N.array(flux1,'f')
	flux2=N.array(flux2,'f')
		    
	pylab.xlabel('24 micron flux')
	pylab.ylabel('24 micron Flux/24 micron Flux error')
	title='24 Micron Flux/24 micron Flux error Versus 24 micron Flux for '+str(name)
	pylab.title(title)
	pylab.plot(flux2,flux1, 'bo')
	pylab.axhline(6.)
	k=6.*19.
	pylab.axhline(3.)
	k=6.*19.
	pylab.axvline(k)
	k=6.*26.
	pylab.axvline(k)
	pylab.axis([0,500,0,20])
	pylab.show()

    def plothist24(self):
	pylab.hist(self.f24,40)
	pylab.axis([0,500,0,100])
	pylab.show()

    def writeteacherfile(self,prefix):
	
	totsfr=0.
	file=str(prefix)+".teachers.dat"
	output=open(file,'w')
	file2=str(prefix)+".teacherswmag.dat"
	output2=open(file2,'w')
	for i in range(len(self.ediscsID)):
	    if self.wmin[i] > .3:
		s="%s %12.8f %12.8f %6.1f %8.2f %6.2f %4.1f %4.1f %8.4f \n"%(self.ediscsID[i],self.ra[i],self.dec[i],self.hstvisnumtype[i],self.f24[i],self.errf24[i],self.membflag[i],self.newmatchflag[i],N.log10(self.nufnu24[i]))
		s2="%s %12.8f %12.8f %6.1f %8.2f %6.2f %4.1f %4.1f %8.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n"%(self.ediscsID[i],self.ra[i],self.dec[i],self.hstvisnumtype[i],self.f24[i],self.errf24[i],self.membflag[i],self.newmatchflag[i],N.log10(self.nufnu24[i]),self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i])
		output.write(s)
		output2.write(s2)
	output.close()
	output2.close()

	file2=str(prefix)+'.reg'
	file3=str(prefix)+'.members.reg'
	output1=open(file2,'w')
	output2=open(file3,'w')
	output1.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	output2.write("global color=blue font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")
	j=0
	for i in range(len(self.f24)):
	    if self.wmin[i] > .3:
		string1 = "circle(%12.8f, %12.8f, 3\") \n"%(self.ra[i],self.dec[i])
		string2 = "circle(%12.8f, %12.8f, 4\") \n"%(self.ra[i],self.dec[i])
		output1.write(string1)
		if (self.membflag[i] > 0.):
		    if (self.matchflagspecediscs[i] > 0.):
			if (self.specmembership[i] > 0.):	
			    output2.write(string2)
			    totsfr=totsfr+self.sfr24[i]
			    j=j+1
			    #print "1, ",self.newmatchflag[i],self.membflag[i],j
		    elif (self.newmatchflag[i] < 0.):
			output2.write(string2)
			totsfr=totsfr+self.sfr24[i]
			j=j+1
			#print "2, ",self.newmatchflag[i],self.membflag[i],j
	output1.close()
	output2.close()
	print "Total SFR = ",totsfr

mod = models()
mod.readDaleSED()

gocl1216=0.
gocl1037=0.
gocl1232=0
gocl1227=0
if gocl1216 > 0:
    z=0.794
    cl1216 = ediscs24()
    cl1216.readmastertable('cl1216mastertable24.dat')
    cl1216.findDale24()
    #cl1216.readGregSED()
    #cl1216.writeteacherfile('cl1037')
    #cl1216.plotsnr('cl1216')
    #cl1216.plothist24()
    #print "index and z of spec members:"
    #for i in range(len(cl1216.membflag)):
	#if abs(cl1216.specz[i]-z) < 0.05:
	#    sz="%6.4f " % (cl1216.specz[i]) 
	#    print i, sz
    #gal=int(sys.argv[1])

    cl1216.plotallseds()
#cl1037.plot1sed(gal)
#cl1037.plot1sedmod(gal)

if gocl1037 > 0:
    z=0.57
    cl1037 = ediscs24()
    cl1037.readmastertable('cl1037mastertable24.dat')
    cl1037.findDale24()
    #cl1037.readGregSED()
    cl1037.writeteacherfile('cl1037')
    cl1037.plotsnr('Cl1037')
    cl1037.plothist24()
    print "index and z of spec members:"
    for i in range(len(cl1037.membflag)):
	if abs(cl1037.specz[i]-z) < 0.05:
	    sz="%6.4f " % (cl1037.specz[i]) 
	    print i, sz
    cl1037.plotallseds()
    #gal=int(sys.argv[1])
#cl1037.plot1sed(gal)
#cl1037.plot1sedmod(gal)

#cl1037.plothavs24('cl1037')
#cl1037.plotEWhavsgimmorph('cl1037')
#cl1037.plotEWhavsvismorph('cl1037')
#cl1037.plotcolormag()
#cl1037.plotcolormagmemb()
#cl1037.plotseds()

if gocl1232 > 0.:
    z=.54
    cl1232 = ediscs24()
    cl1232.readmastertable('cl1232mastertable24.dat')
    cl1232.findDale24()
    #cl1232.readGregSED()
    cl1232.writeteacherfile('cl1232')
    cl1232.plotsnr('Cl1232')
    cl1232.plothist24()
    print "index and z of spec members:"
    for i in range(len(cl1232.membflag)):
	if abs(cl1232.specz[i]-z) < 0.05:
	    sz="%6.4f " % (cl1232.specz[i]) 
	    print i, sz
if gocl1227 > 0.:
    z=.63
    cl1227 = ediscs24()
    cl1227.readmastertable('cl1227mastertable24.dat')
    cl1227.findDale24()
    #cl1227.readGregSED()
    cl1227.writeteacherfile('cl1227')
    cl1227.plotsnr('Cl1227')
    cl1227.plothist24()
    print "index and z of spec members:"
    for i in range(len(cl1232.membflag)):
	if abs(cl1227.specz[i]-z) < 0.05:
	    sz="%6.4f " % (cl1227.specz[i]) 
	    print i, sz


#z=0.036
#prop=ediscs24()
#prop.findDale24()
#prop.plot1sedmodprop()

z=0.35
ddt=ediscs24()
ddt.findDale24()
ddt.plot1sedmodDDT()

#1 	self.ediscsID = []
#2 	self.ra = N.zeros(ngal,'f')
#3 	self.dec = N.zeros(ngal,'f')
#4 	self.EWha = N.zeros(ngal,'f')
#5 	self.EWhaerr = N.zeros(ngal,'f')
#6 	self.SFR = N.zeros(ngal,'f')
#7 	self.SFRerr = N.zeros(ngal,'f')
#8 	self.matchflaghaediscs = N.zeros(ngal,'f')
#9 	self.hstgimtype = N.zeros(ngal,'f')
#10	self.matchflagmorphgimediscs = N.zeros(ngal,'f')
#11 	self.hstvisnumtype = N.zeros(ngal,'f')
#12 	self.matchflagmorphvisediscs = N.zeros(ngal,'f')
#13 	self.f24 = N.zeros(ngal,'f')
#14 	self.errf24 = N.zeros(ngal,'f')
#15 	self.matchflagediscs24 = N.zeros(ngal,'f')
#16 	self.misoV = N.zeros(ngal,'f')
#17 	self.misoeVapsim = N.zeros(ngal,'f')
#18 	self.misoR = N.zeros(ngal,'f')
#19 	self.misoeRapsim = N.zeros(ngal,'f')
#20	self.misoI = N.zeros(ngal,'f')
#21 	self.misoeIapsim = N.zeros(ngal,'f')
#22 	self.misoJ = N.zeros(ngal,'f')
#23 	self.misoeJapsim = N.zeros(ngal,'f')
#24 	self.misoK = N.zeros(ngal,'f')
#25 	self.misoeKapsim = N.zeros(ngal,'f')
#26 	self.magV = N.zeros(ngal,'f')
#27 	self.mageVapsim = N.zeros(ngal,'f')
#28 	self.magR = N.zeros(ngal,'f')
#29 	self.mageRapsim = N.zeros(ngal,'f')
#30	self.magI = N.zeros(ngal,'f')
#31 	self.mageIapsim = N.zeros(ngal,'f')
#32 	self.magJ = N.zeros(ngal,'f')
#33 	self.mageJapsim = N.zeros(ngal,'f')
#34 	self.magK = N.zeros(ngal,'f')
#35 	self.mageKapsim = N.zeros(ngal,'f')
#36 	self.membflag = N.zeros(ngal,'f')
#37 	self.specmembership = N.zeros(ngal,'f')
#38 	self.specz = N.zeros(ngal,'f')
#39 	self.spectype = N.zeros(ngal,'f')
#40	self.specEWOII = N.zeros(ngal,'f')
#41 	self.matchflagspecediscs = N.zeros(ngal,'f')
#42 	self.bestz = N.zeros(ngal,'f')
#43 	self.membflag = N.zeros(ngal,'f')
#44 	self.lowz = N.zeros(ngal,'f')
#45 	self.highz = N.zeros(ngal,'f')
#46 	self.wmin = N.zeros(ngal,'f')
#47 	self.Pclust = N.zeros(ngal,'f')
#48 	self.LUlowzclust = N.zeros(ngal,'f')
#49 	self.LUbestzclust = N.zeros(ngal,'f')
#50	self.LUhighzclust = N.zeros(ngal,'f')
#51 	self.LBlowzclust = N.zeros(ngal,'f')
#52 	self.LBbestzclust = N.zeros(ngal,'f')
#53 	self.LBhighzclust = N.zeros(ngal,'f')
#54 	self.LVlowzclust = N.zeros(ngal,'f')
#55 	self.LVbestzclust = N.zeros(ngal,'f')
#56 	self.LVhighzclust = N.zeros(ngal,'f')
#57 	self.LRlowzclust = N.zeros(ngal,'f')
#58 	self.LRbestzclust = N.zeros(ngal,'f')
#59 	self.LRhighzclust = N.zeros(ngal,'f')
#60	self.LIlowzclust = N.zeros(ngal,'f')
#61 	self.LIbestzclust = N.zeros(ngal,'f')
#62 	self.LIhighzclust = N.zeros(ngal,'f')
#63 	self.UBlowzclust = N.zeros(ngal,'f')
#64 	self.UBbestzclust = N.zeros(ngal,'f')
#65 	self.UBhighzclust = N.zeros(ngal,'f')
#66 	self.BVlowzclust = N.zeros(ngal,'f')
#67 	self.BVbestzclust = N.zeros(ngal,'f')
#68 	self.BVhighzclust = N.zeros(ngal,'f')
#69 	self.UVlowzclust = N.zeros(ngal,'f')
#70	self.UVbestzclust = N.zeros(ngal,'f')
#71 	self.UVhighzclust = N.zeros(ngal,'f')
#72 	self.LJlowzclust = N.zeros(ngal,'f')
#73 	self.LJbestzclust = N.zeros(ngal,'f')
#74 	self.LJhighzclust = N.zeros(ngal,'f')
#75 	self.LKlowzclust = N.zeros(ngal,'f')
#76 	self.LKbestzclust = N.zeros(ngal,'f')
#77 	self.LKhighzclust = N.zeros(ngal,'f')
#78 	self.xcorr = N.zeros(ngal,'f')
#79 	self.ycorr = N.zeros(ngal,'f')
#80	self.starflag = N.zeros(ngal,'f')
