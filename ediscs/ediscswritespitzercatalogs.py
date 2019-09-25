#!/usr/bin/env python

'''
read in 'cluster'mastertable.v7.dat and write out a fits version.  Also include other quantities that were calculated in spitzer1.py readmaster() function.


    Written by Rose Finn, date unknown, updated on July 2, 2013

    PURPOSE: 

      
    CALLING SEQUENCE

      from the terminal command line, type:


    INPUT PARAMETERS
      none
      
    OUTPUT PARAMETERS
      Creates a table
      

    EXAMPLES

    PROCEDURE


    REQUIRED PYTHON MODULES

    ADDITIONAL REQUIRED MODULES

    NOTES
      Required Files


'''
from pylab import *
import pyfits, os
from ediscsCommon import *
import mystuff as my
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

f24IRz=[]
f24IRconv=[]
errf24IRconv=[]
infile2=open('/Users/rfinn/research/clusters/spitzer/Mastertables/fluxconv24toIR.dat','r')
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
f24IRz=array(f24IRz,'f')
f24IRconv=array(f24IRconv,'f')
errf24IRconv=array(errf24IRconv,'f')
infile2.close()

class cluster:
    def __init__(self,prefix):
        self.prefix=prefix
        self.fullname=fullname[self.prefix]
        self.zcl=redshift[self.prefix]
        self.f80=F80[self.prefix]
        self.errf80=ErrorF80[self.prefix]
	dzmin=1000.
	for i in range(len(f24IRz)):
		dz=abs(self.zcl-f24IRz[i])
		if dz < dzmin:
			dzmin=dz
			conv=f24IRconv[i]			
			errconv=errf24IRconv[i]
	self.f24IRconv=conv
	self.errf24IRconv=errconv

	#print self.prefix," at z= %6.4f: flux 24 to IR conv = %4.1f +/- %4.1f (%4.1f )"%(self.zcl,self.f24IRconv,self.errf24IRconv,self.errf24IRconv/self.f24IRconv*100)

        self.readoldmastertable()
        self.writefitstable()
    def readoldmastertable(self):
	mastertable=homedir+'research/clusters/spitzer/MasterTables/'+self.fullname+'mastertable.v7.dat'
        print mastertable
	nmemb=0.
	nspecmemb=0.
        input=open(mastertable,'r')
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

	    #if float(t[35]) > 24.:#keep only galaxies brighter than I=24
	#	    continue
	    #if float(t[42]) < 0.5:#keep members only
		#    continue

	
            ngal=ngal+1
        input.close()


	self.specmembflag = zeros(ngal,'f')
	self.matchflagspec = zeros(ngal,'f')
	self.fflat24 = zeros(ngal,'f')
	self.errfflat24 = zeros(ngal,'f')
	self.matchflagediscsflat24 = zeros(ngal,'f')
	

        self.ediscsID = []
        self.ediscsoldID = []
	self.ra = zeros(ngal,'f')
	self.dec = zeros(ngal,'f')
	self.xcorr = zeros(ngal,'f')
	self.ycorr = zeros(ngal,'f')
	self.starflag = zeros(ngal,'f')
	self.EW = zeros(ngal,'f')#Halpha EW
	self.EWerr = zeros(ngal,'f')
	self.SFR = zeros(ngal,'f')#Halpha SFR
	self.SFRerr = zeros(ngal,'f')
	self.matchflagha = zeros(ngal,'i')#match flag ha
	self.SFflag = zeros(ngal,'f')#if SF 
	self.gimtype = zeros(ngal,'f')
	self.matchflagmorphgimtype = zeros(ngal,'i')	
	self.vistype = zeros(ngal,'f')
	self.matchflagvistype = zeros(ngal,'i')
	self.matchflag24 = zeros(ngal,'i')
	self.flux24 = zeros(ngal,'f')
	self.flux24err = zeros(ngal,'f')
	self.nmatchediscs24 = zeros(ngal,'f')

	self.misoV = zeros(ngal,'f')
	self.misoeVapsim = zeros(ngal,'f')
	self.misoR = zeros(ngal,'f')
	self.misoeRapsim = zeros(ngal,'f')
	self.misoI = zeros(ngal,'f')
	self.misoeIapsim = zeros(ngal,'f')
	self.misoJ = zeros(ngal,'f')
	self.misoeJapsim = zeros(ngal,'f')	
	self.misoK = zeros(ngal,'f')
	self.misoeKapsim = zeros(ngal,'f')
	self.magV = zeros(ngal,'f')
	self.mageVapsim = zeros(ngal,'f')
	self.magR = zeros(ngal,'f')
	self.mageRapsim = zeros(ngal,'f')
	self.magI = zeros(ngal,'f')
	self.mageIapsim = zeros(ngal,'f')
	self.magJ = zeros(ngal,'f')
	self.mageJapsim = zeros(ngal,'f')
	self.magK = zeros(ngal,'f')
	self.mageKapsim = zeros(ngal,'f')
	self.membflag = zeros(ngal,'i')
	self.newspecmatchflag = zeros(ngal,'i')	
	self.defmembflag = zeros(ngal,'i')
	self.specz = zeros(ngal,'f')
	self.spectype = zeros(ngal,'f')
	self.specEWOII = zeros(ngal,'f')
	self.matchflagspecediscs = zeros(ngal,'i')
	self.specEWOIIflag = zeros(ngal,'f')
	self.bestz = zeros(ngal,'f')
	self.lowz = zeros(ngal,'f')
	self.highz = zeros(ngal,'f')
	self.wmin = zeros(ngal,'f')
	self.Pclust = zeros(ngal,'f')
	self.LUlowzclust = zeros(ngal,'f')
	self.LUbestzclust = zeros(ngal,'f')
	self.LUhighzclust = zeros(ngal,'f')	
	self.LBlowzclust = zeros(ngal,'f')
	self.LBbestzclust = zeros(ngal,'f')
	self.LBhighzclust = zeros(ngal,'f')	
	self.LVlowzclust = zeros(ngal,'f')
	self.LVbestzclust = zeros(ngal,'f')
	self.LVhighzclust = zeros(ngal,'f')	
	self.LRlowzclust = zeros(ngal,'f')
	self.LRbestzclust = zeros(ngal,'f')
	self.LRhighzclust = zeros(ngal,'f')	
	self.LIlowzclust = zeros(ngal,'f')
	self.LIbestzclust = zeros(ngal,'f')
	self.LIhighzclust = zeros(ngal,'f')	
	self.LJlowzclust = zeros(ngal,'f')
	self.LJbestzclust = zeros(ngal,'f')
	self.LJhighzclust = zeros(ngal,'f')	
	self.LKlowzclust = zeros(ngal,'f')
	self.LKbestzclust = zeros(ngal,'f')
	self.LKhighzclust = zeros(ngal,'f')
	self.fluxK = zeros(ngal,'f')
	self.UBlowzclust = zeros(ngal,'f')
	self.UBbestzclust = zeros(ngal,'f')
	self.UBhighzclust = zeros(ngal,'f')
	self.BVlowzclust = zeros(ngal,'f')
	self.BVbestzclust = zeros(ngal,'f')
	self.BVhighzclust = zeros(ngal,'f')
	self.UVlowzclust = zeros(ngal,'f')
	self.UVbestzclust = zeros(ngal,'f')
	self.UVhighzclust = zeros(ngal,'f')
	self.matchflagediscsirac = zeros(ngal,'i')
	self.iracf1 = zeros(ngal,'f')
	self.iracf2 = zeros(ngal,'f')
	self.iracf3 = zeros(ngal,'f')
	self.iracf4 = zeros(ngal,'f')
	self.erriracf1 = zeros(ngal,'f')
	self.erriracf2 = zeros(ngal,'f')
	self.erriracf3 = zeros(ngal,'f')
	self.erriracf4 = zeros(ngal,'f')
	self.iracsexflag0 = zeros(ngal,'i')
	self.iracsexflag1 = zeros(ngal,'i')
	self.iracwch1 = zeros(ngal,'f')
	self.iracwch2 = zeros(ngal,'f')
	self.iracwch3 = zeros(ngal,'f')
	self.iracwch4 = zeros(ngal,'f')
	self.iracwmin = zeros(ngal,'f')
	self.nmatchediscsirac = zeros(ngal,'f')
	self.L24 = zeros(ngal,'d')
	self.errL24 = zeros(ngal,'d')
	self.LHa = zeros(ngal,'d')
	self.errLHa = zeros(ngal,'d')
	self.snr24 = zeros(ngal,'d')
	self.imagex24 = zeros(ngal,'d')
	self.imagey24 = zeros(ngal,'d')
	self.fap1 = zeros(ngal,'d')#24 um aperture fluxes
	self.fap2 = zeros(ngal,'d')
	self.fap3 = zeros(ngal,'d')
	self.fap4 = zeros(ngal,'d')#2.6 pixel radius
	self.fap5 = zeros(ngal,'d')
	self.fap6 = zeros(ngal,'d')
	self.fap7 = zeros(ngal,'d')
	self.fap8 = zeros(ngal,'d')
	self.fap9 = zeros(ngal,'d')
	self.fap10 = zeros(ngal,'d')
	self.errfap1 = zeros(ngal,'d')
	self.errfap2 = zeros(ngal,'d')
	self.errfap3 = zeros(ngal,'d')
	self.errfap4 = zeros(ngal,'d')
	self.errfap5 = zeros(ngal,'d')
	self.errfap6 = zeros(ngal,'d')
	self.errfap7 = zeros(ngal,'d')
	self.errfap8 = zeros(ngal,'d')
	self.errfap9 = zeros(ngal,'d')
	self.errfap10 = zeros(ngal,'d')


	self.LVlow = zeros(ngal,'f')
	self.LVbest = zeros(ngal,'f')
	self.LVhigh = zeros(ngal,'f')	
	self.LRlow = zeros(ngal,'f')
	self.LRbest = zeros(ngal,'f')
	self.LRhigh = zeros(ngal,'f')	
	self.UBlow = zeros(ngal,'f')
	self.UBbest = zeros(ngal,'f')
	self.UBhigh = zeros(ngal,'f')
	self.BVlow = zeros(ngal,'f')
	self.BVbest = zeros(ngal,'f')
	self.BVhigh = zeros(ngal,'f')
	self.UVlow = zeros(ngal,'f')
	self.UVbest = zeros(ngal,'f')
	self.UVhigh = zeros(ngal,'f')

	self.ra24 = zeros(ngal,'f')
	self.dec24 = zeros(ngal,'f')
	self.flux80flag = zeros(ngal,'i')


        input=open(mastertable,'r')
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
	    #if float(t[35]) > 24.:#keep only galaxies brighter than I=24
		#    continue

	    #if float(t[42]) < 0.5:#keep members only
		#    continue


	    #print "number of fields in a line = ",len(t)
	    #print i,t[len(t)-1]
	    self.ediscsID.append(t[0])
            self.ediscsoldID.append(t[len(t)-1])
	    for j in range(1,(len(t)-1)):
		t[j]=float(t[j])



            #print i,len(t),len(self.ra)
	    (self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EW[i],self.EWerr[i],self.SFR[i],self.SFRerr[i],self.matchflagha[i],self.SFflag[i],self.gimtype[i],self.matchflagmorphgimtype[i],self.vistype[i],self.matchflagvistype[i],self.matchflag24[i],self.flux24[i],self.flux24err[i],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.specz[i],self.spectype[i],self.specEWOII[i],self.matchflagspecediscs[i],self.specEWOIIflag[i],self.bestz[i],self.lowz[i],self.highz[i],self.wmin[i],self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i],self.nmatchediscsirac[i],self.L24[i],self.errL24[i],self.LHa[i],self.errLHa[i],self.snr24[i],self.imagex24[i],self.imagey24[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i],self.LVlow[i],  self.LVbest[i],  self.LVhigh[i],  self.LRlow[i],  self.LRbest[i],  self.LRhigh[i], self.UBlow[i],  self.UBbest[i],  self.UBhigh[i],  self.BVlow[i],  self.BVbest[i],  self.BVhigh[i],  self.UVlow[i],  self.UVbest[i],  self.UVhigh[i],self.ra24[i],self.dec24[i],self.flux80flag[i])=t[1:(len(t)-1)]
	    #print i,self.defmembflag[i],self.matchflagvistype[i],self.vistype[i],self.matchflag24[i]
	    if (self.defmembflag[i] > 0.) & (self.matchflagvistype[i] > 0.):
		k=k+1
		#print k,self.vistype[i],self.matchflag24[i]
	    i=i+1


	input.close()

	self.photmembflag=zeros(len(self.membflag),'f')
	for i in range(len(self.membflag)):
	    if (self.membflag[i] > 0.) & (self.wmin[i] > 0.4):
		self.photmembflag[i]=1.
	self.supermembflag=zeros(len(self.membflag),'f')
	for i in range(len(self.membflag)):
		if self.newspecmatchflag[i] > 0.:
			self.supermembflag[i]=1.
		elif (self.defmembflag[i] > 0.) & (self.wmin[i] > 0.4):
			self.supermembflag[i]=1.
	    

	pre=self.prefix


	self.MR = 4.28 - 2.5*log10(self.LRbestzclust)-25.+5.*log10(h)
	self.MU = 5.66 - 2.5*log10(self.LUbestzclust)-25.+5.*log10(h)
	#self.errMRlo = self.MR - (4.28 - 2.5*log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*log10(self.Rlumhi) -25.) - self.MR 
	self.MV = 4.82 - 2.5*log10(self.LVbestzclust) -25.+5.*log10(h)
	#self.errMvlo = self.MV - (4.82 - 2.5*log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*log10(self.vlumhi) - 25. ) - self.MV
	self.MB = 5.48 - 2.5*log10(self.LBbestzclust) -25.+5.*log10(h)#check zeropoint here...
	#B-band mass-to-light ratio (Bell et al 2003)
	self.stellmass=10.**(1.737*(self.MB-self.MV)-0.942)*self.LBbestzclust*1.e10/h**2#luminosities are in 10^10 h^-2 Lsol
	#self.stellmass=10.**(1.305*(self.MB-self.MV)-0.628)*self.LVbestzclust*1.e10/h**2#M* from MV (Bell et al 2003)
	self.MRbestz = 4.28 - 2.5*log10(self.LRbest)-25.+5.*log10(h)
	#self.errMRlo = self.MR - (4.28 - 2.5*log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*log10(self.Rlumhi) -25.) - self.MR 
	self.MVbestz = 4.82 - 2.5*log10(self.LVbest) -25.+5.*log10(h)
	self.MUbestz = 5.66 - 2.5*log10(self.LVbest+self.UVbest)-25.+5.*log10(h)
	self.MBbestz = 5.48 - 2.5*log10(self.LVbest+self.BVbest)-25.+5.*log10(h)
	#self.errMvlo = self.MV - (4.82 - 2.5*log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*log10(self.vlumhi) - 25. ) - self.MV

	#self.redflag=(self.MB-self.MV)>(-.022*(self.MV+20)+.65)#true if galaxy is redder than B-V cut
	#	    y=-.09*(x-20)+dm-1.6
	#self.redflag=(-1.*0.09*(self.cI-20.) + (1.9*(self.z-0.42) +3.74)-1.6 - (self.cV-self.cI)) < 0.3
	self.redflag=(-1.*0.09*(self.magI-20.) + myRSOffset[self.prefix]-1.6 - (self.magV-self.magI)) < 0.3
	#if (dmag < .3): #redsequence

	#self.redflag= (self.MU-self.MV)>(1.15-0.31*self.zcl-0.08*(self.MV+20))#true for red galaxies

	self.L24=self.L24/Lsol
	self.errL24=self.errL24/Lsol
	#self.f24c=self.fap4*self.apcor4#corrected to infinite aperture
	#self.errf24c=self.errfap4*self.apcor4
	self.f24c=array(self.flux24,'d')
	self.errf24c=array(self.flux24err,'d')
	#print 'Running into trouble here ',ngal,file,self.f24c
	self.fmin=min(self.f24c)#[where(abs(self.f24c/self.errf24c)>2.5)])
	self.fmax=max(self.f24c)#[where(abs(self.f24c/self.errf24c)>2.5)])
	conversion=1.e-6*1.e-23*(4.*pi*my.dLcm(self.zcl,h)**2)*(3.e8/23.8e-6)#Converts from microJy to erg/s
	Lir=self.f24c*conversion*self.f24IRconv#factor of ten converts get scaling to go from from nu fnu(24) to L(IR)
	self.Lir=array(Lir,'d')/Lsol#convert to solar luminosities
	errLir=self.errf24c*conversion*self.f24IRconv
	self.errLir=array(errLir,'d')/Lsol
	#print 'conversion = ',conversion
	#print self.Lir
	self.Lir80=self.f80*conversion*self.f24IRconv/Lsol#luminosity corresponding to f at 80% comple
	#print self.prefix," log10(Lir80) = %6.3f" %(log10(self.Lir80))


#	#calculate Lir for bestz - best photoz
#
#	f24conv=splev(self.bestz,tck24)
#	self.f24conv=f24conv
#	conversion=1.e-6*1.e-23*(4.*pi*my.dLcm(self.bestz,h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s
#	print len(f24conv),len(conversion)

#	Lir=self.f24c*conversion*f24conv#factor of ten converts get scaling to go from from nu fnu(24) to L(IR)
#	self.Lirbestz=array(Lir,'d')/Lsol#convert to solar luminosities
#	errLir=self.errf24c*conversion*self.f24IRconv
#	self.errLirbestz=array(errLir,'d')/Lsol


	for i in range(len(self.ediscsID)):#set Lir of undetected sources to Lir80, 80% completeness limit
		if (abs(self.matchflag24[i]) < 0.1):
			self.Lir[i]=self.Lir80
			self.errLir[i]=self.Lir80

	sfrconv=bellconv*Lsol
	self.SFRir=self.Lir*sfrconv
	self.SFRirerr=self.errLir*sfrconv
	print "Got ",len(self.membflag)," galaxies"

	#self.SFR80=self.Lir80*sfrconv
	#print self.prefix," SFR80 = %5.2f" %(self.SFR80)

    def writefitstable(self):

	c1=pyfits.Column(name='EDISCS-ID',format='A21',unit='',array=self.ediscsID)
	c2=pyfits.Column(name='EDISCS-ID-OLD',format='A21',unit='deg',array=self.ediscsoldID)
	c3=pyfits.Column(name='RA',format='D',unit='deg',array=self.ra)
	c4=pyfits.Column(name='DEC',format='D',unit='deg',array=self.dec)
	c5=pyfits.Column(name='xcorr',format='E',unit='pix',array=self.xcorr)
	c6=pyfits.Column(name='ycorr',format='E',unit='pix',array=self.ycorr)
	c7=pyfits.Column(name='starflag',format='E',unit='',array=self.starflag)
	c8=pyfits.Column(name='EW',format='E',unit='A',array=self.EW)
	c9=pyfits.Column(name='EWerr',format='D',unit='A',array=self.EWerr)
	c10=pyfits.Column(name='SFR',format='D',unit='Msun/yr',array=self.SFR)
	c11=pyfits.Column(name='SFRerr',format='D',unit='Msun/yr',array=self.SFRerr)
	c12=pyfits.Column(name='matchflaghalpha',format='L',unit='',array=(self.matchflagha > 0.1))
	c12a=pyfits.Column(name='onHaimageflag',format='L',unit='',array=(self.matchflagha > -0.1))
	c13=pyfits.Column(name='sfflag',format='I',unit='',array=self.SFflag)

#	c14=pyfits.Column(name='gimtype',format='D',unit='Msun/yr',array=self.gimtype)
#	c15=pyfits.Column(name='matchflaggimtype',format='L',unit='Msun/yr',array=self.matchflagmorphgimtype)
#	c16=pyfits.Column(name='vistype',format='I',unit='',array=self.vistype)
#	c17=pyfits.Column(name='matchflagvistype',format='L',unit='',array=self.matchflagvistype)

	c18=pyfits.Column(name='matchflag24',format='L',unit='',array=(self.matchflag24>0.1))
	c18a=pyfits.Column(name='on24imageflag',format='L',unit='',array=(self.matchflag24> -0.1))
	c19=pyfits.Column(name='flux24',format='D',unit='micro-Jy',array=self.flux24)
	c20=pyfits.Column(name='flux24err',format='D',unit='micro-Jy',array=self.flux24err)
	c21=pyfits.Column(name='nmatchediscs24',format='I',unit='',array=self.nmatchediscs24)
        c24=pyfits.Column(name='snr24',format='D',unit='',array=self.snr24)     
        c25=pyfits.Column(name='imagex24',format='E',unit='pix',array=self.imagex24)  
        c26=pyfits.Column(name='imagey24',format='E',unit='pix',array=self.imagey24)  
        c27=pyfits.Column(name='ra24',format='D',unit='deg',array=self.ra24)      
        c28=pyfits.Column(name='dec24',format='D',unit='deg',array=self.dec24)     
        c29=pyfits.Column(name='flux80flag',format='I',unit='',array=self.flux80flag)
        c30=pyfits.Column(name='L24',format='D',unit='Lsun',array=self.L24)       
        c31=pyfits.Column(name='L24err',format='D',unit='Lsun',array=self.errL24)    
        c32=pyfits.Column(name='Lir',format='D',unit='Lsun',array=self.Lir)       
        c33=pyfits.Column(name='errLir',format='D',unit='Lsun',array=self.errLir)    
        c34=pyfits.Column(name='SFRir',format='D',unit='Msun/yr',array=self.SFRir)     
        c35=pyfits.Column(name='SFRirerr',format='D',unit='Msun/yr',array=self.SFRirerr)  
        c36=pyfits.Column(name='matchflagediscsirac',format='I',unit='',array=self.matchflagediscsirac)
        c37=pyfits.Column(name='iracf1',format='D',unit='',array=self.iracf1)     
        c38=pyfits.Column(name='iracf2',format='D',unit='',array=self.iracf2)     
        c39=pyfits.Column(name='iracf3',format='D',unit='',array=self.iracf3)     
        c40=pyfits.Column(name='iracf4',format='D',unit='',array=self.iracf4)     
        c41=pyfits.Column(name='erriracf1',format='D',unit='',array=self.erriracf1)  
        c42=pyfits.Column(name='erriracf2',format='D',unit='',array=self.erriracf2)  
        c43=pyfits.Column(name='erriracf3',format='D',unit='',array=self.erriracf3)  
        c44=pyfits.Column(name='erriracf4',format='D',unit='',array=self.erriracf4)  
        c45=pyfits.Column(name='iracsexflag0',format='I',unit='',array=self.iracsexflag0)
        c46=pyfits.Column(name='iracsexflag1',format='I',unit='',array=self.iracsexflag1)
        c47=pyfits.Column(name='iracwch1',format='E',unit='',array=self.iracwch1)   
        c48=pyfits.Column(name='iracwch2',format='E',unit='',array=self.iracwch2)   
        c49=pyfits.Column(name='iracwch3',format='E',unit='',array=self.iracwch3)   
        c50=pyfits.Column(name='iracwch4',format='E',unit='',array=self.iracwch4)   
        c51=pyfits.Column(name='iracwmin',format='E',unit='',array=self.iracwmin)   
        c52=pyfits.Column(name='nmatchediscsirac',format='I',unit='',array=self.nmatchediscsirac)
        c53=pyfits.Column(name='matchflagmorphgimtype',format='I',unit='',array=self.matchflagmorphgimtype)
        c54=pyfits.Column(name='gimtype',format='E',unit='',array=self.gimtype)
        c55=pyfits.Column(name='matchflagvistype',format='I',unit='',array=self.matchflagvistype)
        c56=pyfits.Column(name='vistype',format='I',unit='',array=self.vistype)
        c57=pyfits.Column(name='misoV',format='E',unit='mag',array=self.misoV)
        c58=pyfits.Column(name='misoeVapsim',format='E',unit='',array=self.misoeVapsim)
        c59=pyfits.Column(name='misoR',format='E',unit='mag',array=self.misoR)
        c60=pyfits.Column(name='misoeRapsim',format='E',unit='',array=self.misoeRapsim)
        c61=pyfits.Column(name='misoI',format='E',unit='mag',array=self.misoI)
        c62=pyfits.Column(name='misoeIapsim',format='E',unit='',array=self.misoeIapsim)
        c63=pyfits.Column(name='misoJ',format='E',unit='mag',array=self.misoJ)
        c64=pyfits.Column(name='misoeJapsim',format='E',unit='',array=self.misoeJapsim)
        c65=pyfits.Column(name='misoK',format='E',unit='mag',array=self.misoK)
        c66=pyfits.Column(name='misoeKapsim',format='E',unit='',array=self.misoeKapsim)
        c67=pyfits.Column(name='magV',format='E',unit='mag',array=self.magV)
        c68=pyfits.Column(name='mageVapsim',format='E',unit='',array=self.mageVapsim)
        c69=pyfits.Column(name='magR',format='E',unit='mag',array=self.magR)
        c70=pyfits.Column(name='mageRapsim',format='E',unit='',array=self.mageRapsim)
        c71=pyfits.Column(name='magI',format='E',unit='mag',array=self.magI)
        c72=pyfits.Column(name='mageIapsim',format='E',unit='',array=self.mageIapsim)
        c73=pyfits.Column(name='magJ',format='E',unit='mag',array=self.magJ      )
        c74=pyfits.Column(name='mageJapsim',format='E',unit='',array=self.mageJapsim)
        c75=pyfits.Column(name='magK',format='E',unit='mag',array=self.magK      )
        c76=pyfits.Column(name='mageKapsim',format='E',unit='',array=self.mageKapsim)
        c77=pyfits.Column(name='membflag',format='L',unit='',array=self.membflag  )
        c78=pyfits.Column(name='newspecmatchflag',format='I',unit='',array=self.newspecmatchflag   )
        c79=pyfits.Column(name='defmembflag',format='L',unit='',array=self.defmembflag)
        c80=pyfits.Column(name='photmembflag',format='L',unit='',array=self.photmembflag)
        c81=pyfits.Column(name='supermembflag',format='L',unit='',array=self.supermembflag)      
        c82=pyfits.Column(name='specz',format='E',unit='',array=self.specz)
        c83=pyfits.Column(name='spectype',format='I',unit='',array=self.spectype)  
        c84=pyfits.Column(name='specEWOII',format='E',unit='',array=self.specEWOII)   
        c85=pyfits.Column(name='matchflagspecediscs',format='I',unit='',array=self.matchflagspecediscs)
        c86=pyfits.Column(name='specEWOIIflag',format='I',unit='',array=self.specEWOIIflag)
        c87=pyfits.Column(name='bestz',format='E',unit='',array=self.bestz)
        c88=pyfits.Column(name='lowz',format='E',unit='',array=self.lowz)     
        c89=pyfits.Column(name='highz',format='E',unit='',array=self.highz)     
        c90=pyfits.Column(name='wmin',format='E',unit='',array=self.wmin)   
        c91=pyfits.Column(name='Pclust',format='E',unit='',array=self.Pclust)    
        c92=pyfits.Column(name='MR',format='D',unit='',array=self.MR)          
        c93=pyfits.Column(name='MU',format='D',unit='',array=self.MU)         
        c94=pyfits.Column(name='MV',format='D',unit='',array=self.MV)          
        c95=pyfits.Column(name='MB',format='D',unit='',array=self.MB)          
        c96=pyfits.Column(name='stellmass',format='D',unit='',array=self.stellmass)
        c97=pyfits.Column(name='redflag',format='I',unit='',array=self.redflag)
        c98=pyfits.Column(name='LUlowzclust',format='D',unit='',array=self.LUlowzclust)
        c99=pyfits.Column(name='LUbestzclust',format='D',unit='',array=self.LUbestzclust)
        c100=pyfits.Column(name='LUhighzclust',format='D',unit='',array=self.LUhighzclust) 
        c101=pyfits.Column(name='LBlowzclust',format='D',unit='',array=self.LBlowzclust)
        c102=pyfits.Column(name='LBbestzclust',format='D',unit='',array=self.LBbestzclust) 
        c103=pyfits.Column(name='LBhighzclust',format='D',unit='',array=self.LBhighzclust) 
        c104=pyfits.Column(name='LVlowzclust ',format='D',unit='',array=self.LVlowzclust) 
        c105=pyfits.Column(name='LVbestzclust',format='D',unit='',array=self.LVbestzclust) 
        c106=pyfits.Column(name='LVhighzclust',format='D',unit='',array=self.LVhighzclust) 
        c107=pyfits.Column(name='LRlowzclust',format='D',unit='',array=self.LRlowzclust)
        c108=pyfits.Column(name='LRbestzclust',format='D',unit='',array=self.LRbestzclust) 
        c109=pyfits.Column(name='LRhighzclust',format='D',unit='',array=self.LRhighzclust) 
        c110=pyfits.Column(name='LIlowzclust',format='D',unit='',array=self.LIlowzclust)
        c111=pyfits.Column(name='LIbestzclust',format='D',unit='',array=self.LIbestzclust) 
        c112=pyfits.Column(name='LIhighzclust',format='D',unit='',array=self.LIhighzclust) 
        c113=pyfits.Column(name='LJlowzclust',format='D',unit='',array=self.LJlowzclust)
        c114=pyfits.Column(name='LJbestzclust',format='D',unit='',array=self.LJbestzclust) 
        c115=pyfits.Column(name='LJhighzclust',format='D',unit='',array=self.LJhighzclust) 
        c116=pyfits.Column(name='LKlowzclust',format='D',unit='',array=self.LKlowzclust)
        c117=pyfits.Column(name='LKbestzclust',format='D',unit='',array=self.LKbestzclust) 
        c118=pyfits.Column(name='LKhighzclust',format='D',unit='',array=self.LKhighzclust)           
        
	mastertb=pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c12a,c13,c18,c18a,c19,c20,c21,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,c41,c42,c43,c44,c45,c46,c47,c48,c49,c50,c51,c52,c53,c54,c55,c56,c57,c58,c59,c60,c61,c62,c63,c64,c65,c66,c67,c68,c69,c70,c71,c72,c73,c74,c75,c76,c77,c78,c79,c80,c81,c82,c83,c84,c85,c86,c87,c88,c89,c90,c91,c92,c93,c94,c95,c96,c97,c98,c99,c100,c101,c102,c103,c104,c105,c106,c107,c108,c109,c110,c111,c112,c113,c114,c115,c116,c117,c118])#,c119,c120,c121,c122,c123,c124,c125,c126,c127,c128,c129,c130])
	s=homedir+'research/clusters/spitzer/MasterTables/'+self.fullname+'mastertable.fits'
	mastertb.writeto(s,clobber='yes')


    def othercalculations(self): #keeping some of the other code from readmaster() in case I need it

	#create field sample
	#zmin=0.43
	self.fieldmatchflag24=compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.matchflag24)
	self.fieldz=compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1)  & (self.matchflag24 > -0.1),self.specz)
	self.fieldf24c=compress((abs(self.newspecmatchflag) < 0.1) & (abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.flux24)
	self.fielderrf24c=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.flux24err)
	self.fieldMr=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.MRbestz)
	self.fieldMv=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.MVbestz)
	self.fieldBV=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.BVbest)
	self.fieldUB=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.UBbest)
	self.fieldra=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.ra)
	self.fielddec=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.dec)


	index=arange(len(self.ediscsID))
	self.fieldindex=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),index)

	ids=[]
	for i in index:
		ids.append(self.ediscsID[i])
	self.fieldediscsID=ids
	self.fieldLir=zeros(len(self.fieldz),'d')
	self.fielderrLir=zeros(len(self.fieldz),'d')
	for i in range(len(self.fieldz)):
		f24conv=splev(self.fieldz[i],tck24)
		conversion=1.e-6*1.e-23*(4.*pi*my.dLcm(self.fieldz[i],h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s

		self.fieldLir[i]=self.fieldf24c[i]*conversion*f24conv/Lsol
		self.fielderrLir[i]=self.fielderrf24c[i]*conversion*f24conv/Lsol
		dL=my.dL(self.fieldz[i],.7)

	self.fieldSFRir=self.fieldLir*bellconv*Lsol
	self.fielderrSFRir=self.fielderrLir*bellconv*Lsol
	print "Number of spec field galaxies = ",len(self.fieldz)
	#print "w/redshifts = ",self.fieldz


	nmultimatch=compress((self.matchflag24 > 0.1) & (self.nmatchediscs24 > 1.1),self.nmatchediscs24)
	n24match=len(compress(self.matchflag24 > 0.1,self.matchflag24))
	print "Number of ediscs sources w/more than 1 optical match = ",len(nmultimatch),n24match
	print "number of matches = ",nmultimatch
	self.nmultimatch=len(nmultimatch)
	self.n24match=n24match


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

							output112.write(string1)

		output111.close()
		output112.close()



	self.drdeg=sqrt(((self.ra-self.rac)*cos(self.dec*pi/180.))**2+(self.dec-self.decc)**2)#dr, projected distance from cluster center, in degrees
	self.drarcmin=self.drdeg*60.#dr in arcmin
	self.dr=self.drarcmin/self.r200arcmin#dr/ R200

	self.fielddr=compress((abs(self.newspecmatchflag) < 0.1) &(abs(self.specz-self.zcl) < 0.1) & (self.matchflag24 > -0.1),self.dr)

	#calculate fraction of OII members w/no 24 micron
	n=len(compress((self.newspecmatchflag > 0.1)& ((self.spectype-3) < 1.1)& (abs(self.matchflag24) < .1),self.matchflag24))
	nt=len(compress((self.newspecmatchflag > 0.1)& ((self.spectype-3.)<1.1) & (self.matchflag24 > -.5),self.matchflag24))
	print n,nt, "in calculating OII fraction"
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb OII w/no 24 = %i, Numb OII = %i, frac = %6.3f - %6.3f + %6.3f"%(n,nt,a,b,c)
	self.nOIIno24=n
	self.nOII=nt

	#calculate fraction of 24 micron - detected members w/no OII
	n=len(compress((self.newspecmatchflag > 0.1)& ((self.spectype-3.)<1.1) & (self.matchflag24 > .1),self.matchflag24))
	nt=len(compress((self.newspecmatchflag > 0.1)& (self.matchflag24 > 0.1),self.matchflag24))
	print "printing n,nt ",n,nt
	(a,b,c) = my.ratioerror(n,nt)
	print "Numb 24 spec memb w/emission lines = %i, Numb 24 = %i, frac = %6.3f - %6.3f +%6.3f"%(n,nt,a,b,c)
	self.n24specnoOII=n
	self.n24spec=nt
	self.nspec=len(compress((self.newspecmatchflag > 0.1),self.matchflag24))
	     

	for k in range(len(self.matchflag24)):
		if ((self.newspecmatchflag[k] > 0.1)& (self.specEWOIIflag[k] < 0.1)& (self.matchflag24[k] > .1)):
			s='%s \n'%(self.ediscsID[k])
			vandanaout.write(s)


def doall():
    cl1216=cluster('cl1216')
    cl1354=cluster('cl1354')    
    cl105412=cluster('cl105412') 
    cl1040 =cluster('cl1040')   
    cl105411=cluster('cl105411') 
    cl1227=cluster('cl1227')   
    cl1353=cluster('cl1353')   
    cl1037=cluster('cl1037')   
    cl1232=cluster('cl1232')   
    cl1411=cluster('cl1411')   
    cl1420=cluster('cl1420')   
    cl1301=cluster('cl1301')   
    cl1138=cluster('cl1138')   
    cl1018=cluster('cl1018')   
    cl1059=cluster('cl1059')   
    cl1202=cluster('cl1202')   

# ediscsID ra dec xcorr ycorr starflag EW EWerr SFR SFRerr matchflagha SFflag gimtype matchflagmorphgimtype vistype matchflagvistype matchflag24 flux24 flux24err nmatchediscs24 misoV misoeVapsim misoR misoeRapsim misoI misoeIapsim misoJ misoeJapsim misoK misoeKapsim magV mageVapsim magR mageRapsim magI mageIapsim magJ mageJapsim magK mageKapsim membflag newspecmatchflag defmembflag specz spectype specEWOII matchflagspecediscs specEWOIIflag bestz lowz highz wmin Pclust LUlowzclust LUbestzclust LUhighzclust LBlowzclust LBbestzclust LBhighzclust LVlowzclust LVbestzclust LVhighzclust LRlowzclust LRbestzclust LRhighzclust LIlowzclust LIbestzclust LIhighzclust LJlowzclust LJbestzclust LJhighzclust LKlowzclust LKbestzclust LKhighzclust fluxK UBlowzclust UBbestzclust UBhighzclust BVlowzclust BVbestzclust BVhighzclust UVlowzclust UVbestzclust UVhighzclust matchflagediscsirac iracf1 iracf2 iracf3 iracf4 erriracf1 erriracf2 erriracf3 erriracf4 iracsexflag0 iracsexflag1 iracwch1 iracwch2 iracwch3 iracwch4 iracwmin nmatchediscsirac L24 errL24 LHa errLHa snr24 imagex24 imagey24 fap1 fap2 fap3 fap4 fap5 fap6 fap7 fap8 fap9 fap10 errfap1 errfap2 errfap3 errfap4 errfap5 errfap6 errfap7 errfap8 errfap9 errfap10 LVlow   LVbest   LVhigh   LRlow   LRbest   LRhigh  UBlow   UBbest   UBhigh   BVlow   BVbest   BVhigh   UVlow   UVbest   UVhigh ra24 dec24 flux80flag ediscsIDold




'''
fluxK
UBlowzclust
UBbestzclust
UBhighzclust
BVlowzclust
BVbestzclust
BVhighzclust
UVlowzclust
UVbestzclust
UVhighzclust
LVlow
LVbest
LVhigh
LRlow
LRbest
LRhigh
UBlow
UBbest
UBhigh
BVlow
BVbest
BVhigh
UVlow
UVbest
UVhigh

fap1
fap2
fap3
fap4
fap5
fap6
fap7
fap8
fap9
fap10
errfap1
errfap2
errfap3
errfap4
errfap5
errfap6
errfap7
errfap8
errfap9
errfap10
'''



























