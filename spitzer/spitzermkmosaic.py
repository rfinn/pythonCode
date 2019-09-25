#!/usr/bin/env python
"""
run from /Users/rfinn/clusters/spitzer/GroupLirgs

"""

from pylab import *
from pyraf import iraf
import pyfits
import glob
from matplotlib.backends.backend_pdf import PdfFile

delta=30.#width of cutouts in arcsec

class cluster:
    def __init__(self,ncl):
	self.prefix=prefix[ncl]
	self.image=images[ncl]
	self.image24=images24[ncl]
	self.rimage24=rimages24[ncl]

    def readmastertable(self):
	s='/Users/rfinn/clusters/spitzer/MasterTables/'+self.prefix+'*mastertable.dat'
	file=glob.glob(
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


	self.specmembflag = zeros(ngal,'f')
	self.matchflagspec = zeros(ngal,'f')
	self.fflat24 = zeros(ngal,'f')
	self.errfflat24 = zeros(ngal,'f')
	self.matchflagediscsflat24 = zeros(ngal,'f')
	

        self.ediscsID = []
	self.ra = zeros(ngal,'f')
	self.dec = zeros(ngal,'f')
	self.xcorr = zeros(ngal,'f')
	self.ycorr = zeros(ngal,'f')
	self.starflag = zeros(ngal,'f')
	self.EW = zeros(ngal,'f')#Halpha EW
	self.EWerr = zeros(ngal,'f')
	self.SFR = zeros(ngal,'f')#Halpha SFR
	self.SFRerr = zeros(ngal,'f')
	self.matchflagha = zeros(ngal,'f')#match flag ha
	self.SFflag = zeros(ngal,'f')#if SF 
	self.gimtype = zeros(ngal,'f')
	self.matchflagmorphgimtype = zeros(ngal,'f')	
	self.vistype = zeros(ngal,'f')
	self.matchflagvistype = zeros(ngal,'f')
	self.matchflag24 = zeros(ngal,'f')
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
	self.membflag = zeros(ngal,'f')
	self.newspecmatchflag = zeros(ngal,'f')	
	self.defmembflag = zeros(ngal,'f')
	self.specz = zeros(ngal,'f')
	self.spectype = zeros(ngal,'f')
	self.specEWOII = zeros(ngal,'f')
	self.matchflagspecediscs = zeros(ngal,'f')
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
	self.matchflagediscsirac = zeros(ngal,'f')
	self.iracf1 = zeros(ngal,'f')
	self.iracf2 = zeros(ngal,'f')
	self.iracf3 = zeros(ngal,'f')
	self.iracf4 = zeros(ngal,'f')
	self.erriracf1 = zeros(ngal,'f')
	self.erriracf2 = zeros(ngal,'f')
	self.erriracf3 = zeros(ngal,'f')
	self.erriracf4 = zeros(ngal,'f')
	self.iracsexflag0 = zeros(ngal,'f')
	self.iracsexflag1 = zeros(ngal,'f')
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



	    (self.ra[i],self.dec[i],self.xcorr[i],self.ycorr[i],self.starflag[i],self.EW[i],self.EWerr[i],self.SFR[i],self.SFRerr[i],self.matchflagha[i],self.SFflag[i],self.gimtype[i],self.matchflagmorphgimtype[i],self.vistype[i],self.matchflagvistype[i],self.matchflag24[i],self.flux24[i],self.flux24err[i],self.nmatchediscs24[i],self.misoV[i],self.misoeVapsim[i],self.misoR[i],self.misoeRapsim[i],self.misoI[i],self.misoeIapsim[i],self.misoJ[i],self.misoeJapsim[i],self.misoK[i],self.misoeKapsim[i],self.magV[i],self.mageVapsim[i],self.magR[i],self.mageRapsim[i],self.magI[i],self.mageIapsim[i],self.magJ[i],self.mageJapsim[i],self.magK[i],self.mageKapsim[i],self.membflag[i],self.newspecmatchflag[i],self.defmembflag[i],self.specz[i],self.spectype[i],self.specEWOII[i],self.matchflagspecediscs[i],self.specEWOIIflag[i],self.bestz[i],self.lowz[i],self.highz[i],self.wmin[i],self.Pclust[i],self.LUlowzclust[i],self.LUbestzclust[i],self.LUhighzclust[i],self.LBlowzclust[i],self.LBbestzclust[i],self.LBhighzclust[i],self.LVlowzclust[i],self.LVbestzclust[i],self.LVhighzclust[i],self.LRlowzclust[i],self.LRbestzclust[i],self.LRhighzclust[i],self.LIlowzclust[i],self.LIbestzclust[i],self.LIhighzclust[i],self.LJlowzclust[i],self.LJbestzclust[i],self.LJhighzclust[i],self.LKlowzclust[i],self.LKbestzclust[i],self.LKhighzclust[i],self.fluxK[i],self.UBlowzclust[i],self.UBbestzclust[i],self.UBhighzclust[i],self.BVlowzclust[i],self.BVbestzclust[i],self.BVhighzclust[i],self.UVlowzclust[i],self.UVbestzclust[i],self.UVhighzclust[i],self.matchflagediscsirac[i],self.iracf1[i],self.iracf2[i],self.iracf3[i],self.iracf4[i],self.erriracf1[i],self.erriracf2[i],self.erriracf3[i],self.erriracf4[i],self.iracsexflag0[i],self.iracsexflag1[i],self.iracwch1[i],self.iracwch2[i],self.iracwch3[i],self.iracwch4[i],self.iracwmin[i],self.nmatchediscsirac[i],self.L24[i],self.errL24[i],self.LHa[i],self.errLHa[i],self.snr24[i],self.imagex24[i],self.imagey24[i],self.fap1[i],self.fap2[i],self.fap3[i],self.fap4[i],self.fap5[i],self.fap6[i],self.fap7[i],self.fap8[i],self.fap9[i],self.fap10[i],self.errfap1[i],self.errfap2[i],self.errfap3[i],self.errfap4[i],self.errfap5[i],self.errfap6[i],self.errfap7[i],self.errfap8[i],self.errfap9[i],self.errfap10[i],self.LVlow[i],  self.LVbest[i],  self.LVhigh[i],  self.LRlow[i],  self.LRbest[i],  self.LRhigh[i], self.UBlow[i],  self.UBbest[i],  self.UBhigh[i],  self.BVlow[i],  self.BVbest[i],  self.BVhigh[i],  self.UVlow[i],  self.UVbest[i],  self.UVhigh[i],self.ra24[i],self.dec24[i])=t[1:]
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


	self.MR = 4.28 - 2.5*log10(self.LRbestzclust)-25.
	#self.errMRlo = self.MR - (4.28 - 2.5*log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*log10(self.Rlumhi) -25.) - self.MR 
	self.MV = 4.82 - 2.5*log10(self.LVbestzclust) -25.
	#self.errMvlo = self.MV - (4.82 - 2.5*log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*log10(self.vlumhi) - 25. ) - self.MV
	self.MB = 5.48 - 2.5*log10(self.LBbestzclust) -25.#check zeropoint here...

	self.MRbestz = 4.28 - 2.5*log10(self.LRbest)-25.
	#self.errMRlo = self.MR - (4.28 - 2.5*log10(self.Rlumlo)-25.)
	#self.errMRhi = (4.28 - 2.5*log10(self.Rlumhi) -25.) - self.MR 
	self.MVbestz = 4.82 - 2.5*log10(self.LVbest) -25.
	#self.errMvlo = self.MV - (4.82 - 2.5*log10(self.vlumlo) -25. )
	#self.errMvhi = (4.82 - 2.5*log10(self.vlumhi) - 25. ) - self.MV

	
	self.L24=self.L24/Lsol
	self.errL24=self.errL24/Lsol
	#self.f24c=self.fap4*self.apcor4#corrected to infinite aperture
	#self.errf24c=self.errfap4*self.apcor4
	self.f24c=self.flux24
	self.errf24c=self.flux24err
	conversion=1.e-6*1.e-23*(4.*pi*my.dLcm(self.zcl,h)**2)*(3.e8/23.8e-6)#converts from microJy to erg/s
	Lir=self.f24c*conversion*self.f24IRconv#factor of ten converts get scaling to go from from nu fnu(24) to L(IR)
	self.Lir=array(Lir,'d')/Lsol#convert to solar luminosities
	errLir=self.errf24c*conversion*self.f24IRconv
	self.errLir=array(errLir,'d')/Lsol

	self.Lir80=self.f80*conversion*self.f24IRconv/Lsol#luminosity corresponding to f at 80% comple
	print self.prefix," log10(Lir80) = %6.3f" %(log10(self.Lir80))

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

	#create field sample
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

    def getcoords(self):
	f='/Users/rfinn/clusters/spitzer/MasterTables/'+self.prefix+'.HST.24em.spec.dat'
	infile=open(f,'r')
	ra=[]
	dec=[]
	id=[]
	specmembflag=[]
	mipsflag=[]
	f24=[]
	SFRir=[]
	specz=[]
	f80flag=[]
	for line in infile:
	    t=line.split()
	    ra.append(float(t[0]))
	    dec.append(float(t[1]))
	    id.append(t[2])
	    mipsflag.append(float(t[3]))
	    specmembflag.append(float(t[4]))
	    specz.append(float(t[5]))
	    f24.append(float(t[6]))
	    SFRir.append(float(t[7]))
	    f80flag.append(float(t[8]))
			    
	infile.close()
	self.ra=array(ra,'f')
	self.dec=array(dec,'f')
	self.ediscsID=id
	self.mipsflag=array(mipsflag,'f')
	self.specmembflag=array(specmembflag,'f')
	self.specz=array(specz,'f')
	self.f24=array(f24,'f')
	self.SFRir=array(SFRir,'f')

	s=self.prefix+'.radec'
	self.incoords=s
	out1=open(s,'w')
	for i in range(len(self.ediscsID)):
	    name=id[i]
	    s='%f %f %s\n'%(self.ra[i],self.dec[i],self.ediscsID[i])
	    out1.write(s)
	out1.close()

    def transcoords(self):
	print ncl,self.prefix
	outcoords=str(self.prefix)+'.xy'
	iraf.imcoords.wcsctran(image=self.image,input=self.incoords,output=outcoords,inwcs='world',outwcs='logical')
	self.outcoords=outcoords


    def makecutouts(self):
	cutouts=[]
	iraf.imgets(image=self.image,param='CD2_1')#get x value corresponding to RA 
	xplate=abs(float(iraf.imgets.value))#deg/pixel
	dpix=delta/3600./xplate/2.
	
	infile=open(self.outcoords,'r')
	for line in infile:
	    if line.find('#') > -1:
		continue
	    if len(line)<2:
		continue
	    x,y,id=line.split()
	    x=float(x)
	    y=float(y)

	    xmin=int(round(x-dpix))
	    xmax=int(round(x+dpix))
	    ymin=int(round(y-dpix))
	    ymax=int(round(y+dpix))
	    s=self.image+'[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax)
	    print s
		
	    outim=id+'cutout.fits'
	    print outim
	    if (makecutouts > 0.1):
		iraf.imcopy(s,outim)
	    cutouts.append(outim)
	infile.close()
	self.cutouts=cutouts



    def transcoords24(self):
	outcoords=str(self.prefix)+'.xy24'
	iraf.imcoords.wcsctran(image=self.rimage24,input=self.incoords,output=outcoords,inwcs='world',outwcs='logical')
	self.outcoords24=outcoords

    def makecutouts24(self):
	cutouts24=[]
	cutout24flag=ones(len(self.ra))
	iraf.imgets(image=self.rimage24,param='CD2_1')#get x plate scale on rotated image
	print iraf.imgets.value
	xplate=abs(float(iraf.imgets.value))#deg/pixel
	dpix=delta/3600./xplate/2.
	iraf.imgets(image=self.rimage24,param='naxis1')#get x value corresponding to RA 
	xpixmax=(int(iraf.imgets.value))#deg/pixel
	iraf.imgets(image=self.rimage24,param='naxis2')#get x value corresponding to RA 
	ypixmax=(int(iraf.imgets.value))#deg/pixel
	infile=open(self.outcoords24,'r')
	#print outcoords
        i=-1
	for line in infile:

	    #print images24[i],line,outcoords
	    if line.find('#') > -1:
		continue
	    if len(line)<2:
		continue
            i=i+1
            print i,line,self.ediscsID[i]
	    x,y,id=line.split()
	    x=float(x)
	    y=float(y)
	    #delta=100.
	    xmin=int(round(x-dpix))
	    xmax=int(round(x+dpix))
	    ymin=int(round(y-dpix))
	    ymax=int(round(y+dpix))
	    if xmin < 1:
		xmin=1
	    if ymin < 1:
		ymin=1
	    if xmax > xpixmax:
		xmax=xpixmax
	    if ymax > ypixmax:
		ymax=ypixmax
	    if xmin > xpixmax:
		cutout24flag[i]=0
		print self.ediscsID[i],"pixel value out of range 1"
		continue
	    if ymin > ypixmax:
		cutout24flag[i]=0
		print self.ediscsID[i],"pixel value out of range 2"
		continue
	    if xmax < 5:
		cutout24flag[i]=0
		print self.ediscsID[i]," pixel value out of range 3"
		continue
	    if ymax < 5:
		cutout24flag[i]=0
		print self.ediscsID[i]," pixel value out of range 4"
		continue
	    if (abs(xmax-xmin)<5.):
		cutout24flag[i]=0
		print self.ediscsID[i]," image less than 5 pixels wide"
		continue
	    if (abs(ymax-ymin)<5.):
		cutout24flag[i]=0
		print self.ediscsID[i]," image height less than 5 pixels"
		continue

	    s=self.rimage24+'[%i:%i,%i:%i]'%(xmin,xmax,ymin,ymax)
	    print s
	    #print ra[i],dec[i]
	    outim=id+'cutout24.fits'
	    #print outim
	    
	    try:
		iraf.imcopy(s,outim)
	    except:
		cutout24flag[i]=0
		print self.ediscsID[i]," pixel value out of range 7"
		continue

	    cutouts24.append(outim)
	infile.close()
	self.cutouts24=cutouts24
	self.cutout24flag=cutout24flag
        f=self.prefix+'cutoutfiles.dat'
        output99=open(f,'w')
        for i in range(len(cutout24flag)):
            s='%s %i \n'%(self.ediscsID[i],self.cutout24flag[i])
            output99.write(s)
        output99.close()
    def displaycutouts(self):
        hstfiles=[]
        mipsfiles=[]
        f=self.prefix+'cutoutfiles.dat'
        infile=open(f,'r')
	cutout24flag=[]
        for line in infile:
            t=line.split()
	    cutout24flag.append(float(t[1]))
	    id=t[0]
	    hstfiles.append(id+"cutout.fits")
	    mipsfiles.append(id+"cutout24.fits")
	#hstfiles=self.cutouts
	#mipsfiles=self.cutouts24
	print "number of cutouts = ",len(hstfiles),len(mipsfiles)
	for i in range(len(mipsfiles)):
	    if cutout24flag[i] > 0.1:
		print hstfiles[i],mipsfiles[i]
		print self.ediscsID[i]
		print 'mips flag = ',self.mipsflag[i]
		print '24um flux (SFSR) = ',self.f24[i],self.SFRir[i]
		print 'flux80 flag = ',self.f80flag[i]
		print 'spec memb flag (z) = ',self.specmembflag[i],self.specz[i]
		try:
		    os.system('xpaset -p ds9 frame 2')
		    s='xpaset -p ds9 file '+mipsfiles[i]
		    s='cat '+mipsfiles[i]+' | xpaset ds9 fits '
		    os.system(s)
		    os.system('xpaset -p ds9 frame 1')
		    s='xpaset -p ds9 file '+hstfiles[i]
		    s='cat '+hstfiles[i]+' | xpaset ds9 fits '
		    os.system(s)
		    os.system('xpaset -p ds9 match frames wcs')
		except:
		    print 'problem displaying ',hstfiles[i] 
		#iraf.display(hstfiles[i],frame=1,contrast=0.01,fill='No')
		#iraf.display(mipsfiles[i],frame=2,contrast=.01,fill='No',zscale='no',zrange='no',z1=-.01,z2=.3)
		t=raw_input('hit any key to continue, q to quit \n')
		print t
		t=str(t)
		if t.find('q') > -1:
		    break

    def plotcutouts(self):
	clf()
	cla()
	s=self.prefix+'Cutouts.pdf'
	#pdf=PdfFile('s')
	figure(figsize=(12,16))
	subplots_adjust(left=0.1, right=.95,bottom=.1,top=0.95,wspace=0.001,hspace=0.001)
	nx=6
	ny=8
	nplot=0
	npage=1
	nplotmax=float(nx)*float(ny)
	#hstfiles=glob.glob('EDCSNJ*cutout.fits')
	#mipsfiles=glob.glob('EDCSNJ*cutout24.fits')

        hstfiles=[]
        mipsfiles=[]
        f=self.prefix+'cutoutfiles.dat'
        infile=open(f,'r')
        for line in infile:
            t=line.split()
            if (float(t[1])>.1):#24um image is ok
                id=t[0]
                hstfiles.append(id+"cutout.fits")
                mipsfiles.append(id+"cutout24.fits")
	#hstfiles=self.cutouts
	#mipsfiles=self.cutouts24
	print "number of cutouts = ",len(hstfiles),len(mipsfiles)

	for i in range(len(mipsfiles)):
	    if (self.cutout24flag[i] < 0.5):
		continue
	    if (float(nplot)> (nplotmax-1)):
		f=self.prefix+'AllCutouts'+str(npage)+'.eps'
		npage=npage+1
		#savefig(pdf,format='pdf')
		savefig(f)
		print "starting new plot",npage
		clf()
		cla()
		nplot=0
	    name=hstfiles[i]
	    for j in range(2):
		nplot=nplot+1
		k=nplot
		print k,nplotmax,mipsfiles[i]
		subplot(ny,nx,nplot)
		if (j == 0):
		    fits=pyfits.open(hstfiles[i])
		    name=hstfiles[i]
		    t=name.split('cutout')
		    im=fits[0].data.copy()
		    lim=.07
		    im[where(im>(lim))]=lim
		    nave=3.5
		    im[where(im<0)]=0
		else:
		    fits=pyfits.open(mipsfiles[i])
		    name=mipsfiles[i]
		    im=fits[0].data.copy()
		    im[where(im<0)]=0
		    im[where(im>(.08))]=.08
		    axis([1.,5.,1.,5.])
	
		fits.close()
		axis('equal')
		imshow(-1.*(im),interpolation='nearest',origin='upper',cmap='gray')#,vmin=myvmin,vmax=myvmax)
		ax=gca()
		ax.set_xticklabels(([]))
		ax=gca()
		ax.set_yticklabels(([]))
		col='k'
	    #print col
		t=self.ediscsID[i]
		junk,id=t.split('J')

		text(.5,.82,id,fontsize=10,transform=ax.transAxes,color=col,horizontalalignment='center')
		if (j == 0):
		    if (self.specmembflag[i] > 0.1):
			s="%6.4f"%(self.specz[i])
			text(.6,.75,s,fontsize=10,transform=ax.transAxes,color=col,horizontalalignment='left')
		else:
		    col='r'
		    s="%2.1f"%(self.mipsflag[i])
		    text(.6,.75,s,fontsize=10,transform=ax.transAxes,color=col,horizontalalignment='left')

		    if (self.mipsflag[i] > 0.1):
			s="%5.1f, %6.1f"%(self.f24[i],self.SFRir[i])
			text(.6,.65,s,fontsize=10,transform=ax.transAxes,color=col,horizontalalignment='left')

	f=self.prefix+'AllCutouts'+str(npage)+'.eps'
	savefig(f)
	#savefig(pdf,format='pdf')
	#pdf.close()


    def doall(self):
	self.getcoords()
	#self.transcoords()
	#self.makecutouts()
	#self.transcoords24()
	#self.makecutouts24()
	#self.plotcutouts()
	self.displaycutouts()


def rotate24():

    #for i in range(len(images24)):
    for i in range(3,6):
    #for i in range(1):
	iraf.imgets(image=images24[i],param='CROTA2')#[deg] Orientation of axis 2 (W of N, +=CW)
	print iraf.imgets.value
	rot24=float(iraf.imgets.value)
	iraf.imgets(image=images[i],param='ORIENTAT')#position angle of image y axis (deg. e of n)
	print iraf.imgets.value
	rotacs=float(iraf.imgets.value)
	angle=-1*rot24-rotacs
	outimage=rimages24[i]
	iraf.rotate(input=images24[i],output=outimage,rotation=angle)



images=[]
images24=[]
rimages24=[]
i=0
images=['/Users/rfinn/clusters/ediscs/HST-images/cl1040_drz_gsci.fits',
'/Users/rfinn/clusters/ediscs/HST-images/cl1054-12_drz_gsci.fits',
'/Users/rfinn/clusters/ediscs/HST-images/cl1054-11_drz_gsci.fits',
#'/Users/rfinn/clusters/ediscs/HST-images/cl1103_drz_gsci.fits',
'/Users/rfinn/clusters/ediscs/HST-images/cl1216_drz_gsci.fits',
'/Users/rfinn/clusters/ediscs/HST-images/cl1232_drz_gsci.fits',
'/Users/rfinn/clusters/ediscs/HST-images/cl1354_drz_gsci.fits']
path24='/Users/rfinn/clusters/spitzer/final-images/'
prefix=['cl1040','cl105412','cl105411','cl1216','cl1232','cl1354']
for pre in prefix:
    s=path24+pre+'_final24.fits'
    images24.append(s)
    s=path24+'r'+pre+'_final24.fits'
    rimages24.append(s)
incoords=['cl1040lirgs.radec','cl1103lirgs.radec','cl1054-12lirgs.radec','cl1054-11lirgs.radec']
#prefix=['cl1040','cl1103','cl1054-12','cl1054-11']

#rotate24()
#transcoords()
#cutouts=makecutouts()
#transcoords24()
#cutouts24=makecutouts24(delta)
#plotcutouts()
makecutouts=1

ncl=0
cl1040=cluster(ncl)
cl1040.doall()

ncl=1
cl105411=cluster(ncl)
#cl105411.doall()

ncl=2
cl105412=cluster(ncl)
#cl105412.doall()

ncl=3
cl1216=cluster(ncl)
#cl1216.doall()
