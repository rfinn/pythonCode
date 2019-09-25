#!/usr/bin/env python

import numpy as np
from pylab import *
from astropy.cosmology import WMAP9 as cosmo

dale_directory = '/Users/rfinn/research/dale/'
class models:
    def __init__(self):
        print "model SED templates"

    def readDaleSED(self):
        file=dale_directory+'spectra.dat'
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
        
        self.lam = np.zeros(ngal,'d')
        self.lamerr = np.zeros(ngal,'d')
        self.flux  = np.zeros([64,ngal],'d')
        self.alpha=np.zeros(64,'d')
        self.ircolor=np.zeros(64,'d')
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
                self.flux[j][i] = float(t[int(j+1)]) # in units of nu*f_nu
            i=i+1
        input.close()

        input=open(dale_directory+'alpha.dat','r')
        i=0
        for line in input:
            t=line.split()
            self.alpha[i]=float(t[0])
            self.ircolor[i]=float(t[1])
            i=i+1
        input.close()


        #calculate LTIR
        self.Ltir=np.zeros(len(self.alpha),'d')
        self.sfr=np.zeros(len(self.alpha),'d')
        self.logfnu24=np.zeros(len(self.alpha),'d')
        input=open(dale_directory+'model.0.dat','r')
        (c1,c2,c3)=(1.559,0.7686,1.347) # dale coefficients for z = 0
        (nu1,nu2,nu3)=(3.e8/(23.8e-6),3.e8/(70.e-6),3.e8/(160.e-6)) # frequency of spitzer bands
        #print nu1,nu2,nu3
        i=0
        for line in input:
            t=line.split()
            self.logfnu24[i]=float(t[11])
            self.Ltir[i]=(c1*nu1*10.**float(t[11])+c2*nu2*10.**float(t[12])+c3*nu3*10.**float(t[13]))*1.e-23#luminosity in erg/s
            self.Ltir[i]=self.Ltir[i]/3.9e33#luminosity in solar luminosities
            self.sfr[i]=self.Ltir[i]/(6.31e9)#convert to SFR
            i += 1
        
    def findIndexWavelength(self,wavelength,z):#find index of Dale model that is closest to observed wavelength
        '''
        wavelength = observed wavelength in microns
        z = redshift of source
        '''
        lobs=wavelength/(1.+z)
        min=100000.
        self.iwave=-99
        for i in range(len(self.lam)):
            d=abs(self.lam[i]-lobs)
            if d < min:
                min=d
                self.iwave=i  

    def plotSEDtemplate(self,flux_offset = 0):
        lam=np.array(self.lam,'f')
        flux1=np.array(self.flux[18][:],'f')
        flux2=np.array(self.flux[20][:],'f')
        flux3=np.array(self.flux[25][:],'f')
        flux4=np.array(self.flux[30][:],'f')
        flux5=np.array(self.flux[35][:],'f')

        plt.figure()
        plt.plot(lam,flux1 + flux_offset, 'b-')
        plt.plot(lam,flux2 + flux_offset, 'g-')
        plt.plot(lam,flux3 + flux_offset, 'r-')
        plt.plot(lam,flux4 + flux_offset, 'm-')
        plt.plot(lam,flux5 + flux_offset, 'k-')
        ax=plt.gca()
        ax.set_xscale("log")
        plt.xlabel(r'$\rm{Wavelength \ (\mu m)}$',fontsize=20)
        plt.ylabel(r'$\rm{\nu F_\nu} \ (Jy \ Hz)}$',fontsize=20)
        #plt.title('Dust-only Models')

        #plt.show()


    def scale_models(self,flux): #flux in Jy to scale to
        
        modindex=np.array([0,15,20,25,30],'i')
        model=np.zeros([len(modindex),len(self.flux[0][:])],'d')
        for j in range(len(modindex)):
            model[j][:]=self.flux[modindex[j]][:]
        self.iwave=int(self.iwave)
        scale=np.zeros(len(modindex),'d')
        sfr=np.zeros(len(modindex),'d')
        Ltir=np.zeros(len(modindex),'d')
        #distance_conversion = 4*np.pi*(cosmo.luminosity_distance(z).value*3.08e24)**2*1.e-23
        #print distance_conversion
        for j in range(len(scale)):
            scale[j]=np.log10(flux)-model[j][self.iwave]
            print j,scale[j],model[j][self.iwave]
            model[j][:]=10.**(model[j][:] + scale[j])
            sfr[j]=self.sfr[modindex[j]]*10.**(scale[j]-3.)#*distance_conversion#*4*np.pi*(cosmo.luminosity_distance(z).value*3.08e24)**2 # no idea what -3 is there for, maybe total LIR vs 3-800 microns?
            #print 'sfr = ',sfr[j],self.sfr[modindex[j]],10.**(scale[j]-3.)
            Ltir[j]=np.log10(self.Ltir[modindex[j]]*10.**(scale[j]-3.)*3.9e33)#*distance_conversion#*4*np.pi*(cosmo.luminosity_distance(z).value*3.08e24)**2
        (line1,line2,line3,line4,line5)=plt.plot(self.lam,model[0][:], 'y-',self.lam,model[1][:],'g-',self.lam,model[2][:], 'r-',self.lam,model[3][:], 'm-',self.lam,model[4][:], 'c-')
        s=[]
        for j in range(len(modindex)):
	    #print str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]])
	    #l=str(mod.alpha[modindex[j]])+", %3.2e"%(mod.sfr[modindex[j]]*10.**scale)
            #l="%3.1f, %5.2e, %4.1e"%(self.alpha[modindex[j]],Ltir[j],sfr[j])
            l="SFR = %5.3f"%(sfr[j])
            s.append(l)
        plt.legend((line1,line2,line3,line4,line5),(s),loc=(0.62,0.05),)
if __name__ == "__main__":
    mod = models()
    mod.findIndexWavelength(wavelength=12.,z=0)
    mod.scale_models()
