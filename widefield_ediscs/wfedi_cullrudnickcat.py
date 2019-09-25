#!/usr/bin/env python
import atpy
import os
from ediscsCommon import *
import chary_elbaz_24um as chary
#from pylab import *
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

class catalog:
    def __init__(self):
        # infile=homedir+'research/NSA/nsa_v0_1_2.fits'
        infile=homedir+'research/WifiEdiscs/rudnick_catalog/megacat.v5.4.mips.fits'
        infile=homedir+'research/WifiEdiscs/rudnick_catalog/megacat.topcat.fits'
        self.ndat=atpy.Table(infile,type='fits')
    def clustersubset(self,clustername):
        # create a file containing the subset of NSA galaxies w/in 3 deg of cluster center
        prefix=clustername
        ra=racenter[clustername]
        dec=deccenter[clustername]
        zcluster=redshift[clustername]
        zmin=zcluster-.2
        zmax=zcluster+.2
        distance=sqrt((ra-self.ndat.RA)**2+(dec-self.ndat.Dec)**2)
        mipsflag=self.ndat.flux > 0.
        ldpflag=self.ndat.zLDP > 0.
        redshiftflag=(self.ndat.zLDP > zmin) & (self.ndat.zLDP < zmax)
        keepflag=(distance < catalog_radial_cut) & mipsflag & ldpflag & redshiftflag
        self.keepflag=keepflag
        clustersubset=self.ndat[keepflag]
        self.ndatsub=self.ndat.where(keepflag)

        self.ceLIR=zeros(len(keepflag),'d')
        self.ceSFR=zeros(len(keepflag),'d')
        self.ceLIR_ZCLUST=zeros(len(keepflag),'d')
        self.ceSFR_ZCLUST=zeros(len(keepflag),'d')
        # CALCULATE LIR AND SFR USING REDSHIFT
        self.ceLIR[keepflag],self.ceSFR[keepflag]=chary.chary_elbaz_24um(self.ndat.zLDP[keepflag],self.ndat.flux[keepflag])

        # CALCULATE LIR AND SFR ASSUMING THE GALAXY IS AT THE CLUSTER REDSHIFT
        cl_redshift=ones(len(keepflag),'f')*zcluster
        self.ceLIR_ZCLUST[keepflag],self.ceSFR_ZCLUST[keepflag]=chary.chary_elbaz_24um(cl_redshift[keepflag],self.ndat.flux[keepflag])#*mipsconv_MJysr_to_uJy)

        self.ndatsub.add_column('LIR_ZLDP',self.ceLIR[keepflag],unit='Lsun')
        self.ndatsub.add_column('SFRIR_ZLDP',self.ceSFR[keepflag],unit='Msun/yr')
        self.ndatsub.add_column('LIR_ZCLUST',self.ceLIR_ZCLUST[keepflag],unit='Lsun')
        self.ndatsub.add_column('SFRIR_ZCLUST',self.ceSFR_ZCLUST[keepflag],unit='Msun/yr')
        outfile=homedir+'research/WifiEdiscs/rudnick_catalog/'+clustername+'_megacat.v5.4.fits'
        print outfile, sum(keepflag)
        if os.path.exists(outfile):
            os.remove(outfile)
        self.ndatsub.write(outfile)



edi=catalog()        
catalog_radial_cut=14./60. # radial cut is 14 arcmin, bigger than VIMOS FOV
def getclustersubsets():
    for cl in fieldnames:
        edi.clustersubset(cl)
