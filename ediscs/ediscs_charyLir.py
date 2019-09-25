#!/usr/bin/env python

'''
    Written by Rose Finn, updated on July 2, 2013

    PURPOSE
      This program calculates the total IR luminosity and SFR for the ediscs spec sample
      using the Chary & Elbaz templates.  The program calculates two sets of 
      (LIR,SFR) using distances derived from the galaxy redshift and 
      the cluster redshift.  You should use the values for the cluster redshift
      for cluster members.

      The program also using the completeness estimates that I made 

      /Users/rfinn/research/clusters/spitzer/completeness

      to assign a weight to each galaxy.  
      
    CALLING SEQUENCE

      from the terminal command line, type:
        ediscs_charyLir.py

    INPUT PARAMETERS
      none
      
    OUTPUT PARAMETERS
      Creates a table
      
      outfile= filedir+'ediscs_charyelbaz_lir.fits'

      that contains a line-matched table of the two sets of (LIR,SFR), the weight, 
      and other key info from John's spec table


    PROCEDURE

      Completeness and weight are assigned to each galaxy by interpolating 
      the completeness vs flux data.

    REQUIRED PYTHON MODULES
      numpy
      pylab
      os
      atpy
      scipy

    ADDITIONAL REQUIRED MODULES
      ediscsCommon
      chary_elbaz_24um

    NOTES
      Required Files
      
        John Moustakas's ediscs spec files:
          filedir+'ediscs_spec1d_info.fits'
          filedir+'ediscs_photometry.v23.fits.gz'
        My completeness files:
          homedir+'research/clusters/spitzer/completeness/'+cl+'FracComplvsFlux.dat'

    UPDATES
      2015-02-15:
      - updating directories to reflect paths on my newest laptop;
      - adding measurements from 1103 field
      - reading 24um fluxes from ediscs_ircatalogs_spec.fits

'''

import os
import numpy as np
from scipy.interpolate import interp1d
#import pyfits
import atpy

from ediscsCommon import *
#from LCSReadmasterBaseWithProfileFits import *
import chary_elbaz_24um as chary

mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro or laptop"
    homedir='/Users/rfinn/'
    filedir=homedir+'research/Moustakas/'
    filedir=homedir+'Dropbox/research-macbook/ediscsClusters/'
    completeness_directory=homedir+'Dropbox/research-macbook/clusters/spitzer/completeness/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'
    filedir='/home/share/research/ediscs/'
    completeness_directory=homedir+'research/clusters/spitzer/completeness/'

# READ IN REDSHIFT FROM ediscs_spec1d_info.fits
infile=filedir+'ediscs_spec1d_info.fits'
spec=atpy.Table(infile)
# READ IN F24 FROM ediscs_photometry.v23.fits
infile=filedir+'ediscs_photometry.v23.fits.gz'
ephot=atpy.Table(infile)
infile=filedir+'ediscs_ircatalogs_spec.fits'
phot=atpy.Table(infile)

# track cluster/group vs field galaxies

clusterflag=spec.CLUSTER_Z > 0.



# CALCULATE THE CHARY & ELBAZ L_IR USING SPECTROSCOPIC REDSHIFT

lir_zspec=np.zeros(len(spec.Z),'d')
sfr_zspec=np.zeros(len(spec.Z),'d')
MATCHFLAG24=np.zeros(len(spec.Z),'bool')
flag=(phot.MATCHFLAG24 == 1)
print 'number of 24um matches = ',sum(flag)
MATCHFLAG24[flag]=np.ones(sum(flag))#,'bool')
on24imageflag=np.zeros(len(spec.Z),'bool')
flag=(phot.MATCHFLAG24 > -1)
on24imageflag[flag]=np.ones(sum(flag),'bool')
print 'calculating lir_zspec'
lir_zspec[MATCHFLAG24],sfr_zspec[MATCHFLAG24]=chary.chary_elbaz_24um(spec.Z[MATCHFLAG24],phot.F24[MATCHFLAG24])

print 'calculating lir_zclust'
# CALCULATE CHARY & ELBAZ L_IR USING CLUSTER REDSHIFT
lir_zclust=np.zeros(len(spec.Z),'d')
sfr_zclust=np.zeros(len(spec.Z),'d')
flag=clusterflag & MATCHFLAG24
lir_zclust[flag],sfr_zclust[flag]=chary.chary_elbaz_24um(spec.Z[flag],phot.F24[flag])




# ESTIMATE THE COMPLETENESS BY INTERPOLATING COMPLETENESS VS FLUX DATA FILES
#   NEED TO PROPERLY ASSOCIATE EACH GALAXY WITH THE CORRESPONDING 24-MICRON IMAGE
#   CLUSTER field contains the name of the target field

mips_weight=np.zeros(len(spec.Z),'f')
mips_compl=np.zeros(len(spec.Z),'f')
for cl in fieldnames:
    print cl
    # open data file containing completeness vs flux in micro-Jy
    infile=completeness_directory+cl+'FracComplvsFlux.dat'
    cdat=np.loadtxt(infile)
    flux=cdat[:,0] # flux in MJy/sr
    comp=cdat[:,1]
    comp=comp[0:-2] # late row in file is all zeros, so getting rid of that
    flux=flux[0:-2] # late row in file is all zeros, so getting rid of that


    # INTERPOLATE DATA
    # interpolate completeness for cl
    interpfunc=interp1d(flux,comp,bounds_error=False,fill_value=1.0) # save for later use w/in loops

    # identify galaxies with spec.CLUSTER == cl
    if len(cl) < 8:
        cl=cl+'  '
    membindex=np.where((spec.CLUSTER == cl) & MATCHFLAG24)
    # for galaxies in the cl field, calculate completeness
    t=interpfunc(phot.F24[membindex])
    mips_compl[membindex]=np.array(t,'f')

    #flag=phot.FLUX24[membindex] < min(flux)
    #mips_compl[membindex][flag]=zeros(len(flag),'f')
    #flag=phot.FLUX24[membindex] > max(flux)
    #mips_compl[membindex][flag]=ones(len(flag),'f')

    mips_weight[membindex]=1./mips_compl[membindex]


# WRITE OUT FILE CONTAINING EDISCS_ID, SPECZ, F24,F24_ERR LIR_ZSPEC, SFR_ZSPEC, LIR_ZCLUST, SFR_ZCLUST 
output=atpy.Table()
output.add_column('EDISCS_ID',spec.EDISCS_ID)
output.add_column('GALAXY',spec.GALAXY)
output.add_column('Z',spec.Z)
output.add_column('CLUSTER_MEMBER',clusterflag)
output.add_column('CLUSTER_Z',spec.CLUSTER_Z)
output.add_column('CLUSTER_FULLNAME',spec.CLUSTER_FULLNAME)
output.add_column('MATCHFLAG24',MATCHFLAG24,unit='')
output.add_column('ON24IMAGEFLAG',on24imageflag,unit='')
output.add_column('FLUX24',phot.F24*MATCHFLAG24,unit='micro-Jy')
output.add_column('FLUX24ERR',phot.F24_ERR*MATCHFLAG24,unit='micro-Jy')
#output.add_column('LIR_OLD',phot.LIR*MATCHFLAG24,unit='Lsun')
output.add_column('LOG_LIR_OLD',phot.LOG_LIR*MATCHFLAG24,unit='Lsun')
output.add_column('LIR_ZSPEC',lir_zspec,unit='Lsun')
output.add_column('SFR_ZSPEC',sfr_zspec,unit='Msun/yr')
output.add_column('LIR_ZCLUST',lir_zclust,unit='Lsun')
output.add_column('SFR_ZCLUST',sfr_zclust,unit='Msun/yr')
#output.add_column('MIPS_COMPLETENESS',mips_compl,unit='')
output.add_column('MIPS_WEIGHT',mips_compl,unit='')
output.add_column('SPEC_WEIGHT',ephot.SPEC_WEIGHT,unit='')
output.add_column('GEOM_WEIGHT',ephot.GEOM_WEIGHT,unit='')
#output.name=self.table_name

outfile= filedir+'ediscs_charyelbaz_lir.fits'
if os.path.exists(outfile):
    os.remove(outfile)
output.write(str(outfile))


