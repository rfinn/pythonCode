#!/usr/bin/env python

#  Written by Rose A. Finn, Feb 6, 2013
# 
# program reads in the agcnorth.sav file using idlsave
# and then cuts the agc according to RA, Dec, and redshift required
# for the Local Cluster Survey
# Table : /USERS/RFINN/RESEARCH/NSA/NSA_V0_1_2.FITS#1
# ----------------------------------------------------
# |            Name | Unit |           Type | Format |
# ----------------------------------------------------
# |         IAUNAME |      |           |S19 |  3689s |
# |          SUBDIR |      |           |S27 |  3689s |
# |              RA |      |            >f8 | 25.17e |
# |             DEC |      |            >f8 | 25.17e |
# |           ISDSS |      |            >i4 |    12i |
# |            INED |      |            >i4 |    12i |
# |          ISIXDF |      |            >i4 |    12i |
# |        IALFALFA |      |            >i4 |    12i |
# |           IZCAT |      |            >i4 |    12i |
# |          ITWODF |      |            >i4 |    12i |
# |             MAG |      |            >f4 |  16.8e |
# |               Z |      |            >f4 |  16.8e |
# |            ZSRC |      |            |S7 |  3689s |
# |            SIZE |      |            >f4 |  16.8e |
# |             RUN |      |            >i2 |     5i |
# |          CAMCOL |      |            >i2 |     5i |
# |           FIELD |      |            >i2 |     5i |
# |           RERUN |      |            |S3 |  3689s |
# |            XPOS |      |            >f4 |  16.8e |
# |            YPOS |      |            >f4 |  16.8e |
# |             ZLG |      |            >f4 |  16.8e |
# |           ZDIST |      |            >f4 |  16.8e |
# |       ZDIST_ERR |      |            >f4 |  16.8e |
# |           NSAID |      |            >i4 |    12i |
# |            NMGY |      |   ('>f4',(7,)) |  16.8e |
# |       NMGY_IVAR |      |   ('>f4',(7,)) |  16.8e |
# |              OK |      |            >i2 |     5i |
# |           RNMGY |      |   ('>f4',(7,)) |  16.8e |
# |          ABSMAG |      |   ('>f4',(7,)) |  16.8e |
# |          AMIVAR |      |   ('>f4',(7,)) |  16.8e |
# |      EXTINCTION |      |   ('>f4',(7,)) |  16.8e |
# |        KCORRECT |      |   ('>f4',(7,)) |  16.8e |
# |          KCOEFF |      |   ('>f4',(5,)) |  16.8e |
# |            MTOL |      |   ('>f4',(7,)) |  16.8e |
# |            B300 |      |            >f4 |  16.8e |
# |           B1000 |      |            >f4 |  16.8e |
# |            METS |      |            >f4 |  16.8e |
# |            MASS |      |            >f4 |  16.8e |
# |            XCEN |      |            >f8 | 25.17e |
# |            YCEN |      |            >f8 | 25.17e |
# |           NPROF |      |   ('>f4',(7,)) |  16.8e |
# |        PROFMEAN |      | ('>f4',(105,)) |  16.8e |
# |   PROFMEAN_IVAR |      | ('>f4',(105,)) |  16.8e |
# |         QSTOKES |      | ('>f4',(105,)) |  16.8e |
# |         USTOKES |      | ('>f4',(105,)) |  16.8e |
# |        BASTOKES |      | ('>f4',(105,)) |  16.8e |
# |       PHISTOKES |      | ('>f4',(105,)) |  16.8e |
# |       PETROFLUX |      |   ('>f4',(7,)) |  16.8e |
# |  PETROFLUX_IVAR |      |   ('>f4',(7,)) |  16.8e |
# |       FIBERFLUX |      |   ('>f4',(7,)) |  16.8e |
# |  FIBERFLUX_IVAR |      |   ('>f4',(7,)) |  16.8e |
# |            BA50 |      |            >f4 |  16.8e |
# |           PHI50 |      |            >f4 |  16.8e |
# |            BA90 |      |            >f4 |  16.8e |
# |           PHI90 |      |            >f4 |  16.8e |
# |      SERSICFLUX |      |   ('>f4',(7,)) |  16.8e |
# | SERSICFLUX_IVAR |      |   ('>f4',(7,)) |  16.8e |
# |        SERSIC_N |      |            >f4 |  16.8e |
# |       SERSIC_BA |      |            >f4 |  16.8e |
# |      SERSIC_PHI |      |            >f4 |  16.8e |
# |       ASYMMETRY |      |   ('>f4',(7,)) |  16.8e |
# |          CLUMPY |      |   ('>f4',(7,)) |  16.8e |
# |          DFLAGS |      |   ('>i4',(7,)) |    12i |
# |             AID |      |            >i4 |    12i |
# |             PID |      |            >i4 |    12i |
# |        DVERSION |      |            |S8 |  3689s |
# |       PROFTHETA |      |  ('>f4',(15,)) |  16.8e |
# |      PETROTHETA |      |            >f4 |  16.8e |
# |       PETROTH50 |      |            >f4 |  16.8e |
# |       PETROTH90 |      |            >f4 |  16.8e |
# |     SERSIC_TH50 |      |            >f4 |  16.8e |
# |           OBJNO |      |            >i4 |    12i |
# |           PLATE |      |            >i4 |    12i |
# |         FIBERID |      |            >i4 |    12i |
# |             MJD |      |            >i4 |    12i |
# |           COEFF |      |   ('>f4',(7,)) |  16.8e |
# |           VDISP |      |            >f4 |  16.8e |
# |           D4000 |      |            >f4 |  16.8e |
# |        D4000ERR |      |            >f4 |  16.8e |
# |              FA |      |            >f4 |  16.8e |
# |           FAERR |      |            >f4 |  16.8e |
# |          S2FLUX |      |            >f4 |  16.8e |
# |       S2FLUXERR |      |            >f4 |  16.8e |
# |            S2EW |      |            >f4 |  16.8e |
# |         S2EWERR |      |            >f4 |  16.8e |
# |         S2VMEAS |      |            >f4 |  16.8e |
# |         S2VMERR |      |            >f4 |  16.8e |
# |         S2RATIO |      |            >f4 |  16.8e |
# |          HAFLUX |      |            >f4 |  16.8e |
# |       HAFLUXERR |      |            >f4 |  16.8e |
# |            HAEW |      |            >f4 |  16.8e |
# |         HAEWERR |      |            >f4 |  16.8e |
# |         HAVMEAS |      |            >f4 |  16.8e |
# |         HAVMERR |      |            >f4 |  16.8e |
# |          N2FLUX |      |            >f4 |  16.8e |
# |       N2FLUXERR |      |            >f4 |  16.8e |
# |            N2EW |      |            >f4 |  16.8e |
# |         N2EWERR |      |            >f4 |  16.8e |
# |         N2VMEAS |      |            >f4 |  16.8e |
# |         N2VMERR |      |            >f4 |  16.8e |
# |          HBFLUX |      |            >f4 |  16.8e |
# |       HBFLUXERR |      |            >f4 |  16.8e |
# |            HBEW |      |            >f4 |  16.8e |
# |         HBEWERR |      |            >f4 |  16.8e |
# |         HBVMEAS |      |            >f4 |  16.8e |
# |         HBVMERR |      |            >f4 |  16.8e |
# |          O1FLUX |      |            >f4 |  16.8e |
# |       O1FLUXERR |      |            >f4 |  16.8e |
# |            O1EW |      |            >f4 |  16.8e |
# |         O1EWERR |      |            >f4 |  16.8e |
# |         O1VMEAS |      |            >f4 |  16.8e |
# |         O1VMERR |      |            >f4 |  16.8e |
# |          O2FLUX |      |            >f4 |  16.8e |
# |       O2FLUXERR |      |            >f4 |  16.8e |
# |            O2EW |      |            >f4 |  16.8e |
# |         O2EWERR |      |            >f4 |  16.8e |
# |         O2VMEAS |      |            >f4 |  16.8e |
# |         O2VMERR |      |            >f4 |  16.8e |
# |          O3FLUX |      |            >f4 |  16.8e |
# |       O3FLUXERR |      |            >f4 |  16.8e |
# |            O3EW |      |            >f4 |  16.8e |
# |         O3EWERR |      |            >f4 |  16.8e |
# |         O3VMEAS |      |            >f4 |  16.8e |
# |         O3VMERR |      |            >f4 |  16.8e |
# |           AHGEW |      |            >f4 |  16.8e |
# |        AHGEWERR |      |            >f4 |  16.8e |
# |           AHDEW |      |            >f4 |  16.8e |
# |        AHDEWERR |      |            >f4 |  16.8e |
# |           NE3EW |      |            >f4 |  16.8e |
# |        NE3EWERR |      |            >f4 |  16.8e |
# |           NE5EW |      |            >f4 |  16.8e |
# |        NE5EWERR |      |            >f4 |  16.8e |
# |              AV |      |            >f4 |  16.8e |
# |         S2NSAMP |      |            >f4 |  16.8e |
# |           RACAT |      |            >f8 | 25.17e |
# |          DECCAT |      |            >f8 | 25.17e |
# |       ZSDSSLINE |      |            >f4 |  16.8e |
# |          SURVEY |      |            |S6 |  3689s |
# |     PROGRAMNAME |      |           |S23 |  3689s |
# |    PLATEQUALITY |      |            |S8 |  3689s |
# |            TILE |      |            >i4 |    12i |
# |         PLUG_RA |      |            >f8 | 25.17e |
# |        PLUG_DEC |      |            >f8 | 25.17e |
# ----------------------------------------------------


import atpy, os
from astropy.io import fits
from pylab import *
from LCScommon import *

def cullNSA():
    infile=homedir+'research/NSA/nsa_v0_1_2.fits'
    if homedir.find('home') > -1:
        infile=homedir+'research/NSA/my_nsa_v0_1_2.fits.gz'
    ndat=fits.getdata(infile)
    positionflag=(ndat.RA > ramin) & (ndat.RA < ramax) & (ndat.DEC > decmin) & (ndat.DEC < decmax)
    redshiftflag=(ndat.Z > zmin) & (ndat.Z < zmax)
    keepflag=positionflag & redshiftflag
    outfile=homedir+'research/NSA/NSA_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    fits.writeto(outfile,ndat[keepflag])

    '''
    infile=homedir+'research/NSA/nsa_v0_1_2.fits'
    if homedir.find('home') > -1:
        infile=homedir+'research/NSA/my_nsa_v0_1_2.fits.gz'
    ndat=atpy.Table(infile)
    positionflag=(ndat.RA > ramin) & (ndat.RA < ramax) & (ndat.DEC > decmin) & (ndat.DEC < decmax)
    redshiftflag=(ndat.Z > zmin) & (ndat.Z < zmax)
    keepflag=positionflag & redshiftflag
    nsubset=ndat[keepflag]
    outfile=homedir+'research/NSA/NSA_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    nsubset.write(outfile)
    '''

    infile=homedir+'research/NSA/nsa_v0_1_2_wise.fits.gz'
    wdat=fits.getdata(infile)
    outfile=homedir+'research/NSA/NSA_WISE_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    fits.writeto(outfile,wdat[keepflag])

    infile=homedir+'research/NSA/nsa_v1_2_fsps_v2.4_miles_chab_charlot_sfhgrid01.fits'
    mdat=fits.getdata(infile)
    outfile=homedir+'research/NSA/NSA_JM_MSTAR_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    fits.writeto(outfile,mdat[keepflag])

def cullGalaxyZoo():
    infile=homedir+'research/GalaxyZoo/GalaxyZoo1_DR_table2.fits'
    zdat=atpy.Table(infile,type='fits')
    zRA=zeros(len(zdat.RA),'d')
    zDEC=zeros(len(zdat.RA),'d')
    for i in range(len(zdat.RA)):
        r=zdat.RA[i].split(':')
        zRA[i]=(float(r[0])+float(r[1])/60.+float(r[2])/3600.)*15
        d=zdat.RA[i].split(':')
        zDEC[i]=(float(d[0])+float(d[1])/60.+float(d[2])/3600.)
    positionflag=(zRA > ramin) & (zRA < ramax) & (zDEC > decmin) & (zDEC < decmax)
    keepflag=positionflag

    zsubset=zdat.where(keepflag)
    outfile=homedir+'research/GalaxyZoo/GalaxyZoo1_DR_table2_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    zsubset.write(outfile)

    infile=homedir+'research/GalaxyZoo/GalaxyZoo1_DR_table3.fits'
    zdat=atpy.Table(infile,type='fits')
    zRA=zeros(len(zdat.RA),'d')
    zDEC=zeros(len(zdat.RA),'d')
    for i in range(len(zdat.RA)):
        r=zdat.RA[i].split(':')
        zRA[i]=(float(r[0])+float(r[1])/60.+float(r[2])/3600.)*15
        d=zdat.RA[i].split(':')
        zDEC[i]=(float(d[0])+float(d[1])/60.+float(d[2])/3600.)
    positionflag=(zRA > ramin) & (zRA < ramax) & (zDEC > decmin) & (zDEC < decmax)
    keepflag=positionflag

    zsubset=zdat.where(keepflag)
    outfile=homedir+'research/GalaxyZoo/GalaxyZoo1_DR_table3_LCSregion.fits'
    if os.path.exists(outfile):
        os.remove(outfile)            
    zsubset.write(outfile)

if __name__ == '__main__':
    cullNSA()
#cullGalaxyZoo()
