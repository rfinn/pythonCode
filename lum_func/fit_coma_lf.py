#!/usr/bin/env python

import atpy
import numpy as np
import fit_lf


def comatest(logL,irweight):

    # Bai et al 2009 LF parameters
    alpha = -1.41
    logLstar = 10.5
    phistar = 6.5
    t=fit_lf.fit_schechter(logL,1./irweight,phistar,logLstar,alpha,nbin=8,xmin=7.5,xmax=10.6)
    return t
    

infile='/Users/rfinn/research/LocalClusters/NSAmastertables/CharyElbazTables/Coma_ce_lir.fits'
irdat=atpy.Table(infile)
infile='/Users/rfinn/research/LocalClusters/NSAmastertables/NSAwithAGC/Coma_NSAmastertable_topcat.fits'
ndat=atpy.Table(infile)
keepflag = irdat.MATCHFLAG24 & ndat.MEMBFLAG
logL=np.log10(irdat.LIR[keepflag])
irweight=irdat.MIPS_WEIGHT[keepflag]

    
