#!/usr/bin/env python

# Goal:  
# * match Simard GSW Legacy Catalog to NSA
# 
# Required Files:
# * The NASA Sloan Atlas catalog (nsa_v0_1_2.fits) found here http://www.nsatlas.org/data
# * lsp_gswlc_galex-sdss-wise_multi_x1_multi_v1_cat.fits (https://archive.stsci.edu/prepds/gswlc/)
# 


import numpy as np
from astropy.io import fits
import argparse
from astropy.coordinates import ICRS, SkyCoord
from astropy import units as u
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description ='Match the NSA catalog with the Simard catalogs')
parser.add_argument('--gsw-path', dest = 'gswpath', default = '/Users/rfinn/Dropbox/Research/GSWLC/', action = 'store_true', help = 'path to Salim+2016 GXWLC catalog ')
parser.add_argument('--nsa-path', dest = 'path', default = '/Users/rfinn/research/NSA/', action = 'store_true', help = 'path to NSA catalog nsa_v0_1_2.fits')
args = parser.parse_args()
gswpath = args.gswpath
nsapath = args.path

gswfile ='hlsp_gswlc_galex-sdss-wise_multi_x1_multi_v1_cat.fits'
gswdat = fits.getdata(gswpath+gswfile)
nsadat = fits.getdata(nsapath+'nsa_v0_1_2.fits')


nsacat = SkyCoord(nsadat.RA*u.degree,nsadat.DEC*u.degree,frame='icrs')
gswcat = SkyCoord(gswdat.RA*u.degree,gswdat.DECL*u.degree,frame='icrs')


# match Simard+2011 Table 1 to NSA
index,dist2d,dist3d = nsacat.match_to_catalog_sky(gswcat)


# only keep matches with matched RA and Dec w/in 1 arcsec
matchflag = dist2d.degree < 4./3600


# write out line-matched catalog
outfile=nsapath+'NSA-GSWLC.fits'
matchedarray1=np.zeros(len(nsadat),dtype=gswdat.dtype)
matchedarray1[matchflag] = gswdat[index[matchflag]]
fits.writeto(outfile,matchedarray1,overwrite=True)


