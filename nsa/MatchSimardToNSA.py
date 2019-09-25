
# coding: utf-8

# Goal:  
# Combine the NSAYang catalog with the Simard catalog to include sersic fits
# 
# Required Files:
# * The NASA Sloan Atlas catalog (nsa_v0_1_2.fits) found here http://www.nsatlas.org/data
# * asu.fit
# 
# Obtaining asu.fit:
# 1. http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJS/196/11
# 2. Check all boxes and click "Join selected tables"
# 3. Scroll down and check "ALL col" then uncheck "All", "Sloan", and "DR7"
# 4. Click any of the submit buttons
# 5. On the left in the "Preferences" box change "max" to unlimited, "HTML Table" to "FITS (binary) table, and then click submit
# 
# Notes:
# * All entries in asu.fit are strings. Can be changed to floats using np.astype(np.float32)
# * updated on 10/10/16 by R Finn to streamline matching and preserve data types in matched output catalog.  Changed Grant's instructions to download asu.fits as a binary fits, not an ascii fits table.

# In[1]:

import csv
import numpy as np
from astropy.io import fits
import fnmatch
import time
import argparse
from astropy.coordinates import ICRS, SkyCoord
from astropy import units as u
from matplotlib import pyplot as plt
#get_ipython().magic(u'matplotlib inline')


try:
    parser = argparse.ArgumentParser(description ='Match the NSA catalog with the Simard catalogs')
    parser.add_argument('--simard-path', dest = 'spath', default = '/Users/rfinn/research/SimardSDSS2011/', action = 'store_true', help = 'path to Simard+2011 catalog asu.fit')
    parser.add_argument('--nsa-path', dest = 'path', default = '/Users/rfinn/research/NSA/', action = 'store_true', help = 'path to NSA catalog nsa_v0_1_2.fits')
    args = parser.parse_args()
    simardpath = args.spath
    nsapath = args.path
except:
    simardpath = '/Users/rfinn/research/SimardSDSS2011/'
    nsapath = '/Users/rfinn/research/NSA/'


asu1 = fits.getdata(simardpath+'asu.fit',1)
asu2 = fits.getdata(simardpath+'asu.fit',2)
asu3 = fits.getdata(simardpath+'asu.fit',3)

nsadat = fits.getdata(nsapath+'nsa_v0_1_2.fits')


nsacat = SkyCoord(nsadat.RA*u.degree,nsadat.DEC*u.degree,frame='icrs')
simardcat = SkyCoord(asu1._RA*u.degree,asu1._DE*u.degree,frame='icrs')


# match Simard+2011 Table 1 to NSA
index,dist2d,dist3d = nsacat.match_to_catalog_sky(simardcat)


# only keep matches with matched RA and Dec w/in 1 arcsec
matchflag = dist2d.degree < 2./3600


# write out line-matched catalog
outfile='Simard1ToNSA.fits'
matchedarray1=np.zeros(len(nsadat),dtype=asu1.dtype)
matchedarray1[matchflag] = asu1[index[matchflag]]
fits.writeto(outfile,matchedarray1,clobber=True)


# write out line-matched catalog
outfile='Simard2ToNSA.fits'
matchedarray1=np.zeros(len(nsadat),dtype=asu2.dtype)
matchedarray1[matchflag] = asu2[index[matchflag]]
fits.writeto(outfile,matchedarray1,clobber=True)


# In[19]:

# write out line-matched catalog
outfile='Simard3ToNSA.fits'
matchedarray1=np.zeros(len(nsadat),dtype=asu3.dtype)
matchedarray1[matchflag] = asu3[index[matchflag]]
fits.writeto(outfile,matchedarray1,clobber=True)

