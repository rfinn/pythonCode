#!/usr/bin/env python

import os

prefix=['cl1040','cl105411','cl105412','cl1216','cl1354','cl1037','cl1232','cl1227']
directory=['r13810432','r13810176','r13812480','r13810688','r13810944','r17053184','r17053696','r17053952']
for i in range(len(prefix)):
    s='cp /Users/rfinn/clusters/spitzer/'+str(prefix[i])+'/mips/24/'+str(directory[i])+'/ch1/bcd/pbcd/Combine/mosaic.fits /Users/rfinn/clusters/spitzer/final-24-images/'+str(prefix[i])+'mosaic.fits'
    os.system(s)
