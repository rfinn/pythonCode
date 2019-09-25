#!/usr/bin/env python
from pylab import *
RA=array([197.539174,235.980204,215.158852,247.140420,198.955991,231.599167,237.439958],'f')
Dec=array([29.295241,29.463944,43.725871, 29.321022,43.575129, 43.004556,43.057064],'f')
kissid=['K0225','K0847','K1759','K1038','K1516','K1953','K2042']

for i in  range(len(kissid)):
    s='/Users/rfinn/kiss/RegionsFiles/'+kissid[i]+'.reg'
    outfile3=open(s,'w')
    s='global color=cyan font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source \n'
    outfile3.write(s)
    s='fk5 \n'
    outfile3.write(s)
    s='circle(%12.8f, %12.8f, 10") \n'%(RA[i],Dec[i])
    outfile3.write(s)
    outfile3.close()
