#!/usr/bin/env python
from pylab import *
import numpy as np
lirmin=10.95
lirminloz=10.75
xmin=1.e8
xmax=5.e12
ymin=2.e9
ymax=3.e12
figure()
files=['/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.Ediscs.Hiz.blue.dat','/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.Ediscs.Hiz.red.dat','/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.Ediscs.Liz.blue.dat', '/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.Ediscs.Liz.red.dat']
for i in range(1,5):
    subplot(2,2,i)
    infile=open(files[i-1],'r')
    t=np.loadtxt(infile)#,dtype='double')
    infile.close()
    x=(t[:,0])
    y=(t[:,1])
    plot(x,y,'k.')
    ax=gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

gemsfiles=['/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.GemsHiz.dat','/home/rfinn/research/clusters/spitzer/MasterTables/ForPascale.GemsLoz.dat']
for i in range(len(gemsfiles)):
    infile=open(gemsfiles[i],'r')
    t=np.loadtxt(infile)#,dtype='double')
    infile.close()
    x=(t[:,0])
    y=(t[:,1])
    fl=t[:,2]
    for j in range(2):
        if i < 1:
            nplot=i+j+1

        else:
            nplot=i+j+2
        subplot(2,2,nplot)
        print 'gems subplot = ',nplot,i,j
        flag=abs(fl-j)<.1
        plot(x[flag],y[flag],'.',color='w',markeredgecolor='k')
        
        axis([xmin,xmax,ymin,ymax])
show()
print 'hey'
