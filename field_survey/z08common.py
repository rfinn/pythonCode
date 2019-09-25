#!/usr/bin/env python
import os
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

zmin=0.79
zmax=0.815
clusternames=['RDCSJ1317', 'WARP1350','RXJ1716','RXJ1821' ]
clusterRAh={'RDCSJ1317':13, 'WARP1350':13,'RXJ1716':17,'RXJ1821':18 }
clusterRAm={'RDCSJ1317':17, 'WARP1350':50,'RXJ1716':16,'RXJ1821':21 }
clusterRAs={'RDCSJ1317':21.7, 'WARP1350':48.5,'RXJ1716':49.6,'RXJ1821':32.9 }
clusterDech={'RDCSJ1317':29, 'WARP1350':60,'RXJ1716':67,'RXJ1821':68 }
clusterDecm={'RDCSJ1317':13, 'WARP1350':13,'RXJ1716':17,'RXJ1821':18 }
clusterDecs={'RDCSJ1317':13, 'WARP1350':13,'RXJ1716':17,'RXJ1821':18 }
clusterRAdeg={'RDCSJ1317':(13.+17/60.+21.7/3600)*15,'WARP1350':(13.+50./60.+48.5/3600.)*15 ,'RXJ1716':(17.+16./60+49.6/3600)*15,'RXJ1821': (18.+21./60+32.9/3600)*15}
clusterDecdeg={'RDCSJ1317':+29.+11./60+18./3600, 'WARP1350':+60.+7./60+7./3600,'RXJ1716':+67+8./60+30./3600,'RXJ1821':68.+27./60.+55./3600}    

clustersigma={'RDCSJ1317':1000., 'WARP1350':1000.,'RXJ1716':1000.,'RXJ1821':1000.}

clusterz={'RDCSJ1317':.804, 'WARP1350': .804,'RXJ1716': .81,'RXJ1821': .81}

newfirm_offset={'RDCSJ1317':[90.,360.], 'WARP1350': [210.,-150.],'RXJ1716':[-150.,-150.],'RXJ1821': [-150.,-150.]}
