#!/usr/bin/env python
from pylab import *
#from ediscsReadSpitzerMaster import baseCluster
from spitzer1 import ediscs24

class cluster(ediscs24):
    def plotJK(self):
        plotflag=(self.magJ > -90) & (self.magK > -90)
        figure()
        plot()
