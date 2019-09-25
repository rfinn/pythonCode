#!/usr/bin/env python

'''
converting bootstrap.f to python

'''

#from pylab import *
import numpy as np
from scipy.interpolate import interp1d

def bootstrap(xdat,ndraw=1000):
    xdata=np.array(xdat,'d')
    ndat=len(xdata)
    actual_median=np.median(xdata)
    actual_mean=np.mean(xdata)
    actual_std=np.std(xdata)
    bmean=np.zeros(ndraw,'d')
    bmedian=np.zeros(ndraw,'d')
    # calculate mean and median for ndraw random realizations of data
    for i in range(ndraw):
        xboot=xdata[np.random.random_integers(0,ndat-1,ndat)]
        bmean[i]=np.mean(xboot)
        bmedian[i]=np.median(xboot)

    std_bmean=np.std(bmean)
    std_bmedian=np.std(bmedian)

    # sort array of bootstrap mean/median
    bmean.sort()
    bmedian.sort()
    bindex=np.arange(ndraw)

    # define intervals
    upper_intervals= np.array([.6827,0.9545,0.9973,  0.9999, 0.999994],'f')

  #'S' parameter for these limits
    #S_d = np.array([1.0, 1.282, 1.645, 1.960, 2.00, 2.326, 2.576, 3.000, 3.090, 3.291],'f')

    lower_intervals=1-(upper_intervals)

    func1=interp1d(bindex,bmean)
    upper_range_mean=func1(upper_intervals*(ndraw-1))
    lower_range_mean=func1(lower_intervals*(ndraw-1))


    func1=interp1d(bindex,bmedian)
    upper_range_median=func1(upper_intervals*(ndraw-1))
    lower_range_median=func1(lower_intervals*(ndraw-1))

    return np.mean(bmean),np.median(bmedian),lower_range_mean, upper_range_mean, lower_range_median, upper_range_median 
