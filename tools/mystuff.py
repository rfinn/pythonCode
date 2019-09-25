import sys, os
#import Numeric as N
import numpy as N
import scipy
import weightedstats as ws
#import scipy.stats
import pylab
#from Scientific.Functions.Romberg import romberg
from scipy.integrate import romberg
from math import *
import matplotlib.mlab
#import ppgplot
import binomial
import poisson
OmegaL = 0.7
OmegaM = 0.3
cl=3.e5

def clusterMass(sigma,z,h):
    mcl=1.2e15*(sigma/1000.)**2*(1./sqrt(OmegaL+OmegaM*(1+z)**3))*1./h
    return mcl

                                 
def drawbox(data,style):#feed in center x,y,dx,dy,rotation E of N
    #xcoords of unrotated box, going around CCW
    xl=pylab.array([data[0]-0.5*data[2],data[0]+0.5*data[2],data[0]+0.5*data[2],data[0]-0.5*data[2],data[0]-0.5*data[2]],'d')
    yl=pylab.array([data[1]-0.5*data[3],data[1]-0.5*data[3],data[1]+0.5*data[3],data[1]+0.5*data[3],data[1]-0.5*data[3] ],'d')

    xl=pylab.array([-0.5*data[2],+0.5*data[2],+0.5*data[2],-0.5*data[2],-0.5*data[2]],'d')
    yl=pylab.array([-0.5*data[3],-0.5*data[3],+0.5*data[3],+0.5*data[3],-0.5*data[3] ],'d')

    ang=data[4]*pi/180.*-1.#convert rotation to radians
    #rotate coordinates
    xp=pylab.cos(ang)*xl-pylab.sin(ang)*yl
    yp=pylab.sin(ang)*xl+pylab.cos(ang)*yl

    #put back on absolute scale
    xp=data[0]+xp
    yp=data[1]+yp
    #draw rotated box
    pylab.plot(xp,yp,style)


def ellipse(ra,rb,ang,x0,y0,Nb=50):
    '''ra - major axis length
    rb - minor axis length
    ang - angle
    x0,y0 - position of centre of ellipse
    Nb - No. of points that make an ellipse
    
    based on matlab code ellipse.m written by D.G. Long,
    Brigham Young University, based on the
    CIRCLES.m original
    written by Peter Blattner, Institute of Microtechnology,
    University of
    Neuchatel, Switzerland, blattner@imt.unine.ch
    '''

    xpos,ypos=x0,y0
    radm,radn=ra,rb
    an=ang

    co,si=cos(an),sin(an)
    the=linspace(0,2*pi,Nb)
    X=radm*cos(the)*co-si*radn*sin(the)+xpos
    Y=radm*cos(the)*si+co*radn*sin(the)+ypos
    return X,Y

def getvfromz(z):
    k=(1+z)**2
    v=3.e5*(k-1)/(1+k)#velocity in km/s
    return v

def getzfromv(v):
    z=N.sqrt((1+v/3.e5)/(1-v/3.e5))-1
    return z

def beginsWith(pattern,string):
    if string[0:len(pattern)] == pattern:
	return 1
    else:
	return 0
def binitbinsfix(binsleft,binsright,x,y):#give bins for binning
    nbin=len(binsleft)
    xbin=N.zeros(len(binsleft),'f')
    ybin=N.zeros(len(binsleft),'f')
    ybinerr=N.zeros(len(binsleft),'f')
    y=N.take(y,N.argsort(x))#sort y according to x rankings
    x=N.take(x,N.argsort(x))

    j=-1
    for i in range(len(xbin)):
	xmin=binsleft[i]
	xmax=binsright[i]
	yb=[]
	for j in range(len(x)):
	    if x[j] > xmin:
		yb.append(y[j])
	    if x[j] > xmax:
		yb=N.array(yb,'f')
		xbin[i]=0.5*(xmin+xmax)
		try:
		    ybin[i]=pylab.mean(yb)
		    ybinerr[i]=pylab.std(yb)/N.sqrt(1.*len(yb))
		except ZeroDivisionError:
		    print "warning: ZeroDivision error in binitbinsfix"
		    ybin[i]=0.
		    ybinerr[i]=0.

		break

    return xbin,ybin,ybinerr
def scipyhist2(bins,y):#give bins for binning
    nbin=len(bins)
    xbin=N.zeros(len(bins),'f')
    ybin=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    y=N.take(y,N.argsort(y))
    xmin=min(bins)
    xmax=max(bins)
    xbinnumb=((y-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    for i in range(len(xbin)-1):
	xmin=bins[i]
	xmax=bins[i+1]
	yb=0.
	for j in range(len(y)):
	    if y[j] > xmin:
		yb= yb + 1.
	    if y[j] > xmax:
		ybin[i]=yb
		ybinerr[i]=N.sqrt(yb)
		xbin[i]=0.5*(xmin+xmax)

		break

    return xbin,ybin,ybinerr

def completeness(bins,y,yflag):#y is an array of input fluxes, yflag tells if sources was recovered
    #for i in range(len(y)):
#	print y[i],yflag[i]
    nbin=len(bins)
    xbin=N.zeros(len(bins),'f')
    ybin=N.zeros(len(bins),'f')
    ybin2=N.zeros(len(bins),'f')
    yratio=N.zeros(len(bins),'f')
    ybinerr=N.zeros(len(bins),'f')
    yflag=N.take(yflag,N.argsort(y))
    y=N.take(y,N.argsort(y))
    #for i in range(20):
	#print 'sorted',i,y[i],yflag[i]
    xmin=min(bins)
    xmax=max(bins)
    xbinnumb=((y-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    for i in range(len(xbin)-1):
	xmin=bins[i]
	xmax=bins[i+1]
	yb=0.
	yf=0.
	for j in range(len(y)):
	    if y[j] > xmin:
		yb= yb + 1.
		if yflag[j] > 0.1:
		    yf=yf + 1.
	    if y[j] > xmax:
		ybin[i]=yb
		ybin2[i]=yf
		(yratio[i],ybinerr[i],t)=ratioerror(1.*yf,1.*yb)
		xbin[i]=0.5*(xmin+xmax)
		#print 'completeness',xbin[i],yf,yb,yratio[i],ybinerr[i]
		break


    return xbin,yratio,ybinerr

def verticalhist(y,ymin,ymax,nbin):#give bins for binning
    #nbin=len(bins)
    xbin=N.zeros(nbin,'f')
    ybin=N.zeros(nbin,'f')
    ybinerr=N.zeros(nbin,'f')
    dbin=(ymax-ymin)/(1.*nbin)
    bins=N.arange(ymin,(ymax+dbin),dbin)
    y=N.take(y,N.argsort(y))
    xmin=min(bins)
    xmax=max(bins)
    xbinnumb=((y-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    for i in range(len(xbin)-1):
	xmin=bins[i]
	xmax=bins[i+1]
	yb=0.
	for j in range(len(y)):
	    if y[j] > xmin:
		yb= yb + 1.
	    if y[j] > xmax:
		ybin[i]=yb
		ybinerr[i]=N.sqrt(yb)
		#xbin[i]=0.5*(xmin+xmax)
		xbin[i]=xmax
		break
    drawhist(ybin,xbin)

def horizontalhist(x,xmin,xmax,nbin):#give bins for binning
    #nbin=len(bins)
    xbin=N.zeros(nbin,'f')
    ybin=N.zeros(nbin,'f')
    ybinerr=N.zeros(nbin,'f')
    dbin=(xmax-xmin)/(1.*nbin)
    bins=N.arange(xmin,(xmax+dbin),dbin)
    x=N.take(x,N.argsort(x))
    xmin=min(bins)
    xmax=max(bins)
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    #print "from within horizontal hist"
    #print "bins = ",bins
    for i in range(len(xbin)-1):
	xmin=bins[i]
	xmax=bins[i+1]
	xbin[i]=xmin
	yb=0.
	for j in range(len(x)):
	    if (x[j] > xmin) & (x[j] <= xmax):
		yb= yb + 1.
	ybin[i]=yb
	ybinerr[i]=N.sqrt(yb)
	#print i,len(xbin),xbin[i],xmin,xmax,ybin[i],bins[i]
    xbin[(len(xbin)-1)]=xbin[(len(xbin)-2)]+dbin
    #xbin=bins
    #print "from w/in horizontal hist, ybin = ",ybin
    return xbin,ybin,ybinerr



def binitbinsmabs(xmin,xmax,nbin,x,y):#use equally spaced bins
    y=10.**(y)
    dx=float((xmax-xmin)/(nbin))
    xbin=N.arange(xmin,(xmax),dx)+dx/2.
    ybin=N.zeros(len(xbin),'d')
    ybinerr=N.zeros(len(xbin),'d')
    xbinnumb=N.array(len(x),'d')
    x1=N.compress((x >= xmin) & (x <= xmax),x)
    y1=N.compress((x >= xmin) & (x <= xmax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    j=-1
    for i in range(len(xbin)):
	ydata=N.compress(abs(xbinnumb-float(i))<.5,y)
	try:
	    ybin[i]=N.log10(N.average(ydata))
	    #ybin[i]=pylab.median(ydata)
	    ybinerr[i]=pylab.std(ydata)/N.sqrt(float(len(ydata)))
	except ZeroDivisionError:
	    ybin[i]=-99.
	    ybinerr[i]=-99.

    return xbin,ybin,ybinerr

def binitbins(xmin,xmax,nbin,x,y):#use equally spaced bins
    dx=float((xmax-xmin)/(nbin))
    xbin=N.arange(xmin,(xmax),dx)+dx/2.
    ybin=N.zeros(len(xbin),'d')
    ybinerr=N.zeros(len(xbin),'d')
    xbinnumb=N.array(len(x),'d')
    x1=N.compress((x >= xmin) & (x <= xmax),x)
    y1=N.compress((x >= xmin) & (x <= xmax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    j=-1
    for i in range(len(xbin)):
	ydata=N.compress(abs(xbinnumb-float(i))<.5,y)
	try:
	    #ybin[i]=N.average(ydata)
	    ybin[i]=pylab.median(ydata)
	    ybinerr[i]=pylab.std(ydata)/N.sqrt(float(len(ydata)))
	except ZeroDivisionError:
	    ybin[i]=0.
	    ybinerr[i]=0.
    return xbin,ybin,ybinerr

def binitbins_weighted(xmin,xmax,nbin,x,y,yerr,medianflag=False):#use equally spaced bins
    bins = N.linspace(xmin,xmax,nbin)
    dx = bins[1]-bins[0]
    idx = N.digitize(x,bins)-1
    print idx
    ybin=N.zeros(len(bins),'d')
    ybinerr=N.zeros(len(bins),'d')

    j=-1
    for i in range(nbin):
        try:
            #print 'xrange = ',bins[i],dx/2.
            #print 'x values = ',x[idx == i]
            #print y[idx == i]
            #print 1./yerr[idx == i]
            w=1./yerr[i==idx]
            weight = w/(N.sum(w)/(1.*len(w)))
            if medianflag:
                ybin[i]=ws.weighted_median(y[idx == i].tolist(),weights=weight.tolist())
            else:
                ybin[i]=N.average(y[idx == i],weights=weight)
            print 'weighted median = ',ybin[i]
            #ybinerr[i]=N.sqrt(N.sum(yerr[idx == i]**2))/(1.*(float(len(y[idx == i]))))
            ybinerr[i] = N.sqrt(N.sum(yerr[i==idx]**2))/(1.*len(w))
        except ZeroDivisionError:
            ybin[i]=0.
            ybinerr[i]=0.
    return bins+dx/2.,ybin,ybinerr

def binitbinsmedian(xmin,xmax,nbin,x,y):#use equally spaced bins
    dx=float((xmax-xmin)/(nbin))
    xbin=N.arange(xmin,(xmax),dx)+dx/2.
    #print "within binitbins"
    #print "xbin = ",xbin
    #print "dx = ",dx
    #print "xmax = ",xmax
    #print "xmin = ",xmin
    ybin=N.zeros(len(xbin),'d')
    ybinerr=N.zeros(len(xbin),'d')
    xbinnumb=N.array(len(x),'d')
    x1=N.compress((x >= xmin) & (x <= xmax),x)
    y1=N.compress((x >= xmin) & (x <= xmax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    j=-1
    for i in range(len(xbin)):
        ydata=y[abs(xbinnumb-float(i))<.5]
        xdata=x[abs(xbinnumb-float(i))<.5]
        try:
            #ybin[i]=N.average(ydata)
            ybin[i]=pylab.median(ydata)
            xbin[i]=pylab.median(xdata)
            ybinerr[i]=pylab.std(ydata)/N.sqrt(float(len(ydata)))
        except ZeroDivisionError:
            ybin[i]=0.
            xbin[i]=0
            ybinerr[i]=0.
    return xbin,ybin,ybinerr
def vofz(z):#recession velocity corresponding to redshift z
    v = cl*((1+z)**2-1)/((1+z)**2+1)
    return v #velocity in km/s
def relsumv(v1,v2):#relativistic sum of velocities
    v = (v1 + v2)/(1+v1*v2/cl**2)
    return v #velocity in km/s
def reldiffv(v1,v2):#relativistic sum of velocities
    v = (v1 - v2)/(1-v1*v2/cl**2)
    return v #velocity in km/s

def zofv(v):#get redshift corresponding to recession vel v
    z=N.sqrt((1+v/cl)/(1-v/cl))-1.
    return z

def psplotinit(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(8.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(6)  #line width

def psplotinitlarge(output):
    file=output+"/vcps"
    ppgplot.pgbeg(file,1,1)
    ppgplot.pgpap(10.,1.)
    ppgplot.pgsch(1.7) #font size
    ppgplot.pgslw(6)  #line width

def putlabelpt(x,y,dx,dy,label,symb):
    ppgplot.pgtext(x,y,label)
    xlin = N.array([x-dx],'f')
    ylin = N.array([y+dy],'f')
    ppgplot.pgpt(xlin,ylin,symb)

def putlabelptc(x,y,dx,dy,label,symb,colors):
    defsci = ppgplot.pgqci()
    ppgplot.pgtext(x,y,label)
    xlin = N.array([x-dx],'f')
    ylin = N.array([y+dy],'f')
    ppgplot.pgsci(colors)
    ppgplot.pgpt(xlin,ylin,symb)
    ppgplot.pgsci(defsci)

def pgpnts(x,y,symbols):
    for i in range(len(x)):
        #print i,x[i],y[i],symbols[i]
        xt = N.array([x[i]],'f')
        yt = N.array([y[i]],'f')
        symb = symbols[i]
        ppgplot.pgpt(xt,yt,symb)
def pgpntsc(x,y,symbols):
    #print "pgpnts",len(x),len(y),len(symbols)
    for i in range(len(x)):
        #print i,x[i],y[i],symbols[i]
        xt = N.array([x[i]],'f')
        yt = N.array([y[i]],'f')
        symb = symbols[i]
        ppgplot.pgpt(xt,yt,symb)

def dospearold(x,y):
    #(a,b)=scipy.stats.spearmanr(x,y)
    #print "rank correl = %6.5f %10.9f" % (a,b)
    #(a,b)=scipy.stats.kendalltau(x,y)
    #print "kendall tau = %6.3f %6.5f" % (a,b)
    x=x

def dospear(x,y):
    print "Spearman rank test"
    outfile=open('speardata','w')
    for i in range(len(x)):
        outfile.write("%8.3f %8.3f \n" % (x[i],y[i]))
    outfile.close()
    os.system('spearrank')
    infile=open('spear.out','r')
    for line in infile:
        (dspear, zd, probd, rs, probrs)=line.split()
    #print "sum-squared difference of ranks"
    #print "d_spear = ",dspear
    print "number of standard deviations by which d deviates from null hypothesis"
    print "zd = ", zd
    print "two-sided significance level of this deviation probd = ", probd
    print "Spearman's rank correlation rs = ", rs
    print "two-sided significance of its dev from zero probrs = ", probrs
    print "small value of probd or probs indicates sig"
    print "of correlation (rs > 0) or anticorrelation (rs < 0)"
    print "***************************************************"
    print "***************************************************"
    
#def dolinregress(x,y):
#    (a,b,c,d,e)=scipy.stats.linregress(x,y)#
#
#    print "linear regression slope,inter,r,prob,std of est = %6.5f %8.7f %8.7f %8.7f %8.7f" % (a,b,c,d,e)
#    #(a,b)=pylab.polyfit(x,y,1)
#    print "pylab polyfit = %6.3f %6.5f" % (a,b)
#    #(a,b)=scipy.stats.kendalltau(x,y)
#    #print "kendall tau = %6.3f %6.5f" % (a,b)
#    #x=x
#    return a,b
def xysimple(x,y,xlabel,ylabel):
    xmax=max(x)
    xmin=min(x)
    ymax=max(y)
    ymin=min(y)
    ppgplot.pgbox("",0.0,0,"",0.0,0)
    ppgplot.pgenv(xmin,xmax,ymin,ymax,0)
    ppgplot.pglab(xlabel,ylabel,"")
    ppgplot.pgpt(x,y,3)

def binit(x,y,n):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    y=N.take(y,N.argsort(x))#sort y according to x rankings
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        xbin[i]=pylab.median(x[nmin:nmax])
        ybin[i]=pylab.median(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #print y[nmin:nmax]
        ybinerr[i]=pylab.std(y[nmin:nmax])/sqrt(1.*len(y[nmin:nmax]))
    return xbin, ybin, ybinerr
def binitave(x,y,n):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    y=N.take(y,N.argsort(x))#sort y according to x rankings
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        #xbin[i]=pylab.median(x[nmin:nmax])
        #ybin[i]=pylab.median(y[nmin:nmax])
        xbin[i]=N.average(x[nmin:nmax])
        ybin[i]=N.average(y[nmin:nmax])
        ybinerr[i]=pylab.std(y[nmin:nmax])/sqrt(1.*len(y[nmin:nmax]))
    return xbin, ybin, ybinerr

def biniterr(x,y,n):#bin arrays x, y into n equally-populated bins, returning xbin,ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        #xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        #ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        xbin[i]=N.average(x[nmin:nmax])
        ybin[i]=N.average(y[nmin:nmax])
        ybinerr[i]=pylab.std(y[nmin:nmax])/N.sqrt(1.*(nmax-nmin))
    return xbin, ybin, ybinerr

def getbins(x,y,n):#get left side of equally populated bins 
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbinright=N.zeros(n,'f')
    xbinleft=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        #xbin[i]=scipy.stats.stats.median(x[nmin:nmax])
        #ybin[i]=scipy.stats.stats.median(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        xbinleft[i]=(x[nmin])
        xbinright[i]=(x[nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=pylab.std(y[nmin:nmax])#/N.sqrt(1.*(nmax-nmin))
    return xbinleft,xbinright

def sortwindex(x):#sort array, return sorted array and array containing indices for unsorted array
    nx=len(x)
    y=N.arange(0,nx,1)#index of array x
    y=N.take(y,N.argsort(x))#sort y according to x rankings
    xs=N.take(x,N.argsort(x))
    return xs, y#xs=sorted array, y=indices in unsorted array

def findmatch(x,ysort,delta):#find match to x in sorted array ysort, where N.sqrt((x-ysort[i])**2) < tolerance
    #return list matches containing indices of ysort that meet tolerance
    ny=len(ysort)
    ytest=int(round(ny/2.))#index to test
    dy=int(round(ytest/2.))
    d=10000.
    none=0
    flag=1
    #print ny,ytest,d
    while d > delta:
        #print ytest,dy,x,len(ysort),ysort[ytest],d,delta
        d=dist(x,ysort[ytest])
        
        if d > delta:
            #print "d > delta", ysort[ytest],x
            if ysort[ytest] > x:
                #print "ysort > x"
                ytest=ytest-dy
                dy=int(round(dy/2.))
                if ytest < 0:
                    ytest = 0
            else:
                #print "ysort < x"
                ytest=ytest+dy
                dy=int(round(dy/2.))
                if ytest > (ny-1):
                    ytest=ny-1
        if dy == 1:
            none=none+1
        if (dy < 1) | (none > 3):
            # print "no match!!!!"
            match=-999
            flag=0
            break
    else:
        imatch=ytest#matched index
        #find min index that satisfies tolerance
        imin=imatch-1
        d=dist(x,ysort[imin])
        while (d < delta):#
            imin=imin-1
            if (imin < 0):
                break
            d=dist(x,ysort[imin])
        imin=imin+1
        #find max index that satisfies tolerance
        imax=imatch+1
        if imax > (ny - 1):
            imax= (ny-1)
        else:
            d=dist(x,ysort[imax])
            while (d < delta):#
                imax=imax+1
                if (imax > (ny-1)):
                    break
                d=dist(x,ysort[imax])
            imax=imax-1
        if imax > imin:
            match=N.arange(imin,(imax+1),1)#imin to imax, inclusive
        if imax == imin:
            match = [imatch]
        if imax < imin:
            print "Error in findmatch, imax < imin"
    return match,flag#indices of ysort that meet delta

def dist(x,y):
    diff=N.sqrt((x-y)**2)
    return diff

def binitmin(x,y,n):#bin arrays x, y into n bins, returning xbin, and min value of y points, ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    #ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        xbin[i]=pylab.median(x[nmin:nmax])
        ybin[i]=min(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=scipy.stats.std(y[nmin:nmax])
    return xbin, ybin#, ybinerr

def binitsumequal(x,y,n):#bin arrays x, y into n bins, returning xbin,ybin
    nx=len(x)
    #x=N.array(x,'f')
    #y=N.array(y,'f')
    y=N.take(y,N.argsort(x))
    x=N.take(x,N.argsort(x))
    xbin=N.zeros(n,'f')
    ybin=N.zeros(n,'f')
    #ybinerr=N.zeros(n,'f')
    for i in range(n):
        nmin=i*int(float(nx)/float(n))
        nmax=(i+1)*int(float(nx)/float(n))
        xbin[i]=N.average(x[nmin:nmax])
        ybin[i]=N.sum(y[nmin:nmax])
        #xbin[i]=N.average(x[nmin:nmax])
        #ybin[i]=N.average(y[nmin:nmax])
        #ybinerr[i]=scipy.stats.std(y[nmin:nmax])
    return xbin, ybin#, ybinerr


def drawbinned(x,y,nbin):
    xbin,ybin=binit(x,y,nbin)
    #ppgplot.pgsci(2)
    ppgplot.pgline(xbin,ybin)

def drawhist(x,y): #draw histogram of binned data, x=left side of bins, y=number per bin
    #print "from within drawhist, x,y  = ",x,y
    x1=[]
    y1=[]
    n=len(x)
    dx=x[1]-x[0]
    for i in range(n):
        x1.append(float(x[i]))
        x1.append(float(x[i]))
    x1.append(x[(n-1)]+dx)
    x1.append(x[(n-1)]+dx)
    #x1.append(max(x)+dx)
    #x1.append(max(x)+dx)
    y1.append(0.)
    for i in range(n):
        y1.append(float(y[i]))
        y1.append(float(y[i]))
    y1.append(0.)
    #print "Within drawhist"
    #for i in range(len(x1)):
#	print i,x1[i],y1[i], len(x1),len(y1)

    x2=N.array(x1)
    y2=N.array(y1)
    ppgplot.pgline(x2,y2)
def drawhistpylab(x,y,style): #draw histogram of binned data, x=left side of bins, y=number per bin
    x1=[]
    y1=[]
    n=len(x)
    dx=x[1]-x[0]
    for i in range(n):
        x1.append(x[i])
        x1.append(x[i])
    x1.append(x[(n-1)]+dx)
    x1.append(x[(n-1)]+dx)
    y1.append(0.)
    for i in range(n):
        y1.append(y[i])
        y1.append(y[i])
    y1.append(0.)
    #print len(x1),x1
    #print len(y1),y1
    xp=N.array(x1,'f')
    yp=N.zeros(len(y1),'f')
    for i in range(len(yp)):
	yp[i]=float(y1[i])

    #c=color
    #s=style
    pylab.plot(xp,yp,style)

def plotcumulative(x,Nbin,style,magplot=0):
    Nbin=len(x)
    p,bins,t=pylab.hist(x,Nbin,normed=True)
    db=bins[1]-bins[0]
    cdf=N.cumsum(p*db)

    x1=[]
    y1=[]
    n=len(bins)
    for i in range(n):
        x1.append(bins[i])
        x1.append(bins[i])
    x1.append(bins[(n-1)]+db)
    x1.append(bins[(n-1)]+db)
    y1.append(0.)
    for i in range(n):
        y1.append(cdf[i])
        y1.append(cdf[i])
    y1.append(1.)
    xp=N.array(x1,'f')
    yp=N.zeros(len(y1),'f')
    for i in range(len(yp)):
	yp[i]=float(y1[i])

    if magplot < 0.1:
	pylab.plot(xp-0.5*db,yp,style)
    else: #plot is magnitudes, so flip y axis
	pylab.plot(xp-0.5*db,1.-yp,style)

def plotnormhist(x,Nbin,style,magplot=0):
    #Nbin=len(x)
    p,bins,t=pylab.hist(x,Nbin,normed=True)
    db=bins[1]-bins[0]

    x1=[]
    y1=[]
    n=len(bins)
    for i in range(n):
        x1.append(bins[i])
        x1.append(bins[i])
    x1.append(bins[(n-1)]+db)
    x1.append(bins[(n-1)]+db)
    y1.append(0.)
    for i in range(n):
        y1.append(p[i])
        y1.append(p[i])
    y1.append(0.)
    xp=N.array(x1,'f')
    yp=N.zeros(len(y1),'f')
    for i in range(len(yp)):
	yp[i]=float(y1[i])

    if magplot < 0.1:
	pylab.plot(xp-0.5*db,yp,style)
    else: #plot is magnitudes, so flip y axis
	pylab.plot(xp-0.5*db,1.-yp,style)
def drawhistpylab2(x,y,style,lcolor): #draw histogram of binned data, x=left side of bins, y=number per bin
    x1=[]
    y1=[]
    n=len(x)
    dx=x[1]-x[0]
    for i in range(n):
        x1.append(x[i])
        x1.append(x[i])
    x1.append(x[(n-1)]+dx)
    x1.append(x[(n-1)]+dx)
    y1.append(0.)
    for i in range(n):
        y1.append(y[i])
        y1.append(y[i])
    y1.append(0.)
    #print len(x1),x1
    #print len(y1),y1
    xp=N.array(x1,'f')
    yp=N.zeros(len(y1),'f')
    for i in range(len(yp)):
	yp[i]=float(y1[i])

    #c=color
    #s=style
    pylab.plot(xp,yp,ls=style,color=lcolor)

def drawhistpylabnoplot(x,y,style): #draw histogram of binned data, x=left side of bins, y=number per bin
    x1=[]
    y1=[]
    n=len(x)
    dx=x[1]-x[0]
    for i in range(n):
        x1.append(x[i])
        x1.append(x[i])
    x1.append(x[(n-1)]+dx)
    x1.append(x[(n-1)]+dx)
    y1.append(0.)
    for i in range(n):
        y1.append(y[i])
        y1.append(y[i])
    y1.append(0.)
    #print len(x1),x1
    #print len(y1),y1
    xp=N.array(x1,'f')
    yp=N.zeros(len(y1),'f')
    for i in range(len(yp)):
	yp[i]=float(y1[i])

    #c=color
    #s=style
    return xp,yp


def cumulative(input):
    x=N.sort(input)
    n=len(input)
    y=N.arange(0,1,(1./n))
    return x,y

def median(x):
    x=N.sort(x)
    n=float(len(x))
    nmed=int(n/2.)
    return x[nmed]

def erroryp(x,y,err):
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsah(2,180.,1)
    for i in range(len(x)):
        ppgplot.pgarro(x[i],y[i],x[i],(y[i]+err[i]))
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)

def errorym(x,y,err):
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsah(2,180.,1)
    for i in range(len(x)):
        ppgplot.pgarro(x[i],y[i],x[i],(y[i]-err[i]))
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)

def errory(x,y,err):
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    ppgplot.pgsah(2,180.,1)
    for i in range(len(x)):
        ppgplot.pgarro(x[i],y[i],x[i],(y[i]-err[i]))
        ppgplot.pgarro(x[i],y[i],x[i],(y[i]+err[i]))
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)
    ppgplot.pgslw(ls)

def errorylog(x,y,err):
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    ppgplot.pgsah(2,180.,1)
    for i in range(len(x)):
	if (y[i] > 0.):
	    if ( (y[i] - err[i]) <= 0):
		ppgplot.pgarro(x[i],N.log10(y[i]),x[i],-99.)
		ppgplot.pgarro(x[i],N.log10(y[i]),x[i],N.log10(y[i]+err[i]))
	    else:
		ppgplot.pgarro(x[i],N.log10(y[i]),x[i],N.log10(y[i]-err[i]))
		ppgplot.pgarro(x[i],N.log10(y[i]),x[i],N.log10(y[i]+err[i]))
            
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)
    ppgplot.pgsls(ls)

def errorylogpt(x,y,err):#for one point only
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsls(1)
    ppgplot.pgsah(2,180.,1)
    if (y > 0.):
	if ( (y - err) <= 0):
	    ppgplot.pgarro(x,N.log10(y),x,-99.)
	    ppgplot.pgarro(x,N.log10(y),x,N.log10(y+err))
	else:
	    ppgplot.pgarro(x,N.log10(y),x,N.log10(y-err))
	    ppgplot.pgarro(x,N.log10(y),x,N.log10(y+err))
            
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)
    ppgplot.pgsls(ls)

def errorx(x,y,err):
    sch = ppgplot.pgqch()
    lw = ppgplot.pgqlw()
    ls = ppgplot.pgqls()
    ppgplot.pgsch(.5)
    ppgplot.pgslw(1)
    ppgplot.pgsah(2,180.,1)
    for i in range(len(x)):
        ppgplot.pgarro(x[i],y[i],(x[i]-err[i]),y[i])
        ppgplot.pgarro(x[i],y[i],(x[i]+err[i]),y[i])
    ppgplot.pgsch(sch)
    ppgplot.pgslw(lw)

def ratioerrorold(a,b):
    a=1.*a#make sure it's a real number, not integer
    b=1.*b
    try:
	len(a)#continue of a in an array
	ratio=[]
	err=[]
	for i in range(len(a)):
	    try:
		ratio.append(a[i]/b[i])
		if a[i] > 0.:
		    err.append(a[i]/b[i]*N.sqrt(1./a[i] + 1./b[i]))
		else:
		    err.append(0.)
	    except ZeroDivisionError:
		ratio.append(-1)
		err.append(0.)
    except TypeError:#a is a number rather than an array
	ratio=a/b
	try:
	    err = a/b*N.sqrt(1./a + 1./b)
	except:
	    err = 0.
    return ratio,err

def ratioerror(a,b):
    a=1.*a#make sure it's a real number, not integer
    b=1.*b
    try:
	len(a)#continue of a in an array
	ratio=[]
	errup=[]
	errdown=[]
	for i in range(len(a)):
	    try:
		ratio.append(a[i]/b[i])
		(edown,eup)=binomial.Binomial_lim(a[i],(b[i]-a[i]))
		errdown.append(edown)
		errup.append(eup)
	    except ZeroDivisionError:
		if a[i] < 0.1:
		    ratio.append(0.)
		    errup.append(0.)
		    errdown.append(0.)
		else:
		    ratio.append(-1)
		    errup.append(0.)
		    errdown.append(0.)

    except TypeError:#a is a number rather than an array
	try:
	    ratio=a/b
	    (edown,eup)=binomial.Binomial_lim(a,(b-a))
	    errup=eup
	    errdown=edown
	except ZeroDivisionError:
		ratio=-1
		errup=0.
		errdown=0.
    try:
	errdown=ratio-errdown
	errup=errup-ratio
    except:
	ratio=N.array(ratio,'f')
	errdown=N.array(errdown,'f')
	errup=N.array(errup,'f')
	errdown=ratio-errdown
	errup=errup-ratio

    return ratio,errdown,errup

def E(z):
    E = N.sqrt(OmegaL + OmegaM*(1.+z)**3)
    return E

def func(x):
    return 1.0/(N.sqrt(OmegaM*(1+x)**3 + OmegaL))

def H(z):
    H0=100
    Hz=H0*E(z)
    return Hz

def r200(sigma,z,h):
    #r200=1.73*sigma/1000/(E(z))*1./h
    r200=1.41*sigma/1000/(E(z))*1./h
    return r200

def DAold(z,h):#uses clunky integrator
    zstep=100
    dz=z/(1.*zstep)
    c=3.*10**5
    H0=100.*h

    t=0
    s=0
    zsteps=N.arange(0.,z+dz,dz)
    for zi in zsteps:
        #zi=dz*(i-1)
        t=t+dz/E(zi)/(1+zi)
        s=s+dz/E(zi)
    DA=c/H0/(1+z)*s #(Mpc/radian)
    DA=DA*1000./206264#(kpc/arcsec)    
    return DA

def DA(z,h):#uses numerical integration technique gaussian quadrature
    c=3.e5
    H0=100.*h
    #fields=scipy.integrate.quad(func,0,z)
    #s=fields[0]
    #errs=fields[1]
    s=romberg(func,0.,z)#,tol=1.e-6)
    #s=romberg(func,(0.,z))#,tol=1.e-6)
    DA=c/H0/(1+z)*s #(Mpc/radian)
    DA=DA*1000./206264#(kpc/arcsec)    
    return DA

def dLold(z,h):
    zstep=10000
    dz=z/(1.*zstep)
    c=3.*10**5
    H0=100.*h
    t=0
    s=0
    tH=9.78*h
    zsteps=N.arange(0.,z+dz,dz)
    #for i in range(zstep+1):
    for zi in zsteps:
        #zi=dz*(i-1)
        t=t+dz/E(zi)/(1+zi)
        s=s+dz/E(zi)
    DA=c/H0/(1+z)*s #(kpc/arcsec)
    DA=DA*1000./206264
    DL=c/H0*(1+z)*s #(Mpc/h)
    t=t*tH #lookbackt Gyr
    #print "DA = ",DA
    #print "DL = ",DL
    #print "lookback t = ",t
    #distmod=5*log10(DL*1.d6/10)
    #  print "Distance Modulus = ",distmod

    return DL

def dL(z,h):
    c=3.e5
    H0=100.*h
    try:#multiple objects
	n=len(z)
	DL=pylab.zeros(n,'d')
	for i in range(n):
	    s=romberg(func,0.,z[i])#,tol=1.e-6)
	    DL[i]=c/H0*(1+z[i])*s #(Mpc/h)
    except TypeError:
	s=romberg(func,0.,z)#,tol=1.e-6)
	DL=c/H0*(1+z)*s #(Mpc/h)
    return DL


def dLcm(z,h):
    DL=dL(z,h)*1.e6*3.08568025e18#convert to cm
    return DL

def lookbackt(z,h):
    c=3.e5
    H0=100.*h
    tH=9.78*h
    t=romberg(funct,0.,z)
    t=t*tH #time in Gyr
    return t,2*t

def funct(x):#for lookback time
    return 1.0/(N.sqrt(OmegaM*(1+x)**3 + OmegaL)*(1+x))


def cumulative(input):
    #print "max x = ",max(input)
    x=pylab.sort(input)
    #print "max x = ",max(x)
    n=len(input)
    y=pylab.arange(0.,1.,(1./float(n)))
    return x,y


def kcorr(z,band):#calculates the k-correction given a galaxy type t and redshift
    k=N.zeros(4,'f')#galaxy types E, Sbc, Scd, Im
    if(band == 1):
        infile=open('/Users/rfinn/bin/kcorr_F606W_Rc.fit','r')
    elif (band == 2):
        infile=open('/Users/rfinn/bin/kcorr_F675W_Rc.fit','r')
    elif (band == 3):
        infile=open('/Users/rfinn/bin/kcorr_F702W_Rc.fit','r')
    elif (band == 4):
        infile=open('/Users/rfinn/bin/kcorr_Rc_Rc.fit','r')
    elif (band == 5):
        infile=open('/Users/rfinn/bin/kcorr_V_Rc.fit','r')
    elif (band == 6):
        infile=open('/Users/rfinn/bin/kcorr_Ic_Rc.fit','r')
    elif (band == 7):
        infile=open('/Users/rfinn/bin/kcorr_F555W_Rc.fit','r')
    elif (band == 8):
        infile=open('/Users/rfinn/bin/kcorr_J_Rc.fit','r')

    for i in range(4):#evaluate k(i) at z, where i=t+1
        k[3-i]=0.
	sum=0
        for j in range(6):
            c=float(infile.readline())
            #k[3-i]=k[3-i]+c*(z**(j-1))
	    #print i,j,k[3-i],c,z,c*(z**(1.*j))
	    sum=sum+c*(z**(1.*j))
            #k[3-i] = k[3-i]+ c*(z**(1.*j))
            #k[3-i] += c*(z**(1.*j))
	    #print j,k[3-i],c,z
	k[3-i]=sum
    infile.close()


    #c..now solve for k(t,z) by interpolating between closest t values with known 
    #c..k-corr
    #c..n=1 => early type galaxy=E
    #if(n .eq. 1)then
    #kcorr=k(1) + (k(2)-k(1))*t
    #c..n=2 => intermediate type galaxy, Sbc
    #elseif(n .eq. 2)then
    #if(t .lt. 1)then
    #kcorr=k(1)+ (k(2)-k(1))*(t)
    #else
    #kcorr=k(2)+ (k(3)-k(2))*(t-1)
    #endif
    #c..Late type galaxy, Scd
    #elseif(n .eq. 3)then
    #if(t .lt. 2)then
    #kcorr=k(2)+ (k(3)-k(2))*(t-1)
    #elseif(t .lt. 2.5)then
    #kcorr=k(3)+ (k(4)-k(3))*(t-2)
    #c..late type galaxy, Im
    #else
    #kcorr=k(4)+ (k(4)-k(3))*(t-3)
    #endif
    #endif
    #for i in range(len(k)):
    #    print i, k[i]
    return k

def binyfixedbins(bins,x,y):
    ybin=pylab.zeros(len(bins),'f')

    dx=bins[1]-bins[0]

    xbinnumb=(x-bins[0])/dx#calculate x  bin number for each point 
    xbinnumb=pylab.array(xbinnumb,'f')
    for i in range(len(x)):
	j=int(xbinnumb[i])

	if j > (len(ybin)-1):
	    continue
	if j < 0:
	    continue
	try:
	    ybin[j]=ybin[j]+y[i]
	except IndexError:
	    print i,len(y),j,len(ybin)
    return ybin


def contourold(x,y,xmin,xmax,ymin,ymax):
    ncontbin=20.
    dx=float(abs((xmax-xmin)/ncontbin))
    dy=float(abs((ymax-ymin)/ncontbin))
    xcontour=N.arange(xmin,(xmax),dx)
    ycontour=N.arange(ymin,(ymax),dy)
    A=N.zeros((int(ncontbin),int(ncontbin)),'f')
    xbinnumb=N.array(len(x),'f')
    ybinnumb=N.array(len(x),'f')
    x1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),x)
    y1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*ncontbin/(xmax-xmin))#calculate x  bin number for each point 
    ybinnumb=((y-ymin)*ncontbin/(ymax-ymin))#calculate x  bin number for each point 
    binID = zip(xbinnumb,ybinnumb)
    Amax=0
    for (i1,j1) in binID:
        i=int(i1)
        j=int(j1)
	if i > (ncontbin-1):
	    continue
	if j > (ncontbin-1):
	    continue
	#print i,j,A
        A[j,i]=A[j,i]+1#have to switch i,j in A[] to make contours look right
    #pylab.figure()
    #print A
    for i in range(int(ncontbin)):
	for j in range(int(ncontbin)):
	    if A[j,i] > Amax:
		Amax=A[j,i]
    #Amax=max(max(A))
    print "max of A = ",Amax
    V=N.array([.025,.1,.25,.5,.75,.95],'f')
    V=V*Amax
    con=pylab.contour((xcontour+dx),(ycontour),A,V,alpha=5.,colors='r',linewidth=5.,hold='on')

def contour(x,y,xmin,xmax,ymin,ymax,color,nbin):
    ncontbin=nbin
    dx=float(abs((xmax-xmin)/ncontbin))
    dy=float(abs((ymax-ymin)/ncontbin))
    xcontour=pylab.arange(xmin,(xmax),dx)
    ycontour=pylab.arange(ymin,(ymax),dy)
    A=pylab.zeros((int(ncontbin),int(ncontbin)),'f')
    xbinnumb=pylab.array(len(x),'f')
    ybinnumb=pylab.array(len(x),'f')
    x1=pylab.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),x)
    y1=pylab.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),y)
    x=x1
    y=y1
    xbinnumb=((x-xmin)*ncontbin/(xmax-xmin))#calculate x  bin number for each point
    ybinnumb=((y-ymin)*ncontbin/(ymax-ymin))#calculate x  bin number for each point
    binID = zip(xbinnumb,ybinnumb)
    Amax=0
    for (i1,j1) in binID:
        i=int(i1)
        j=int(j1)
        if i > (ncontbin-1):
            continue
        if j > (ncontbin-1):
            continue
        #print i,j,A
        A[j,i]=A[j,i]+1#have to switch i,j in A[] to make contours look right
    #pylab.figure()
    #print A
    for i in range(int(ncontbin)):
        for j in range(int(ncontbin)):
            if A[j,i] > Amax:
                Amax=A[j,i]
    #Amax=max(max(A))
    print "max of A = ",Amax
    V=pylab.array([.025,.1,.25,.5,.75,.95],'f')
    V=V*Amax
    con=pylab.contour((xcontour+dx),(ycontour),A,V,alpha=5.,colors=color,linewidth=.5,hold='on')

def contourf(x,y,xmin,xmax,ymin,ymax):
    ncontbin=20.
    dx=float(abs((xmax-xmin)/ncontbin))
    dy=float(abs((ymax-ymin)/ncontbin))
    xcontour=N.arange(xmin,(xmax),dx)
    ycontour=N.arange(ymin,(ymax),dy)
    A=N.zeros((int(ncontbin),int(ncontbin)),'f')
    xbinnumb=N.array(len(x),'f')
    ybinnumb=N.array(len(x),'f')
    x1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),x)
    y1=N.compress((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*ncontbin/(xmax-xmin))#calculate x  bin number for each point 
    ybinnumb=((y-ymin)*ncontbin/(ymax-ymin))#calculate x  bin number for each point 
    binID = zip(xbinnumb,ybinnumb)
    Amax=0
    for (i1,j1) in binID:
        i=int(i1)
        j=int(j1)
	if i > (ncontbin-1):
	    continue
	if j > (ncontbin-1):
	    continue
	#print i,j,A
        A[j,i]=A[j,i]+1#have to switch i,j in A[] to make contours look right
    #pylab.figure()
    #print A
    for i in range(int(ncontbin)):
	for j in range(int(ncontbin)):
	    if A[j,i] > Amax:
		Amax=A[j,i]
    #Amax=max(max(A))
    print "max of A = ",Amax
    V=N.array([.01,.1,.25,.5,.75,.95],'f')
    V=V*Amax
    con=pylab.contourf((xcontour+dx),(ycontour),A,V,hold='on')
    #con.set_zorder(5.)



def findnearest(x1,y1,x2,y2,delta):
	dmin = 100000000000000000000000000000000000.
	matchflag=1
	nmatch=0
	for i in range(len(x2)):
		d = N.sqrt((x1-x2[i])**2 + (y1-y2[i])**2)
		if d < delta:
			nmatch=nmatch+1
		if d < dmin:
			dmin = d
			imatch = i

	
	if dmin > delta:
		imatch = 0
		matchflag = 0
	return imatch, matchflag,nmatch

def findnearestalt(x1,y1,x2,y2,delta):#use where command
	dmin = 100000000000000000000000000000000000.
	matchflag=1
	nmatch=0
	d= N.sqrt((x1-x2)**2 + (y1-y2)**2)#x2 and y2 are arrays
	t=pylab.where(d<delta)
	matches=t[0]
	if len(matches) > 0:
		nmatch=len(matches)
		z=d[matches]
		imatch=matches[pylab.where(z == z.min())]
		
	else:
		imatch = 0
		matchflag = 0


	return imatch, matchflag,nmatch
