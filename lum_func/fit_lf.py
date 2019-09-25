#!/usr/bin/env python

import pylab as pl
import numpy as np
import atpy
from scipy.optimize import leastsq
#from mpfit import mpfit
import mpfitexpr
# try to fit w/scipy optimize


#schechter = lambda p, x,dx:p[0]*(x/p[1])**(-1*p[2])*exp(-1*x/p[1])#*dx/p[1]

#schechter_errfunc_weighted = lambda p, x, dx, y, erry,weight: weight*(y - schecter(p,x,dx))/erry



#fitfuncdpl=lambda p, x: log(10.)*10.**x*(p[2]/10.**x+p[3]/(10.**x+10.**p[1]))*(p[0]*(10.**(x-p[1]))**(-p[2])*(1.+10.**(x-p[1]))**(-p[3]))
#errfuncdpl=lambda p, x, y, yerr: (y - fitfuncdpl(p, x))**2/yerr**2

weightedErrFuncDpl=lambda p, x, y, yerr, weights: weights*((y - fitfuncdpl(p, x))**2/yerr**2)


def schechter(p,x):
    phistar = p[0]
    Lstar = 10.**p[1]
    alpha = p[2]
    y=phistar*(x/Lstar)**alpha*np.exp(-1.*x/Lstar)
    return y
def schechterlogL(p,x):
    phistar = p[0]
    Lstar = 10.**p[1]
    alpha = p[2]
    y=phistar*(10.**x/Lstar)**(alpha+1)*np.exp(-1.*10.**x/Lstar)
    return y

#schechterlogL_errfunc = lambda p, x, y, erry: (y - schechterlogL(x,p[0],p[1],p[2]))/erry
schechterlogL_string='p[0]*(10.**x/10.**p[1])**(p[2]+1)*numpy.exp(-10.**x/10.**p[1])'
schechterlogL_fixed_alpha_string='p[0]*(10.**x/10.**p[1])**(-0.41)*numpy.exp(-10.**x/10.**p[1])'

schechterlogL_fixed_LCSalpha_string='p[0]*(10.**x/10.**p[1])**(-0.56)*numpy.exp(-10.**x/10.**p[1])'

schechterlogL_bai_string='p[0]*(10.**(x-10.5))**(-0.41)*numpy.exp(-10.**(x-10.5))'
bai_lf = lambda p,x: p[0]*((10.**(x-10.5))**(-0.41)*np.exp(-10.**(x-10.5)))
normalize_bai = lambda p, x, y, erry: (y - bai_lf(p,x))/erry
schechterlogL_errfunc = lambda p, x, y, erry: (y - schechterlogL(p,x))/erry

def binlogL(x,w,xmin,xmax,nbin,sigmaflag=0,sigma_weight=None):#use equally spaced bins
    ''' 
    bin x array into nbin equally-spaced bins that range 
    from xmin to xmax

    weights (e.g. spectroscopic completeness) are passed in using w

    can also normalize by sigma using sigmaflag=1 and sigmaweight = 1/sigma^3 or 1/area

    '''
    #print 'weights from within binlogL', w
    if len(x) != len(w):
        print 'array and weights must be same length'
        return
    #print 'xmin, xmax, nbin = ',xmin,xmax,nbin
    dx=float((xmax-xmin)/(1.*nbin))
    
    xbin=np.arange(xmin,(xmax),dx)+dx/2.
    
    ybin=np.zeros(len(xbin),'d')
    ybinerr=np.zeros(len(xbin),'d')
    xbinnumb=np.array(len(x),'d')
    
    #x=x[((x >= xmin) & (x <= xmax))]
    #print 'x = ',x

    if sigmaflag:
        w=w*sigma_weight
        #print 'sigma_weight'
        #print 'in sigmaflag ',sigma_weight
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))-1#calculate x  bin number for each point
    #print 'xbinnumb = ',xbinnumb
    for i in range(len(xbin)):
        xbindata=x[(abs(xbinnumb-float(i))<.5)]
        xbindata=np.ones(len(x[(abs(xbinnumb-float(i))<.5)]),'f')
        wbindata=w[(abs(xbinnumb-float(i))<.5)]
        #print xbin[i],' wbindata = ',wbindata
        xscaled=xbindata*wbindata
        if sigmaflag:
            #ybinerr[i]=np.sqrt(sum(xscaled))/sigma_weight[(abs(xbinnumb-float(i))<.5)])
            xscaled=xscaled*sigma_weight[(abs(xbinnumb-float(i))<.5)]
            ybin[i]=sum(xscaled)
            ybinerr[i]=np.sqrt(sum((xbindata*wbindata*sigma_weight[(abs(xbinnumb-float(i))<.5)])**2))
        #print 'bin = ',xbin[i],' scale = ',sum(xscaled)/sum(xbindata)
        else:
            ybin[i]=sum(xscaled)
            ybinerr[i]=np.sqrt(sum(xbindata**2*wbindata))
            
    print 'from w/in binlogL ',ybin
    return xbin,ybin,ybinerr

def fit_double_power_unbinned(logL,phiStar,logLstar,alpha,beta):
    ### Make histogram to find bin heights for fit
    hs1=hist(logL,bins=15) ###HISTOGRAM
    clf()
    ### Retrieve bin heights and centers from histogram data
    yvals=hs1[0]
    xcenters=[]
    for i in range(len(yvals)):
        xcenters.append((hs1[1][i]+hs1[1][1+i])/2.)
    logy=log10(yvals)
    yerr=sqrt(yvals) ###Still need to decide if this is the best way to do this
    ### Calculate best fit papamaters using double power law
    p0=[phiStar,logLstar,alpha,beta] ###[phiStar, log10(Lstar), alpha, betta]
    boot,success=fit_double_power_binned(xcenters,yvals,yerr,p0[0],p0[1],p0[2],p0[3])
    return boot, success
        
def fit_double_power_binned(x,y,yerr,phiStar,logLstar,alpha,beta):
    p0=[phiStar,logLstar,alpha,beta] ###[phiStar, log10(Lstar), alpha, betta]
    #Remove zero values
    flagZero=(y>0)
    y=y[where(flagZero)]
    x=array(x).T
    x=x[where(flagZero)]
    yerr=yerr[where(flagZero)]
    #boot,success=leastsq(errfuncdpl,p0[:],args=(array(x),y,yerr))
    boot,success,a,b,c=leastsq(errfuncdpl,p0[:],args=(array(x),y,yerr),full_output=1)
    return boot, success

def fit_double_power_with_weights(x,y,yerr,phiStar,logLstar,alpha,beta,weights):
    """
    INSERT DEBBIE/MIKE'S WEIGHTED HISTOGRAM CODE HERE:
    Outputs -
    yvals (counts - bin heights)
    xcenters (luminosities - bin centers)
    yerr (error in y)
    """
    p0=[phiStar,logLstar,alpha,beta]
    flagZero=(y>0)
    yvals=yvals[where(flagZero)]
    xcenters=xcenters[where(flagZero)]
    yerr=yerr[where(flagZero)]
    boot,success,a,b,c=leastsq(errfuncdpl,p0[:],args=(array(xcenters),yvals,yerr),full_output=1) #currently this is using the error func that doesn't take weights into account
    return boot, success
    
def fit_schechter(logL,weight,phistar,logLstar,alpha,sigmaflag=0,sweight=None,nbin=5,xmin=8.5,xmax=10.5,plotsingle=1,baiflag=1,prefix='LCS '):
    ''' takes unbinned data '''
    #print logL,weight
    #print xmin,xmax,nbin
    #print 'weight from w/in fit_schechter',weight
    xmax=max(logL)

    xbin,ybin,ybinerr=binlogL(logL,weight,xmin,xmax,nbin,sigmaflag=sigmaflag,sigma_weight=sweight)
    if plotsingle:
        pl.figure()
    sxbin,sybin,sybinerr=binlogL(logL,np.ones(len(logL),'f'),xmin,xmax,nbin,sigmaflag=sigmaflag,sigma_weight=sweight)
    if prefix.find('EDisCS') > -1:
        pl.plot(sxbin,sybin,'bo',label=prefix+'raw counts',markersize=12,mfc='None',mec='b')
        pl.plot(xbin,ybin,'bo',label=prefix+'corrected counts')
        pl.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='b',markersize=9)
    elif prefix.find('LCS') > -1:
        pl.plot(sxbin,sybin,'ro',label=prefix+'raw counts',markersize=12,mfc='None',mec='r')
        pl.plot(xbin,ybin,'ro',label=prefix+'corrected counts')
        pl.errorbar(xbin,ybin,yerr=ybinerr,fmt=None,ecolor='r',markersize=9)

    pl.gca().set_yscale('log')
    print 'DATA TO FIT:'
    for i in range(len(xbin)): print xbin[i],ybin[i]

    p0=[phistar,logLstar,alpha] ###[phiStar, log10(Lstar), alpha, betta]
    params,success,a,b,c=leastsq(schechterlogL_errfunc,p0[:],args=(np.array(xbin),ybin,ybinerr),full_output=1)    
    #print params,success
    #pl.plot(xbin,schechterlogL(params,xbin)*1.05,label='leastsq',color='r')

    dx=.01
    xl=np.arange(min(xbin),max(xbin)+dx,dx)
    #print 'xl = ',xl
    #xl=np.arange(7.75,10.75+dx,dx)
    t=mpfitexpr.mpfitexpr(schechterlogL_string, xbin, ybin, ybinerr , p0)
    params=t[0]
    s='mpfit (%5.2f, %5.2f, %5.2f)'%(params[0],params[1],params[2])
    if prefix.find('LCS') > -1:
        pl.plot(xl,schechterlogL(t[0],xl),color='r',label=prefix+s )
    #elif prefix.find('EDisCS') > -1:
    #    pl.plot(xl,schechterlogL(t[0],xl),color='r',label=prefix+s )
    #pl.plot(xbin,t[1],label=s,color='r')

    p1=[phistar,logLstar]
    t=mpfitexpr.mpfitexpr(schechterlogL_fixed_alpha_string, xbin, ybin, ybinerr , p1)
    params=t[0]
    s='mpfit (%5.2f, %5.2f,[-1.41])'%(params[0],params[1])
    pplot=[params[0],params[1],alpha]
    pl.plot(xl,schechterlogL(pplot,xl),color='k',ls='--',label=prefix+s )

    if prefix.find('EDisCS') > -1:
        p1=[phistar,logLstar]
        t=mpfitexpr.mpfitexpr(schechterlogL_fixed_LCSalpha_string, xbin, ybin, ybinerr , p1)
        params=t[0]
        s='mpfit (%5.2f, %5.2f,[-1.56])'%(params[0],params[1])
        pplot=[params[0],params[1],alpha]
        pl.plot(xl,schechterlogL(pplot,xl),color='b',ls='-',label=prefix+s )


    if baiflag:
        bai_alpha = -1.41
        bai_logLstar = 10.5
        bai_phistar = 6.5

        pbai_normalization=[bai_phistar]#,logLstar,alpha]
        params,success,a,b,c=leastsq(normalize_bai,pbai_normalization[:],args=(np.array(xbin),ybin,ybinerr),full_output=1)    

        s='Bai+09 (%5.2f, %5.2f, %5.2f)'%(params[0],logLstar,alpha)
        pbai=[params[0],logLstar,alpha]

        pl.plot(xl,schechterlogL(pbai,xl),label=s,color='c')

    pl.legend(numpoints=1,loc='lower left',prop={'size':12})
    pl.xlabel('$ log_{10}(L_{IR}/L_\odot) $',fontsize=20)

    return params,success,t
