#!/usr/bin/evn python

from pylab import *

def betai(a,b,x):
#      REAL betai,a,b,x
#CU    USES betacf,gammln
#      REAL bt,betacf,gammln
    if ((x < 0) or (x >1.)):
	print 'bad argument x in betai'
	return
    if ((x == 0) or (x == 1)):
        bt=0.
    else:
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
    if (x < (a+1.)/(a+b+2.)):
        betai=bt*betacf(a,b,x)/a
        return betai
    else:
        betai=1.-bt*betacf(b,a,1.-x)/b
        return betai


def betacf(a,b,x):
    MAXIT=100
    EPS=3.e-7
    FPMIN=1.e-30
    qab=a+b
    qap=a+1.
    qam=a-1.
    c=1.
    d=1.-qab*x/qap
    if(abs(d)< FPMIN):
	d=FPMIN
    d=1./d
    h=d
    for m in range(1,(MAXIT+1)):
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d) < FPMIN):
	    d=FPMIN
	c=1.+aa/c
        if(abs(c) < FPMIN):
	    c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d) < FPMIN):
	    d=FPMIN
        c=1.+aa/c
        if(abs(c) < FPMIN):
	    c=FPMIN
        d=1./d
        delt=d*c
        h=h*delt
        if(abs(delt-1.) < EPS):
	    betacf=h
	    return betacf
    print 'a or b too big, or MAXIT too small in betacf'
    raw_input('Press Enter to exit')


def gammln(xx):
#      REAL gammln,xx
#      INTEGER j
#      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
#      SAVE cof,stp
      cof=array([76.18009172947146e0,-86.50532032941677e0,24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2,-.5395239384953e-5],'d')
      stp=2.5066282746310005e0
      x=xx
      y=x
      tmp=x+5.5e0
      tmp=(x+0.5e0)*log(tmp)-tmp
      ser=1.000000000190015e0
      for j in range(6):
        y=y+1.e0
        ser=ser+cof[j]/y
      gammln=tmp+log(stp*ser/x)
      return gammln



def ks2d2s(x1,y1,n1,x2,y2,n2):
#    INTEGER n1,n2
#    REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
#    CU    USES pearsn,probks,quadct
#    INTEGER j
#    REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen,probks
    d1=0.0
    for j in range(n1):

	(fa,fb,fc,fd)=quadct(x1[j],y1[j],x1,y1,n1)


	(ga,gb,gc,gd)=quadct(x1[j],y1[j],x2,y2,n2)

	if (fa > ga):
	    fa += 1.0/n1
	if (fb > gb):
	    fb += 1.0/n1
	if (fc > gc): 
	    fc += 1.0/n1
	if (fd > gd):
	    fd += 1.0/n1


	if (ga > fa):
	    ga += 1.0/n2
	if (gb > fb):
	    gb += 1.0/n2
	if (gc > fc): 
	    gc += 1.0/n2
	if (gd > fd):
	    gd += 1.0/n2


	d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
    d2=0.0
    for j in range(n2):
        (fa,fb,fc,fd)=quadct(x2[j],y2[j],x1,y1,n1)
        (ga,gb,gc,gd)=quadct(x2[j],y2[j],x2,y2,n2)
        d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))

    d=0.5*(d1+d2)
    sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
    (r1,dum,dumm)=pearsn(x1,y1,n1)
    (r2,dum,dumm)=pearsn(x2,y2,n2)
    rr=sqrt(1.0-0.5*(r1**2+r2**2))
    prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
    return d,prob



def pearsn(x,y,n):
#      INTEGER n
#      REAL prob,r,z,x(n),y(n),TINY
      TINY=1.e-20
#    USES betai
#      INTEGER j
#      REAL ax,ay,df,sxx,sxy,syy,t,xt,yt,betai
      ax=0.
      ay=0.
      for j in range(n):
	  ax=ax+x[j]
	  ay=ay+y[j]
      ax=ax/n
      ay=ay/n
      sxx=0.
      syy=0.
      sxy=0.
      for j in range(n):
        xt=x[j]-ax
        yt=y[j]-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
      r=sxy/sqrt(sxx*syy)
      z=0.5*log(((1.+r)+TINY)/((1.-r)+TINY))
      df=n-2
      t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
      prob=betai(0.5*df,0.5,df/(df+t**2))
#     prob=erfcc(abs(z*sqrt(n-1.))/1.4142136)
      return r, prob,z


def probks(alam):
#      REAL probks,alam,EPS1,EPS2

#      INTEGER j
#      REAL a2,fac,term,termbf
#      PARAMETER (EPS1=0.001, EPS2=1.e-8)
    EPS1=0.001
    EPS2=1.e-8
    a2=-2.*alam**2
    fac=2.
    probks=0.
    termbf=0.
    for j in range(1,101):
        term=fac*exp(a2*j**2)
        probks=probks+term
        if((abs(term) < EPS1*termbf) or (abs(term) < EPS2*probks)):
	    return probks
        fac=-fac
        termbf=abs(term)
    probks=1.
    return probks

def quadct(x,y,xx,yy,nn):
#      INTEGER nn
#      REAL fa,fb,fc,fd,x,y,xx(nn),yy(nn)
#      INTEGER k,na,nb,nc,nd
#      REAL ff
    na=0
    nb=0
    nc=0
    nd=0
    for k in range(nn):
        if(yy[k] > y):
	    if(xx[k] > x):
		na=na+1
	    else:
		nb=nb+1

        else:
	    if(xx[k]  > x):
		nd=nd+1
	    else:
		nc=nc+1
    ff=1.0/nn
    fa=ff*na
    fb=ff*nb
    fc=ff*nc
    fd=ff*nd
    return fa,fb,fc,fd

