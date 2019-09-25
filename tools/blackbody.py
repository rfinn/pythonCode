#!/scisoft/bin/python
#generate two blackbody curves for astro exam question
import ppgplot
import Numeric as N
import mystuff as my


c=3.e10#cm/s
h=6.626e-26#erg s
k=1.2807e-16#erg K

l=N.arange(2000.,50000.,1000.)#wavelength in A
l=l*1.e-8#convert to cm

T=40000.#K
B=2.*h*c**2/(l**5)/(N.exp(h*c/(l*k*T))-1)/1.e14
l=l*1.e4
xmin=1.15*min(l)
xmax=max(l)
ymin=min(B)
ymax=1.2*max(B)
my.psplotinit("blackbody.ps")
ppgplot.pgbox("",0.0,0,"",0.0,0)
ppgplot.pgenv(xmin,xmax,ymin,ymax,0,0)
ppgplot.pglab("Wavelength","Energy Output/second","")
ppgplot.pgsci(4)
ppgplot.pgline(l,B)
ppgplot.pgtext(.6,3.,'Star A')
ppgplot.pgtext(.6,-.5,'Blue')
#T=20000.#K
#B=2.*h*c**2/(l**5)/(N.exp(h*c/(l*k*T))-1)
l=(l*1.e-4+20000e-8)*1.e4
ppgplot.pgsci(2)
ppgplot.pgline(l,B)
ppgplot.pgtext(2.6,3.,'Star B')
ppgplot.pgtext(3.6,-.5,'Red')
ppgplot.pgsci(1)

ppgplot.pgend()

    
