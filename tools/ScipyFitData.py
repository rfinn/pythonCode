from pylab import *
from scipy import optimize

from pylab import *
from scipy import *

# define function for calculating powerlaw



##########
# Fitting the data -- Least Squares Method
##########

# Power-law fitting is best done by first converting
# to a linear equation and then fitting to a straight line.
#
#  y = a * x^b
#  log(y) = log(a) + b*log(x)
#
fitfunc = lambda p, x: p[0] + p[1]*x
errfunc = lambda p, x, y, err: (y-fitfunc(p,x))/err
pinit=[1.0,-1.0]

out=optimize.leastsq(errfunc,pinit,args=(logx,logy,logyerr),full_output=1)
