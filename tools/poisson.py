#!/usr/bin/env python
import numpy as N

###########################################################################
#for a given number of events calculate the Poisson upper and lower
#limits for the 1-sided confidence interval.  Formulas are (10) and
#((14) from Gehrels, 1986, ApJ, 303, 336
#0.8413, 0.9772, and 0.9987 correspond to the usual 1-sided intervals
#for 1, 2, and 3 sigma gaussian probabilities
def Poisson_lim(n1):

  n_i = n1#shift #ask greg what these are doing
  conflim_d = .84#shift         

  #my (imatch, iconf)

  #my (diff_d, conf_up, conf_dn)

  #allow only these limits
  conf_d = N.array([0.8413, 0.90, 0.95, 0.975, 0.9772, 0.990, 0.995, 0.9987, 0.999, 0.9995],'f')

  #'S' parameter for these limits
  S_d = N.array([1.0, 1.282, 1.645, 1.960, 2.00, 2.326, 2.576, 3.000, 3.090, 3.291],'f')

  #other parameters for lower limit
  beta_d = N.array([0.0, 0.010, 0.031, 0.058, 0.062, 0.103, 0.141, 0.222, 0.241, 0.287],'f')
  gamma_d = N.array([0.0, -4.0, -2.50, -2.22, -2.19, -2.07, -2.0, -1.88, -1.85, -1.80],'f')

  #in this version just use the value of @conf_d that's closest to the
  #input confidence limit
  #find closest value
  diff_d = 1.e2
  for iconf in range(len(conf_d)):
      delta = abs(conf_d[iconf]-conflim_d)
      if (delta < diff_d):
	  diff_d = delta
	  imatch = iconf

#  for(iconf = 0, diff_d = 1.e2 iconf < @conf_d iconf++) {#ask Greg what this is doing
#    if(abs(conf_d[iconf] - conflim_d) < diff_d) {
#      diff_d = abs(conf_d[iconf] - conflim_d)
#      imatch = iconf
#    }
#  }

  conf_up = n_i + S_d[imatch] * N.sqrt(n_i + 1) + (S_d[imatch]**2 + 2) / 3
  try:
      t=len(n_i)
      conf_dn=N.zeros(len(n_i))
      for i in range(len(n_i)):
	  if (n_i[i] < 1):
	      conf_dn[i] = 0
	  else:
	      conf_dn[i] = n_i[i] * (1 - 1 / (9 * n_i[i]) - S_d[imatch] / (3 *N.sqrt(n_i[i])) + beta_d[imatch] * n_i[i]**gamma_d[imatch])**3

  except TypeError:
      if (n_i < 1):
	  conf_dn = 0
      else:
	  conf_dn = n_i * (1 - 1 / (9 * n_i) - S_d[imatch] / (3 *N.sqrt(n_i)) + beta_d[imatch] * n_i**gamma_d[imatch])**3
  
  return conf_dn, conf_up

