#!/usr/bin/env python
#title           :GalaxyDensityPlot.py
#description     :This will create a couple of plots for Rose
#author          :Graziano Vernizzi, Rose Finn (Siena college)
#date            :Oct 13, 2016
#version         :0.1
#notes           :Copyright (C) 2016, 2017 Graziano Vernizzi and Rose Finn
#                 All rights reserved
#python_version  :  
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

#to read "fits" files
from astropy.io import fits

##### INPUT PARAMETERS
Nbins_phi=15  #number of bins
Nbins_theta=15
Nbins_RS=3  #     
Pick_RS_bin=1  #(labels from 0,1,2,...Nbins_RS-1)


catalogpath = '/Users/rfinn/research/LocalClusters/NSAmastertables/'
#Read the data file
infile=catalogpath+'Coma_NSA.fits'
coma = fits.getdata(infile)

# positions
phi=coma.RA # projection of longitude, angle 
theta=coma.DEC # projection of latitude, angle above equator
rs=coma.Z # redshift

###FOR A QUICK 3D scatter plot
fig = plt.figure()
ax = fig.gca(projection='3d')

#Spherical coordinates (unit sphere)
psi=np.pi/2-theta
X=np.cos(phi)*np.sin(psi)
Y=np.sin(phi)*np.sin(psi)
Z=np.cos(psi)

X = coma.RA
Y = coma.DEC
Z = coma.Z
#scatter plot
ax.scatter(X,Y,Z)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Quick projection onto the unit sphere')
plt.show()
##### End of  Quick 3D scatter plot


#slice the dataset along RS (according to Nbins_RS)
step_rs=(max(rs)-min(rs))/Nbins_RS

From_rs=step_rs*Pick_RS_bin+min(rs)
To_rs=From_rs+step_rs

#redefine phi and theta accordingly
phi2=[]; theta2=[]
for i in range(len(rs)):
    if ( rs[i]<To_rs and rs[i]>=From_rs):
        phi2.append(phi[i])
        theta2.append(theta[i])

    
#Count how many galaxies are in each rectangle of the grid
glx, x, y= np.histogram2d(phi2,theta2, bins=[Nbins_phi,Nbins_theta])
glx2=glx.T[::-1]  #from matrix indices (row, colummn)to x,y axis format

###Compute the volume of each grid elements
step_x=x[2]-x[1]
step_y=y[2]-y[1]
#Create a meshgrid for the densityplot
X,Y=np.meshgrid(x[:-1]+step_x/2,y[:-1]+step_y/2)
Volumes=abs(cos(Y/180*pi)*step_x*step_y*step_rs)

#The density is the number of galaxies in each volume
density=glx2/Volumes

#you can normalize the density to 1 with
density=density/sum(density)

##Plot everything
fig = plt.figure()
## Color map: jet, gray, binary, spectral, hsv, hot, cool, cooper, bone, autumn, winter, summer, spring, etc.
im = imshow(density, extent=[min(phi),max(phi),min(theta),max(theta)], cmap=cm.hot)
im.set_interpolation('bilinear')
xlabel('RA (degrees)'); ylabel('DEC (degrees)'); title('Number of galaxies/volume')
grid(True)
colorbar(im)

# To save the figure .  Possible extensions : png, pdf, ps, eps and svg :
plotfilename = 'densityFITS.png'
savefig(plotfilename,dpi=300)
#or show() if you don't want to save the figure

##Same plot as a surface in 3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
fig = plt.figure()
ax = fig.gca(projection='3d')

#for a surf plot
surf = ax.plot_surface(X, Y, density, rstride=1, cstride=1,cmap=cm.jet,linewidth=0, antialiased=True)
#ax.set_zlim(0, 120)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_xlabel('RA'); ax.set_ylabel('DEC'); ax.set_zlabel('number density')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
