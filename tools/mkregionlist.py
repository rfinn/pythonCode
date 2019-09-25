#!/scisoft/bin/python
"""
useage is
mkregionlist.py inputfile outputfile

inputfile lists ediscs id name
"""
import sys, glob
import Numeric as N
from math import *

ra=[]
dec=[]
infile=(sys.argv[1])#index 1 b/c python filename is 0
outfile=(sys.argv[2])#index 1 b/c python filename is 0
input=open(infile,'r')
output=open(outfile,'w')
output.write("global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\nfk5\n")

for line in input:
    if line.find('#') > -1:
        continue
    t=line.split()
    shortname=t[0]
    #ra=((float(shortname[0:2])+float(shortname[2:4])/60.+float(shortname[4:7])/10./3660.)*15.)#convert from hours to degrees
    #ra=((float(shortname[0:2])+float(shortname[2:4])/60.+(float(shortname[4:7])/10.+.63)/3660.)*15.)#to align w/mips
    #print float(shortname[7:10]),float(shortname[10:12]),float(shortname[12:])
    #if float(shortname[7:10]) < 0:
    #    dec=float(shortname[7:10])-float(shortname[10:12])/60.-(float(shortname[12:])/10.-15.)/3660#to align w/vlt
        #dec=float(shortname[7:10])-float(shortname[10:12])/60.-(float(shortname[12:])/10.-0.)/3660#to align w/mips
    #else:
    #    dec=float(shortname[7:10])+float(shortname[10:12])/60.+float(shortname[12:])/10./3660
    ra=float(t[1])
    dec=float(t[2])
    if (float(t[4]) > 100.):
	string = "circle(%12.8f, %12.8f, 3\") \n"%(ra,dec)
	output.write(string)
    #if float(t[16]) > 0:
    #    string = "circle(%12.8f, %12.8f, 2\") # color=blue\n"%(ra,dec)
    #    output.write(string)
    #    string = "# text(%12.8f, %12.8f) text={"%(ra,dec)
    #    string=string+str(t[14])+"} \n"
    #    output.write(string)
input.close()
output.close()

