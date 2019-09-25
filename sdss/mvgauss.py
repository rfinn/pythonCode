#!/usr/bin/env python
import sys, os
i=0
input = open('gaussian-clusters','r')
for line in input:
    if line.find('#') > -1:
        continue
    fields=line.split()
    name=fields[0]
    i += 1
    cluster="cluster_DR3_"+str(name)+".dat"
    command = "cp "+cluster+" /home/rfinn/SDSS/balogh1DR3/ \n"
    os.system(command)
input.close()
print "number of gaussian clusters = ",i
