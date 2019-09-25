#!/usr/bin/env python
'''
writing code to do a quick analysis of the depth of the 09/2009 data

RXJ1821 looks pretty good

RXJ1716 is compromised by two bright stars on the field

'''
from z08common import *
from pylab import *
import os
import asciidata
import atpy


RAlist=['13:17','17:16','13:50']
def renamearchivestacks():
    os.system('gethead tu*.fits OBJRA OBJDEC FILTER NCOMBINE PRODTYPE DATE-OBS MAGZERO MAGZSIG SEEING SKYMAG > header_info')
    infile=open('header_info','r')
    for line in infile:
        t=line.split()
        r=t[6].split(':')
        fmt=t[0].split('.')
        if t[1].find('13:17') > -1:
            clname='cl1317'
        elif t[1].find('13:50') > -1:
            clname='cl1350'
        elif t[1].find('17:16') > -1:
            clname='cl1716'

        elif t[1].find('10:00') > -1:
            clname='cosmos'
        elif t[1].find('18:21') > -1:
            clname='cl1821'


        if len(fmt) == 2:
            suffix=fmt[1]
        elif len(fmt) == 3:
            suffix=str(fmt[1])+'.'+str(fmt[2])
        s='ln -s '+t[0]+' '+clname+'-'+t[3]+'-'+t[5]+'-'+r[0]+'h.'+suffix
        print s
        os.system(s)

renamearchivestacks()
