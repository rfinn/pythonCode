#!/usr/bin/env python

import os

def getphot24():
	names=['K0225','K0847','K1038','K1516','K1759','K1953','K2042']
	srcid=['33','30','29','32','37','47','48']
	os.system('rm /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.long.dat')#prints entire line from apex output
	os.system('rm /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.dat')
	for i in range(len(names)):
		s="awk '{ if ($1 == "+srcid[i]+") print $0}' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips24/apex_output/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.long.dat"
		os.system(s)
		#s="awk '{ if ($1 == "+srcid[i]+") print $1, $3, $5, $13, $14 }' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips24/apex_output/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.dat"
                s="awk '{ if ($1 == "+srcid[i]+") print $13, $14 }' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips24/apex_output/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.dat"
		os.system(s)



def getphot70():
	names=['K0225','K0847','K1038','K1516','K1759','K1953','K2042']
	srcid=['13','-99','21','13','11','14','-99']
	os.system('rm /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.long.dat')#prints entire line from apex output
	os.system('rm /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.dat')
	for i in range(len(names)):
            if names[i].startswith('K0847'):
                        #s="echo -999 -9999999999 -999999999 6.400e+00 6.40e+00 >>  /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
                        s="echo 4.000e+03 4.00e+03 >>  /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.dat"
                        os.system(s)
            elif names[i].startswith('K2042'):
                        #s="echo -999 -9999999999 -999999999 6.400e+00 6.40e+00 >>  /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
                        s="echo 1.100e+03 1.10e+03 >>  /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.dat"
                        os.system(s)
            else:

		s="awk '{ if ($1 == "+srcid[i]+") print $0}' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips70/ch2/bcd/fpbcd/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.long.dat"
		os.system(s)
		#s="awk '{ if ($1 == "+srcid[i]+") print $1, $3, $5, $13, $14 }' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips24/apex_output/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips24.dat"
                s="awk '{ if ($1 == "+srcid[i]+") print $13, $14 }' < /Users/rfinn/research/KISS/data/"+names[i]+"/mips70/ch2/bcd/fpbcd/mosaic_extract.tbl >> /Users/rfinn/research/KISS/MIPSPhot/KISSmips70.dat"
		os.system(s)


def getphotirac():
	names=['K0225','K0847','K1038','K1516','K1759','K1953','K2042']
	srcidch1=['681','837','891','802','794','933','832']
	srcidch2=['1082','1078','1035','1103','1006','1104','1006']
	srcidch3=['1062','1056','1012','1050','1058','1053','1066']
	srcidch4=['1102','1153','1198','1120','1043','1089','0']


	for j in range(1,5):
		s='rm /Users/rfinn/research/KISS/IracPhot/KISSiracCh'+str(j)+'.long.dat'
		os.system(s)
		s='rm /Users/rfinn/research/KISS/IracPhot/KISSiracCh'+str(j)+'.dat'
		os.system(s)
		if j == 1:
			srcid=srcidch1
		if j == 2:
			srcid=srcidch2
		if j == 3:
			srcid=srcidch3
		if j == 4:
			srcid=srcidch4
		for i in range(len(names)):
                    s="awk '{ if ($1 == "+srcid[i]+") print $0}' < /Users/rfinn/research/KISS/data/"+names[i]+"/irac/ch"+str(j)+"/bcd/output/extract.tbl >> /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".long.dat"
                    os.system(s)

                    if srcid[i].startswith('0'):
                        #s="echo -999 -9999999999 -999999999 6.400e+00 6.40e+00 >>  /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
                        s="echo 6.400e+00 6.40e+00 >>  /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
                        os.system(s)
                    else:
			#s="awk '{ if ($1 == "+srcid[i]+") print $1,$4,$6,$14,$15}' < /Users/rfinn/research/KISS/data/"+names[i]+"/irac/ch"+str(j)+"/bcd/output/extract.tbl >> /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
                        s="awk '{ if ($1 == "+srcid[i]+") print $14,$15}' < /Users/rfinn/research/KISS/data/"+names[i]+"/irac/ch"+str(j)+"/bcd/output/extract.tbl >> /Users/rfinn/research/KISS/IracPhot/KISSiracCh"+str(j)+".dat"
			os.system(s)


getphot24()
getphot70()
getphotirac()
infile=open('/Users/rfinn/research/KISS/IracPhot/KISSiracCh1.dat','r')
ch1=infile.readlines()
infile.close()
infile=open('/Users/rfinn/research/KISS/IracPhot/KISSiracCh2.dat','r')
ch2=infile.readlines()
infile.close()

infile=open('/Users/rfinn/research/KISS/IracPhot/KISSiracCh3.dat','r')
ch3=infile.readlines()
infile.close()

infile=open('/Users/rfinn/research/KISS/IracPhot/KISSiracCh4.dat','r')
ch4=infile.readlines()
infile.close()

infile=open('/Users/rfinn/research/KISS/MIPSPhot/KISSmips24.dat','r')
ch24=infile.readlines()
infile.close()

infile=open('/Users/rfinn/research/KISS/MIPSPhot/KISSmips70.dat','r')
ch70=infile.readlines()
infile.close()

outfile=open('/Users/rfinn/research/KISS/KISSallphot.dat','w')
for i in range(len(ch1)):
    outline=ch1[i].rstrip()+' '+ch2[i].rstrip()+' '+ch3[i].rstrip()+' '+ch4[i].rstrip()+' '+ch24[i].rstrip()+' '+ch70[i]
    outfile.write(outline)
outfile.close()
