#!/usr/bin/env python
from pylab import *
from matplotlib.backends.backend_pdf import PdfFile

pdf=PdfFile('testmultipage.pdf')

#figure(1)
x=arange(0.,10.)
y=x
plot(x,y**2,'bo')
#savefig('test.eps')
#savefig(pdf,format='pdf')
savefig(pdf)
#plot(x,y**3,'ro')
#savefig(pdf)#,format='pdf')
pdf.close()



