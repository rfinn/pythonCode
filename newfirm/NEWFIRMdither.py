#!/usr/bin/env python

'''
creating a script to visualize the dither sequence that Janice and Chun
developed to provide equal illumination when observing w/NEWFIRM.
This instrument is made of 4 arrays and has gaps between the arrays.

'''
from pylab import *
from z08common import *

def drawbox(data,style):#feed in center x,y,dx,dy,rotation E of N
    #xcoords of unrotated box, going around CCW
    xl=array([data[0]-0.5*data[2],data[0]+0.5*data[2],data[0]+0.5*data[2],data[0]-0.5*data[2],data[0]-0.5*data[2]],'d')
    yl=array([data[1]-0.5*data[3],data[1]-0.5*data[3],data[1]+0.5*data[3],data[1]+0.5*data[3],data[1]-0.5*data[3] ],'d')

    xl=array([-0.5*data[2],+0.5*data[2],+0.5*data[2],-0.5*data[2],-0.5*data[2]],'d')
    yl=array([-0.5*data[3],-0.5*data[3],+0.5*data[3],+0.5*data[3],-0.5*data[3] ],'d')

    ang=data[4]*pi/180.*-1.#convert rotation to radians
    #rotate coordinates
    xp=cos(ang)*xl-sin(ang)*yl
    yp=sin(ang)*xl+cos(ang)*yl

    #put back on absolute scale
    xp=data[0]+xp
    yp=data[1]+yp
    #draw rotated box
    plot(xp,yp,style)

xcenter=0
ycenter=0
RAoffset=array([0,-15.,0,15,-7.5,0,0,0])
DECoffset=array([0,15.,-15.,15,-7.5,15,-15,15])
dither_sequence=array([0,0,0,0,0,1,2,3])
dither_step=30. # in arcsec

RAstep=array([-30, 30,  0, -30, -30,  0, 30, 30,  -30],'f')
DECstep=array([-30 ,0,  30,  0,   0,  30, 0,  0,  -30],'f')

def plotdithers():
    figure()
    clf()
    axis([-1000,1000,-1000,1000])
    for i in range(len(RAoffset)):
        for j in range(len(RAstep)):
            drawNewfirm(RAoffset[i]+RAstep[j],DECoffset[i]+DECstep[j])
            title('test')

def drawNewfirm(xcenter,ycenter):
    # draw four arrays given center position of the instrument (which falls in the gap)
    # rectangle needs (x,y) of lower left corner, height, width
    width=819. # length of individual array in arcsec
    gap = 35. # gap b/w detectors, in arcsec
    # draw first array
    x_lower_left=array([xcenter-width-gap/2.,xcenter-width-gap/2., xcenter+gap/2.,  xcenter+gap/2.])
    y_lower_left=array([ycenter+gap/2.,ycenter-width-gap/2., ycenter+gap/2.,  ycenter-width-gap/2.])
    for i in range(len(x_lower_left)):
        #print x_lower_left[i],y_lower_left[i]
        rect=Rectangle((x_lower_left[i],y_lower_left[i]),width,width,color='0.9',alpha=0.1)
    
        gca().add_patch(rect)
        show()


def makeds9RegionFile(xcenter,ycenter,outfile, color='cyan'):
    # xcenter, ycenter in degrees
    # outfile is the name of the regions file
    out1=open(outfile,'w')
    s='global color=green width=2 \n'
    out1.write(s)
    out1.write('fk5 \n')
    # draw four arrays given center position of the instrument (which falls in the gap)
    # rectangle needs (x,y) of lower left corner, height, width
    width=822. # length of individual array in arcsec
    gap = 35. # gap b/w detectors, in arcsec
    # draw first array

    # convert from arcsec to deg
    width=width/3600.
    gap=gap/3600.

    # calculate positions of each field center that will project
    # correctly in ds9.  Ds9 puts the cos(dec) term back in so I have to
    # take it out here.
    RAwidth=width/cos(ycenter/180.*pi)
    RAgap=gap/cos(ycenter/180.*pi)
    print xcenter,ycenter

    xcenter_box=array([xcenter-RAwidth/2.-RAgap/2.,xcenter-RAwidth/2.-RAgap/2., xcenter+RAwidth/2.+RAgap/2.,  xcenter+RAwidth/2.+RAgap/2.])
    ycenter_box=array([ycenter+width/2.+gap/2.,ycenter-width/2.-gap/2., ycenter+width/2.+gap/2.,  ycenter-width/2.-gap/2.])

    for i in range(len(xcenter_box)):
        #print x_lower_left[i],y_lower_left[i]
        s='box(%12.8f,%12.8f, %12.8f, %12.8f, 0) # color=%s \n'%(xcenter_box[i],ycenter_box[i],width,width,color )
        out1.write(s)

    out1.close()

class cluster:
    def __init__(self,prefix):
        self.prefix=prefix
        self.RAdeg=clusterRAdeg[self.prefix]
        self.Decdeg=clusterDecdeg[self.prefix]

        self.RAoffset=newfirm_offset[self.prefix][0] # in arcsec
        print 'RAoffset = ',self.RAoffset
        self.Decoffset=newfirm_offset[self.prefix][1] # in arcsec


        self.getnewfirmds9reg()
    def getnewfirmds9reg(self):

        RAoffset=self.RAoffset/3600./cos(self.Decdeg/180.*pi) # south 2.5 arcmin
        Decoffset=self.Decoffset/3600. # west 2.5 arcmin
        ds9file='/Users/rfinn/research/z08clusters/RegionsFiles/'+self.prefix+'_newfirm.reg'
        makeds9RegionFile(self.RAdeg,self.Decdeg,ds9file)
        ds9file='/Users/rfinn/research/z08clusters/RegionsFiles/'+self.prefix+'_newfirm_offset.reg'
        makeds9RegionFile(self.RAdeg+RAoffset,self.Decdeg+Decoffset,ds9file, color='blue')



rdcsj1317=cluster('RDCSJ1317')
warp1350=cluster('WARP1350')
rxj1716=cluster('RXJ1716')
rxj1821=cluster('RXJ1821')
