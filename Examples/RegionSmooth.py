
import cv2 as cv
import numpy as np
from scipy import signal
import scipy as sp
import scipy.ndimage
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap
#from multiprocessing import Pool
#import gsw                                                                                                                                                                                             \

import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import sys
import os
import cmocean
import G3D_class
import calendar
import N3D_class
import numpy.ma as ma

ymlfile='./local_NEMO_diags.yml'
## ## ## ## ## ##### ## ## ## ## ###  ## ## ## ## ## ###         

mlist = N3D_class.FullLoad(YAML_FILE =ymlfile, begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31))
flist2D = ['BenthicOxygenFlux']

for v in flist2D:
    for mm in mlist[:]:
        Ni  = N3D_class.N3D(mm,ymlfile)
        NmaskDS= ~(Ni.bat.mask)  # Mask should be True where masked                                                                                                                                 
        Ni.testvar(v)
        
        if mm==mlist[0]:
            N=Ni
        else:
            N.dates      = ma.append(N.dates   , Ni.dates,0)
            N.time       = ma.append(N.time    , Ni.time,0)
            setattr(N,v,ma.append(getattr(N,v), getattr(Ni,v),0))

    del Ni

    N.timeclean(begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31))

reg=N.loadregionmask('regio_from_O2diags.nc','regio')



def gkern(kernlen=21, std=3, height=5):
    """Returns a 2D Gaussian kernel array."""
    gkern1d = signal.gaussian(kernlen, std=std).reshape(kernlen, 1)
    gkern2d = np.outer(gkern1d, gkern1d)
    gkern2d = np.round(gkern2d*height)
    return gkern2d

def ApplyKernel(infield, kernel) : 
    '''
    :param infield: input masked array
    :param kernal: input kernel
    '''
    import statistics as stats
    output=infield.copy()
    ks= int(np.floor(kernel.shape[0]/2))

    for i in np.arange(ks+1,infield.shape[0]-ks-1,1):
        for j in np.arange(ks+1,infield.shape[1]-ks-1,1):
            if infield.mask[i,j]:
                continue
            inloc = infield[(i-ks):(i+ks), (j-ks):(j+ks)]
            myvec = ma.empty(int(kernel.sum()))
            cc=0
            for ik in range(2*ks+1):
                for jk in range(2*ks+1):
                    n = int(kernel[ik,jk])
                    for n in range(n):
                        myvec[cc]=inloc[ik,jk]
                        cc+=1
            output[i,j]=stats.mode(myvec.compressed())
        # get surrounding subfield
    return output

kernel=gkern(kernlen=21, std=5,height=3)
oo=ApplyKernel(reg,kernel)







