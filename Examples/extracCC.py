import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap  
#from multiprocessing import Pool
#import gsw                                                                                                                                                                                                         
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import sys
import os

import G3D_class

G1a     = G3D_class.G3D('../Out/PACK4.nc')
#G1b     = G3D_class.G3D('../Out/PACK2.nc')
#G1c     = G3D_class.G3D('../Out/PACK3.nc')


G1a.dates = [dt.datetime(2006,11,12)+dt.timedelta(days=int(t*7)) for t in range(len(G1a.time))]
G1a.time  = [t for t in range(len(G1a.time))]

maskDS= (G1a.bat<120) | (G1a.bat.mask)

G1a.cc  = G1a.avgspatial('CCCn', maskin=maskDS.squeeze())

G1a.gstore('cc')
