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
import cmocean 
import G3D_class
import calendar
import N3D_class

firstyear=2016
lastyear =2016

# Compile
mlistP = [ 'BS_1d_'+str(y)+'0101_'+str(y)+'1231_grid_T_'+str(y)+format(m,'02')+'-'+str(y)+format(m,'02')+'.nc' for y in range(firstyear,lastyear+1) for m in range(1,13)]

for mi,mm in enumerate(mlistP):

    G1P  = N3D_class.N3D(mlistP[mi],'local_NEMO_NRT.yml') 
#    G1P.diagfile =  './diags_NRT/'+mm
    ttt=G1P.testvar('DEN', doload=False)
    if not ttt:
        G1P.testvar('TEM')
        G1P.TTT=G1P.TEM
        G1P.gstore('TTT')

        ttt=G1P.testvar('DEN')
        G1P.gstore('DEN')

    ttt=G1P.testvar('Z14_5', doload=False)
    if not ttt:
        ttt=G1P.testvar('Z14_5')
        G1P.gstore('Z14_5')
    
    if mi==0:
        G=G1P
    else:
        G.dates     = ma.append(G.dates, G1P.dates   ,0)
        for v in (['Z14_5']):
           exec('G.'+v+'= ma.append(G.'+v+', G1P.'+v+',0)')
    del G1P


G.mapStrip('Z14_5',daysbetween=8)
