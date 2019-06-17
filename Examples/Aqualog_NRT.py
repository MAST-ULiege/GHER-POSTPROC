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
mlistB = [ 'BS_1d_'+str(y)+'0101_'+str(y)+'1231_ptrc_T_'+str(y)+format(m,'02')+'-'+str(y)+format(m,'02')+'.nc' for y in range(firstyear,lastyear+1) for m in range(1,13)]
mlistP = [ 'BS_1d_'+str(y)+'0101_'+str(y)+'1231_grid_T_'+str(y)+format(m,'02')+'-'+str(y)+format(m,'02')+'.nc' for y in range(firstyear,lastyear+1) for m in range(1,13)]
mlistU = [ 'BS_1d_'+str(y)+'0101_'+str(y)+'1231_grid_U_'+str(y)+format(m,'02')+'-'+str(y)+format(m,'02')+'.nc' for y in range(firstyear,lastyear+1) for m in range(1,13)]
mlistV = [ 'BS_1d_'+str(y)+'0101_'+str(y)+'1231_grid_V_'+str(y)+format(m,'02')+'-'+str(y)+format(m,'02')+'.nc' for y in range(firstyear,lastyear+1) for m in range(1,13)]

vvarP = ['SAL','TEM']
vvarU  = ['U']
vvarV  = ['V'] 
vvarB = ['DOX','CHL','NPP','PAR','POC','NOS','NHS','PHO']

c1 = 37.9144666667
c2 = 44.5395333333

for mi,mm in enumerate(mlistB):

    G1B  = N3D_class.N3D(mm,'local_NEMO_NRT.yml')
    for v in vvarB:
        exec('l'+v+',zp = G1B.profileatxy(v, c1=c1, c2=c2)')
    del G1B

    G1U  = N3D_class.N3D(mlistU[mi],'local_NEMO_NRT.yml')
    lU,zp = G1U.profileatxy(varname = 'vozocrtx', c1=c1, c2=c2)    
    del G1U

    G1V  = N3D_class.N3D(mlistV[mi],'local_NEMO_NRT.yml')
    lV,zp = G1V.profileatxy(varname = 'vomecrty', c1=c1, c2=c2)
    del G1V

    G1P  = N3D_class.N3D(mlistP[mi],'local_NEMO_NRT.yml') 
    for v in vvarP:
        exec('l'+v+',zp = G1P.profileatxy(v , c1=c1, c2=c2)')
        print('l'+v+',zp = G1P.profileatxy(v, c1=c1, c2=c2)')
    
    if mi==0:
        G=G1P
        for v in (vvarP+vvarB+vvarU+vvarV):
            exec('G.l'+v+'=l'+v)
            print('G.l'+v+'=l'+v)

    else:
        G.dates     = ma.append(G.dates, G1P.dates   ,0)
        for v in (vvarP+vvarB+vvarU+vvarV):
            exec('G.l'+v+'= ma.append(G.l'+v+', l'+v+',0)')
    del G1P


for v in (vvarP+vvarB+vvarU+vvarV):
    G.plotprofile(varname='l'+v,z=-zp, zlim=[-300,0])

G.diagfile = 'Aqualogfrom3D.nc'
for v in (vvarP+vvarB+vvarU+vvarV):
    G.gstore('l'+v,ztab=zp)

