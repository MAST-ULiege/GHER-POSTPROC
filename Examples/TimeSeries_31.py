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

ymlfile1='local_NEMO_31.yml'

plist = ['SST','MLD']
blist = ['VOX','Z20','NPPint'] 

vlist = plist+blist

for gi,ymlfile in enumerate([ymlfile1]):
    print (ymlfile)
    mlist = N3D_class.FullLoad(ymlfile,begdate=dt.datetime(2019,1,1), enddate=dt.datetime(2019,12,31))

#    mlist=mlist[:-6]
    for mm in mlist:
        G1  = N3D_class.N3D(mm.replace('ptrc','grid'),ymlfile)

        maskDS = (G1.bat>60) | (G1.bat.mask)
#        reg = G1.loadregionmask('./regions.nc','kopelevich')
#        G1.m5 = (reg != 5.0) & ~ (G1.bat.mask)

        for v in plist:
            G1.testvar(v)
            exec('G1.'+v+"_mean=G1.avgspatial(v,maskin=maskDS)")
 
        G1b = N3D_class.N3D(mm,ymlfile)

        for v in blist:
             G1b.testvar(v)
             exec('G1b.'+v+"_mean=G1b.avgspatial(v,maskin=maskDS)")

        if mm==mlist[0]:
            G1a=G1
            for v in blist:
                exec('G1a.'+v+"_mean=G1b."+v+"_mean")                        
        else:
            G1a.dates = ma.append(G1a.dates,    G1.dates   ,0)
            for v in plist:
                exec('G1a.'+v+'_mean=ma.append(G1a.'+v+'_mean,G1.'+v+'_mean)') 
            for v in blist:
                exec('G1a.'+v+'_mean=ma.append(G1a.'+v+'_mean,G1b.'+v+'_mean)')
        del G1,G1b

        G1a.timeclean()

    for v in vlist:
        G1a.gstore(v+'_mean')
        G1a.plotseries(v+'_mean')


'''
        exec('GG'+str(gi+1)+'=G1a')

for v in vlist: 
    G1a.plotseriesMultipleRuns([v+'Bflu_r5'],G3Ds=(GG1,GG2,GG3), Glabs = ['try 1','try 2','try 3'],figout='comp_3TRY_r5_Bflu'+v,title='Benthic Flux '+v+' - region5')
    G1a.plotseriesMultipleRuns([v+'Surf_r5'],G3Ds=(GG1,GG2,GG3), Glabs = ['try 1','try 2','try 3'],figout='comp_3TRY_r5_Surf'+v,title='Surface Conc '+v+' - region5')
    G1a.plotseriesMultipleRuns([v+'Bott_r5'],G3Ds=(GG1,GG2,GG3), Glabs = ['try 1','try 2','try 3'],figout='comp_3TRY_r5_Bott'+v,title='Bottom Conc  '+v+' - region5')
'''
