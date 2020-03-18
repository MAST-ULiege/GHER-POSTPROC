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

## ## ## ## ## ###  ## ## ## ## ## ###  ## ## ## ## ## ### 
## Arguments Management ##

import  getopt,sys
 
arglist = sys.argv[1:]
unixOptions = "y:"
gnuOptions = ["yml="]
try:
    arguments, values = getopt.getopt(arglist, unixOptions, gnuOptions)
except getopt.error as err:
    # output error, and return with an error code                                                                                                                                                                   
    print (str(err))
    sys.exit(2)

ymlfile   = None

for a,v  in arguments:
        if a in ("-y",'--yml'):
            ymlfile = v

if ymlfile is None:
    print( 'Dude, I need some -y or --yml argument for the .yml file')
    sys.exit()

## ## ## ## ## ##### ## ## ## ## ###  ## ## ## ## ## ###         
 
mlist = N3D_class.FullLoad(ymlfile)

vlist = ['botfluxDOX','botfluxNOS','botfluxNHS','botfluxDIC','fCSED','sCSED','NCrSED','fSSED','sSSED']#, 'pnit','pdenit', 'panox']

for mm in mlist:
    G1  = N3D_class.N3D(mm.replace('ptrc','btrc'),ymlfile)
    
    reg = G1.loadregionmask('./regions.nc','kopelevich')

    G1.m4 = (reg != 4.0) & ~ (G1.bat.mask)
    G1.m5 = (reg != 5.0) & ~ (G1.bat.mask)

    print(G1.m4.shape)

    for v in vlist:
        G1.testvar(v)#,k='surface')
        exec('G1.'+v+"_r4=G1.avgspatial(v,maskin=G1.m4)")
        exec('G1.'+v+"_r5=G1.avgspatial(v,maskin=G1.m5)")

    if mm==mlist[0]:
        G1a=G1
    else:
        G1a.dates = ma.append(G1a.dates,    G1.dates   ,0)
        for v in vlist:
            exec('G1a.'+v+'_r4=ma.append(G1a.'+v+'_r4,G1.'+v+'_r4,0)') 
            exec('G1a.'+v+'_r5=ma.append(G1a.'+v+'_r5,G1.'+v+'_r5,0)') 
    del G1

G1a.timeclean()

for v in vlist:
    if (v in ['botfluxDOX','botfluxNOS','botfluxNHS','botfluxDIC']):
            exec('G1a.'+v+'_r4=G1a.'+v+'_r4*86400')
            exec('G1a.'+v+'_r5=G1a.'+v+'_r5*86400')

for v in vlist: 
    G1a.plotseries([v+'_r5',v+'_r4'] ,figout='compRB_'+v, title=v)
