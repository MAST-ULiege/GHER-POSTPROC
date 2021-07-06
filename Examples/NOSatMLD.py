import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
                                                                                                                                                                                                                   
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
import pandas as pd

ymlfile='./local_NEMO_MYP.yml'

mlist = N3D_class.FullLoad(YAML_FILE =ymlfile,begdate=dt.datetime(2014,1,1), enddate=dt.datetime(2014,12,31))

flist2D = ['NOSatMLD']

for mm in mlist[:]:
    # This one use the lists generated for btrc_T file, but substitutes parts of the filename to reach somld_bs
    Pi  = N3D_class.N3D(mm.replace('ptrc','grid'),ymlfile)
    Ni  = N3D_class.N3D(mm,ymlfile)
    # A mask for a spatial average over the Deep Sea region. 
    # Mask should be True where masked
    maskDS= (Ni.bat.mask) | (Ni.bat<120) 
    # This loads somld_bs.. it's from the isntance constructed wiuth the Physics file. 
    Pi.testvar('somld_bs')
    # we jsut copy the filed to the btrc_T instance
    Ni.somld_bs=Pi.somld_bs
    # ok the loop on v isn't needed here, just inherited from previous script, and still usefull to maintain the structure for later
    for v in flist2D:
        # loads 'NOSatMLD',  using the dedicated function in G3D.class . The first thing this function does is checking if somld_bs is available. 
        Ni.testvar(v)
        setattr(Ni,v+'DS',Ni.avgspatial(v,maskin=maskDS))


    # The following builds up a sacked version of NOSatMLD along loop iterations. 
    # Note that we only stack NOSatMLD (and dates) because we don't need MLD anymore. 
    if mm==mlist[0]:
        N=Ni
    else:
        N.dates      = ma.append(N.dates   , Ni.dates,0)
        N.time       = ma.append(N.time    , Ni.time,0)
        for v in flist2D:
            setattr(N,v,ma.append(getattr(N,v), getattr(Ni,v),0))
            setattr(N,v+'DS',ma.append(getattr(N,v+'DS'), getattr(Ni,v+'DS'),0))

    del Ni

# This ensure there is no repeated time instance. 
N.timeclean()

for v in flist2D:
    # This generate monthly Maps
    N.SeasonStrip(v, extend='both')
    # this generate a netcdf with the desired diagnostics. 
    N.gstore(v)

