# We only import librairies needed for plotting
# Other librairies are imported in the class definition file, G3D_class.py,
# which contains all process and variables function definition.

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import numpy.ma as ma
import N3D_class
import G3D_class

# We instantiate an object of the class G3D, just by giving the path to the netcdf file to work with
# Up to now I'm working with 4D netcdf files containing several variables. 
# Outputs from different files can be merged easily, as can be seen in other examples
YML= 'local_NEMO_004.yml'

for mm in range(1,13):
    Ni  = N3D_class.N3D('BS_1d_20100101_20171231_ptrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc',YML)
#    Nb  = N3D_class.N3D('BS_1d_20100101_20171231_btrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc',YML)

    Ni.testvar('NOS',k='surface')
    Ni.testvar('SIO',k='surface')
    Ni.testvar('PHO',k='surface')
    Ni.viNPP=Ni.vertint('NPP')

    if mm==1:
        N=Ni
    else:
        N.dates        = ma.append(N.dates       , Ni.dates,0)
        N.time         = ma.append(N.time        , Ni.time,0)

        N.NOSksurface  = ma.append(N.NOSksurface       , Ni.NOSksurface,0)
        N.SIOksurface  = ma.append(N.SIOksurface       , Ni.SIOksurface,0)
        N.PHOksurface  = ma.append(N.PHOksurface       , Ni.PHOksurface,0)
        N.viNPP        = ma.append(N.viNPP             , Ni.viNPP,0)
    del Ni

db=31
 
N.mapStrip('NOSksurface', figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=db)
N.mapStrip('PHOksurface', figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=db)
N.mapStrip('SIOksurface', figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=db)
N.mapStrip('viNPP'      , figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=db)

N.mapStrip('NOSksurface', figsuffix='WHOLE',cmapname='haline', subdomain=None,daysbetween=db)
N.mapStrip('NOSksurface', figsuffix='WHOLEb',cmapname='haline', subdomain=None,daysbetween=db, Clim=[0,6])
N.mapStrip('PHOksurface', figsuffix='WHOLE',cmapname='haline', subdomain=None,daysbetween=db)
N.mapStrip('SIOksurface', figsuffix='WHOLE',cmapname='haline', subdomain=None,daysbetween=db)
N.mapStrip('viNPP'      , figsuffix='WHOLE',cmapname='haline', subdomain=None,daysbetween=db)

N.mapStrip('NOSksurface', figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=db)
N.mapStrip('PHOksurface', figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=db)
N.mapStrip('SIOksurface', figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=db)
N.mapStrip('viNPP'      , figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=db)

