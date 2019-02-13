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
#    Ni  = N3D_class.N3D('BS_1d_20100101_20171231_ptrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc',YML)
    Nb  = N3D_class.N3D('BS_1d_20100101_20171231_btrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc',YML)

#    Ni.testvar('totN',k='bottom')
    Nb.testvar('sCSED')
    Nb.testvar('fCSED')
    Nb.testvar('NCrSED')

    if mm==1:
        N=Nb
    else:
        N.dates        = ma.append(N.dates       , Nb.dates,0)
        N.time         = ma.append(N.time        , Nb.time,0)
#        N.totNkbottom  = ma.append(N.totNkbottom, Ni.totNkbottom,0)
        N.sCSED        = ma.append(N.sCSED       , Nb.sCSED,0)
        N.fCSED        = ma.append(N.fCSED       , Nb.fCSED,0)
        N.NCrSED       = ma.append(N.NCrSED      , Nb.NCrSED,0)
    del Nb
 
#N.makeclim('NOSbottom', climspan=2)
#N.mapStrip('totNkbottom',figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=10)
N.mapStrip('sCSED', figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=15)
N.mapStrip('fCSED', figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=15)
N.mapStrip('NCrSED',figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=15)

#N.mapStrip('totNkbottom',figsuffix='SHELF',cmapname='haline', subdomain="NWS",daysbetween=10)                                                                                                                    
N.mapStrip('sCSED', figsuffix='SHELF_D',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True)
N.mapStrip('sCSED', figsuffix='SHELF_D10',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True, difflag=10)
N.mapStrip('fCSED', figsuffix='SHELF_D',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True)
N.mapStrip('fCSED', figsuffix='SHELF_D10',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True, difflag=10)
N.mapStrip('NCrSED',figsuffix='SHELF_D',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True)
N.mapStrip('NCrSED',figsuffix='SHELF_D10',cmapname='haline', subdomain="NWS",daysbetween=15,diff=True,difflag=10)


#N.mapStrip('totNkbottom',figsuffix='WHOLE',cmapname='haline',daysbetween=10)
N.mapStrip('sCSED', figsuffix='WHOLE',cmapname='haline',daysbetween=15)
N.mapStrip('fCSED', figsuffix='WHOLE',cmapname='haline',daysbetween=15)
N.mapStrip('NCrSED',figsuffix='WHOLE',cmapname='haline',daysbetween=15)

N.mapStrip('sCSED', figsuffix='WHOLE_D' ,cmapname='balance' ,daysbetween=15,diff=True)
N.mapStrip('fCSED', figsuffix='WHOLE_D' ,cmapname='balance' ,daysbetween=15,diff=True)
N.mapStrip('NCrSED',figsuffix='WHOLE_D' ,cmapname='balance' ,daysbetween=15,diff=True)


N.mapStrip('sCSED', figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=15)
N.mapStrip('fCSED', figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=15)
N.mapStrip('NCrSED',figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=15,Clim=[0.05,0.2])

#N.mapStrip('totNkbottom',figsuffix='BOSP',cmapname='haline', subdomain="BOSP",daysbetween=10)                                                                                                                    
N.mapStrip('sCSED', figsuffix='BOSP_D',cmapname='haline', subdomain="BOSP",daysbetween=15,diff=True)
N.mapStrip('fCSED', figsuffix='BOSP_D',cmapname='haline', subdomain="BOSP",daysbetween=15,diff=True)
N.mapStrip('NCrSED',figsuffix='BOSP_D',cmapname='haline', subdomain="BOSP",daysbetween=15,diff=True,Clim=[-.05,.05])

