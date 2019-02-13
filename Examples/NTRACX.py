# We only import librairies needed for plotting
# Other librairies are imported in the class definition file, G3D_class.py,
# which contains all process and variables function definition.

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

import N3D_class
import G3D_class
import numpy.ma as ma
# We instantiate an object of the class G3D, just by giving the path to the netcdf file to work with
# Up to now I'm working with 4D netcdf files containing several variables. 
# Outputs from different files can be merged easily, as can be seen in other examples

for mm in range(1,13):
    Ni  = N3D_class.N3D('BS_1d_20100101_20171231_ptrc_T_2010'+format(mm,'02')+'-2010'+format(mm,'02')+'.nc','local_NEMO_004.yml')

    NmaskDS= (Ni.bat<120 ) & ~(Ni.bat.mask)  # Mask should be True where masked
    Ni.testvar('totN')
    Ni.vitN  = Ni.vertint('totN')
    Ni.aptN  = Ni.avgprofileSIGMA(varname='totN',maskin=NmaskDS)
    Ni.viNOS = Ni.vertint('NOS')
    Ni.apNOS = Ni.avgprofileSIGMA(varname='NOS',maskin=NmaskDS)
    Ni.viNHS = Ni.vertint('NHS')
    Ni.apNHS = Ni.avgprofileSIGMA(varname='NHS',maskin=NmaskDS)
    Ni.viDEN = Ni.vertint('DENITRIFICATION')
    Ni.apDEN = Ni.avgprofileSIGMA(varname='DENITRIFICATION',maskin=NmaskDS)
    Ni.viAMA = Ni.vertint('AMANOX')
    Ni.apAMA = Ni.avgprofileSIGMA(varname='AMANOX',maskin=NmaskDS)
    Ni.viOXN = Ni.vertint('OXIDATIONBYNOS')
    Ni.apOXN = Ni.avgprofileSIGMA(varname='OXIDATIONBYNOS',maskin=NmaskDS)

    if mm==1:
        N=Ni
    else:
        N.dates      = ma.append(N.dates   , Ni.dates,0)
        N.time       = ma.append(N.time    , Ni.time,0)
        N.aptN       = ma.append(N.aptN    , Ni.aptN,0)
        N.apNOS      = ma.append(N.apNOS    , Ni.apNOS,0)
        N.apNHS      = ma.append(N.apNHS    , Ni.apNHS,0)
        N.apDEN      = ma.append(N.apDEN    , Ni.apDEN,0)
        N.apAMA      = ma.append(N.apAMA    , Ni.apAMA,0)
        N.apOXN      = ma.append(N.apOXN    , Ni.apOXN,0)
        N.vitN       = ma.append(N.vitN     , Ni.vitN,0)
        N.viNOS      = ma.append(N.viNOS     , Ni.viNOS,0)
        N.viNHS      = ma.append(N.viNHS     , Ni.viNHS,0)
        N.viAMA      = ma.append(N.viAMA     , Ni.viAMA,0)
        N.viDEN      = ma.append(N.viDEN     , Ni.viDEN,0)
        N.viOXN      = ma.append(N.viOXN     , Ni.viOXN,0)
        

    del Ni


for v in ['tN','NOS','NHS','DEN','AMA','OXN']:
    N.plotprofile('ap'+v,z=-N.z[0,:,0,0])
    N.makeclim('vi'+v)
    N.mapMonthlyClim('vi'+v,figsuffix='SHELF',cmapname='haline', subdomain="NWS")
    N.mapMonthlyClim('vi'+v,figsuffix='WHOLE',cmapname='haline')
    N.mapStrip('vi'+v,figsuffix='_D',daysbetween=30,diff=True)

