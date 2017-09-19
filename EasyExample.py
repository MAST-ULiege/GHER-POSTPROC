# We only import librairies needed for plotting
# Other librairies are imported in the class definition file, G3D_class.py,
# which contains all process and variables function definition.

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

import G3D_class

# We instantiate an object of the class G3D, just by giving the path to the netcdf file to work with
# Up to now I'm working with 4D netcdf files containing several variables. 
# Outputs from different files can be merged easily, as can be seen in other examples

G  = G3D_class.G3D('../Out_001/1982.nc')

# All loaded variables are attributes of the G3D instance.
# For instance the variable "bat" is defined directly when the object is instantiated. 
# Other are loaded only when needed.
# Variables are python Masked_array, so they have an attribute mask which is an arreay of booleans 

# Here we want to define a mask based on bathymetry 
maskDS= (G.bat<50 ) & ~(G.bat.mask)  # Mask should be True where masked
maskSH= (G.bat>=50) & ~(G.bat.mask)  # Mask should be True where masked                                                                                                                                         

# All processing functions are called as function of the G3D instance. 
# Variable name is given as an argument. Some functions allows more argument.

# This would give the basin averaged time series of salinity
T1 = G.avgspatial('SAL')

# The avgspatial function enables an optional mask argument
# Note also , that we can use a variable name that is not defined in the netcdf file. 
# In this case the toolbox will automatically look for the function "instance_SSS"

sssDS=G.avgspatial('SSS',maskDS)
sssSH=G.avgspatial('SSS',maskSH)

# The following is general python plotting .. 
# the "dates" attributes is also loaded automatically 

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               

locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 2)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(G.dates,sssSH, label = 'Average surface salinity on the shelf')
plt.plot(G.dates,sssDS, label = 'Average surface salinity in the open sea')
plt.title('Salinity')
plt.ylabel('Salinity - [p.s.u.]')
fig.savefig(G.figoutputdir+'Simple.png')

