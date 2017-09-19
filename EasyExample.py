import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

import G3D_class

G  = G3D_class.G3D('../Out_001/1982.nc')
T1 = G.avgspatial('SAL')

locator = mdates.AutoDateLocator()
formator = mdates.AutoDateFormatter(locator)

####################                                                                                                                                                                                               
# 1st figure :  
####################                                                                                                                                                                                               

fig=plt.figure(figsize=(15, 8))

ax=plt.subplot(1, 1, 1)
ax.xaxis_date()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formator)
plt.plot(G.dates,T1)
plt.title('Average Salinity')
plt.ylabel('Salinity - [p.s.u.]')
fig.savefig(G.figoutputdir+'Simple.png')
