# GHER model postprocessing toolbox
A. Capet - Sept 2017

Main class is [G3D_class.py](G3D_class.py).
It contains the G3D class, and a serie of 

* variable definition : for variables derived from model outputs
* processing definition : integrate, average, along various dimensions, using masks, etc ... 
* utilities : eg. ad derived variables to the model ouput files, etc ... 

Examples : 

* [CentralVertProperties.py](CentralVertProperties.py)

Post if you want a working examples with real data.


## Currently defined process function : 

- **.integratespatial(varname)** : Full Spatial Integration 

   -> Returns 1D (time) 
- **.avgspatial(varname)**       : Full Spatial Average     

   -> Returns 1D (time) 
- **.avgprofile(varname[,ztab])**: Horizontally-averaged profile in z coordinate, output z-levels may be specified in argument

   -> Returns 2D (time,zlevels) 
- **.avgprofileSIGMA(varname)**  : Horizontally-averaged profile in sigma coordinate

   -> Returns 2D (time,sigma) 
- **.profileatxy(varname,c1,c2)**: Vert profile at coord C1,C2 (may be indexes (if integers) or lon lat (if floats) 

   -> Returns 2D (time,zlevels) 
   
- **.vertint(varname[,zsup,zinf])**: Vertical integration

   -> Returns 2D (lon,lat) 
