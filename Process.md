## .integratespatial(varname)

Full Spatial Integration 

-> Returns 1D (time) 

## .avgspatial(varname)

Full Spatial Average     

   -> Returns 1D (time) 
   
## .avgprofile(varname[,ztab])

Horizontally-averaged profile in z coordinate, output z-levels may be specified in argument

   -> Returns 2D (time,zlevels) 
   
## .avgprofileSIGMA(varname)

Horizontally-averaged profile in sigma coordinate

   -> Returns 2D (time,sigma) 
   
## .profileatxy(varname,c1,c2)

Vert profile at coord C1,C2 (may be indexes (if integers) or lon lat (if floats) 

   -> Returns 2D (time,zlevels) 
   
## .vertint(varname[,zsup,zinf])

Vertical integration

   -> Returns 2D (lon,lat) 
