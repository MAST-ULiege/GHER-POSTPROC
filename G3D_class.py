import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt

class G3D: 
    '''let's test python classes'''

######################################################################

    def __init__(self,infile):
        self.infile=infile
        if os.path.isfile(self.infile):
            print(self.infile + ' -> OK')
        else:
            print(self.infile + ' can not be found')
       
        try:
            YAML_FILE = 'local.yml'
        except Exception:
            print("".join(("\n A file called local.yml should be present",
                       "'\n")))
        print("\nLaunching with YAML file: %s" % YAML_FILE)


    # Read yaml configuration file
        with open(YAML_FILE, 'r') as stream:
            config = yaml.load(stream)
        
    # this should be read somewhere
        self.batfile      = config['BATFILE']
        self.dx           = config['DX']
        self.dy           = config['DY']
    # to get intervals in the double sigma context of GHER
        self.sigI         = config['SIGI']
        self.sigII        = config['SIGII']
        self.hlim         = config['HLIM']
        self.verbose      = config['VERBOSE']
        self.figoutputdir = config['PLOTDIR']
        self.instance_bat()

######################################################################
# VARIABLE DEF : BAT

    def instance_bat(self):
        if os.path.isfile(self.batfile):
            print(self.batfile + ' -> OK')
        else:
            print(self.batfile + ' can not be found')

        with Dataset(self.batfile,'r') as nc:
            self.bat = nc.variables['bat'][:]
            
######################################################################
# VARIABLE DEF : Z

    def instance_z(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        with Dataset(self.infile,'r') as nc:
            self.z = -1* nc.variables['depth'][1] # Let it be 3D for now           
            self.lon = nc.variables['longitude'][:]
            self.lat = nc.variables['latitude'][:]

######################################################################
# VARIABLE DEF : ZI
            
    def instance_zi(self):
        # For now let's build a constant dz,
        # dynamic dz considering ETA can be done later
        self.zi=ma.zeros(self.z.shape+np.array([1,0,0]))
        for k in xrange(0,12):
            self.zi[k]= ma.masked_where( self.bat<=self.hlim,  (( self.bat - self.hlim ) * (1-self.sigII[k]) ) + self.hlim )
        for k in xrange(12,self.zi.shape[0]):
            self.zi[k]= np.minimum(self.bat,self.hlim) * (1-self.sigI[k-12+1])
            
######################################################################
# VARIABLE DEF : DZ

    def instance_dz(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        self.dz = ma.abs(self.zi[1:]-self.zi[:-1]) 

######################################################################
# VARIABLE DEF : DENSITY

#    def instance_density(self):
    
#        p=gsw.p_from_z(,ttlat)
#    SA = gsw.SA_from_SP(S[t,:,:,:],p,ttlon,ttlat)
#    # I considered in situ temperature and ractical salinity from the model .. absolute salinity/practical salinity Temperature/Potential temperature ?
#    Dloc=gsw.rho(SA,T[t,:,:,:],p)
#        # For now let's build a constant z,
#        # dynamic z considering ETA can be done later
#        self.dz = ma.abs(self.zi[1:]-self.zi[:-1])

######################################################################
# UTILITY : LOAD
        
    def gload(self,varname):
        with Dataset(self.infile,'r') as nc:
            try:
                exec('self.'+varname+ '= nc.variables[''varname''][:]')
                print( 'Just loaded %s'%(varname))
            except: 
                print( '%s not found in %s'%(varname,self.infile))
                print( '-> try to compute')
             #   try:
                exec('self.instance_'+varname+'()')
               # except:
               #     print('self.instance_'+varname+' is not defined.')
                    
                
######################################################################
# UTILITY : STORE
#        
#    def gstore(self,varname):
#        with Dataset(self.infile,'a') as nc:
#            try:
#                nc.createVariable(varname, np.float32, ('time','level', 'latitude', 'longitude'),zlib=True)
#            exec('nc.variables[varname][:]=self.'+varname)
#                
######################################################################
# PROCESS : FULL INTEGRATION

    def integratespatial (self,varname):
    # TODO render dz optionnal, allows to use contant depth, and dz instead
    # TODO allows for a optional mask, that could be 2D, 3D or 4D. 
    #    if dz.shape!=field.shape:
    #        print("Wrong Shapes")
    #        return
    #    if len(field.shape)==3:
    #        print("3D field")
    #    
        self.testz()
        self.testvar(varname)
        self.testtime()

        integrated=ma.empty(self.time.shape)
        exec('loc=self.'+varname)
        print('dz: %s  and field: %s'%( len(self.dz.shape),len(loc.shape)))
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):
            print("4D")
            for t in xrange(self.time.shape[0]):
                bi=loc[t]*self.dz*self.dy*self.dx
                integrated[t] = ma.sum(bi)
        # return should be 1D (time)    
        return(integrated)
        
######################################################################
# PROCESS : FULL AVERAGE

    def avgspatial (self,varname):
        # in the sense of volumetric mean
        self.testz()
        self.testvar(varname)        
                        
        avg=ma.empty(self.time.shape)
        exec('loc=self.'+varname)
        print('dz: %s  and field: %s'%( len(self.dz.shape),len(loc.shape)))
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):
            print("4D")
            for t in xrange(self.time.shape[0]):
                bi=loc[t]*self.dz*self.dy*self.dx
                vol=ma.masked_where(bi.mask,self.dz*self.dy*self.dx)
                avg[t] = ma.sum(bi)/ma.sum(vol)
        else:
            print('this dimensional case is not condidered yet')
                
        # return should be 1D (time)    
        return(avg)
        
######################################################################
        
    def testz(self):
        try:
            self.dz
        except :
            print('dz not found -> loading')
            self.instance_z()
            self.instance_zi()
            self.instance_dz()
            
######################################################################
            
    def testvar(self,varname):
        try:
            exec('self.'+varname)
        except:
            print('%s not found -> loading'%(varname))
            self.gload(varname)
            
######################################################################
            
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time')
            self.dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in self.time]
        
######################################################################
# PROCESS : Horizontally-averaged profile in z coordinate, taking sigma coordinate in consideration

    def avgprofile(self,varname,
            ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)])
             ):
        # The idea is to get an average profile as a function of time
        # return is 2D : [z,time]
        # sigma-space makes it a bit complicate
        #ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)])
        #ztab=-1*np.arange(0,100,5)
        
        self.testz()
        self.testvar(varname)
        self.testtime()
        
        avg = ma.empty((self.time.shape[0],ztab.shape[0]-1))
        exec('loc=self.'+varname)
        
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):            
            gridZU = -self.zi[1:]
            gridZD = -self.zi[:-1]
            for k in xrange(ztab.shape[0]-1):
                print('%s / %s'%(k+1,ztab.shape[0]-1))
                dzloc= ma.maximum(ma.zeros(self.dz.shape), np.minimum(gridZU, ztab[k])-np.maximum(gridZD, ztab[k+1]))
                for t in xrange(self.time.shape[0]):
                    vol=ma.masked_where(loc[t].mask,dzloc*self.dy*self.dx)
                    bi=loc[t]*vol
                    avg[t,k]= ma.sum(bi)/ma.sum(vol)
                    if (t==10):
                        print('k %s : from %s to %s' %(k,ztab[k],ztab[k+1]))
                        print('total volume considered for this layer :  %s km^3' %(ma.sum(vol)/1e9))
                        print('mean age :  %s' %(avg[t,k]))
                    
        zforplot=(ztab[:-1]+ztab[1:])/2
        return avg, zforplot
        

######################################################################
# TEST : spatial coordinates

    def test_coord(self,c1,c2):
        if (type(c1) is int) & (type(c2) is int):
           #considered as indexes                                                                                                                                   
            i=c1
            j=c2
        else:
           #considered coordinates   
            X=c1-self.lon
            i=X.index(np.min(np.abs(X)))
            Y=c2-self.lat
            j=Y.index(np.min(np.abs(X)))
        return i,j
        
######################################################################    
    
    def getlonlat(self,c1,c2):
        i,j=self.test_coord(c1,c2)
        lon=self.lon[i]
        lat=self.lat[j]
        return lon,lat


######################################################################
#
# PROCESS : Depth profile at given coordinate
#
#  inputs : 
#          * varname 
#          * indexes : integer will be interpretated as grid indexes, floats as coordinates (lon,lat) 

    def profileatxy(self,varname,c1,c2):

        self.testz()
        self.testvar(varname)
        self.testtime()
        i,j=self.test_coord(c1,c2)

        zforplot=self.z[:,j,i]
#        zforplot=zforplot[~zforplot.mask]
        exec('loc=self.'+varname)

        print loc.shape
        lloc=loc[:,:,j,i]
        print lloc.shape

        return lloc,zforplot
