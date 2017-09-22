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
        print(' *******  \n')
        if os.path.isfile(self.infile):
            print(self.infile + ' -> OK')
            self.found=True
        else:
            print(self.infile + ' can not be found')
            self.found=False
            return
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
        self.testtime()

######################################################################
# VARIABLE : BAT

    def instance_bat(self):
        if os.path.isfile(self.batfile):
            print(self.batfile + ' -> OK')
        else:
            print(self.batfile + ' can not be found')

        with Dataset(self.batfile,'r') as nc:
            self.bat = nc.variables['bat'][:]
            
######################################################################
# VARIABLE : Z

    def instance_z(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        with Dataset(self.infile,'r') as nc:
            self.z   = nc.variables['depth'][1] # Let it be 3D for now           
            self.lon = nc.variables['longitude'][:]
            self.lat = nc.variables['latitude'][:]

######################################################################
# VARIABLE : ZI
            
    def instance_zi(self):
        # For now let's build a constant dz,
        # dynamic dz considering ETA can be done later
        self.zi=ma.zeros(self.z.shape+np.array([1,0,0]))
        for k in xrange(0,12):
            self.zi[k]= - ma.masked_where( self.bat<=self.hlim,  (( self.bat - self.hlim ) * (1-self.sigII[k]) ) + self.hlim )
        for k in xrange(12,self.zi.shape[0]):
            self.zi[k]= - np.minimum(self.bat,self.hlim) * (1-self.sigI[k-12+1])
            
######################################################################
# VARIABLE : DZ

    def instance_dz(self):
        # For now let's build a constant z,
        # dynamic z considering ETA can be done later
        self.dz = ma.abs(self.zi[1:]-self.zi[:-1]) 

######################################################################
# UTILITY : LOAD
        
#    def gload(self,varname):
#        with Dataset(self.infile,'r') as nc:
#            try:
#                exec('self.'+varname+ '= nc.variables[''varname''][:]')
#                print( 'Just loaded %s'%(varname))
#            except: 
#                print( '%s not found in %s'%(varname,self.infile))
#                print( '-> try to compute')
#             #   try:
#                exec('self.instance_'+varname+'()')
#               # except:
#               #     print('self.instance_'+varname+' is not defined.')

######################################################################
# UTILITY : LOAD LOCAL 
        
    def gload(self,varname,ti=None,k=None,i=None,j=None):
        with Dataset(self.infile,'r') as nc:
            try:
                if (ti is None) and (k is None) and (i is None) and (j is None):
                    exec('self.'+varname+ '= nc.variables[varname][:]')
                    print( 'Just loaded '+ (varname) + ' for full')
                elif (ti is None) and (k is None) and (i is not None) and (j is not None):
                    print( 'Loading '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                    exec('self.'+varname+'i'+str(i)+'j'+str(j)+'= nc.variables[varname][:,:,j,i]')
                    print( 'Just loaded '+ (varname) + ' for i:'+str(i)+' and j:'+str(j))
                elif (ti is None) and (k is not None) and (i is None) and (j is None):
                    exec('self.'+varname+'k'+str(k)+'= nc.variables[varname][:,k]')
                    print( 'Just loaded '+ (varname) +' for ki:'+str(k))
                else:
                    print(' Stange case encountered in  G3D_class.py : def gload ')

            except Exception as e: 
                print(e)
                print( '\n %s not found in %s'%(varname,self.infile))
                print( '-> Calling')
                #try:
                print('     self.instance_'+varname+'(i='+str(i)+',j='+str(j)+',k='+str(k)+')')
                exec('self.instance_'+varname+'(i=i,j=j,k=k)')
                #except:
                #    print('self.instance_'+varname+' is not defined.')
                    
                
######################################################################
# UTILITY : STORE
#        
    def gstore(self,varname):
        exec('ndim=len(self.'+varname+'.shape)')
        print('\n Storing now '+varname+' ('+ str(ndim)+' dimensions) on '+self.infile)
        exec('print(self.'+varname+'.shape)')

        with Dataset(self.infile,'a') as nc:
            try:
                if ndim==4:      # assuming here : time, level, lat,lon
                    nc.createVariable(varname, np.float32, ('time','level', 'latitude', 'longitude'),zlib=True)
                elif ndim == 3:  # assuming here : time, lat,lon 
                    nc.createVariable(varname, np.float32, ('time', 'latitude', 'longitude'),zlib=True)
            except:
                print ('Seems that '+varname+' already exists on '+self.infile+'. \n I overwrites')

            exec('nc.variables[varname][:]=self.'+varname)
                
######################################################################
# PROCESS : FULL INTEGRATION

    def integratespatial (self,varname):
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

    def avgspatial (self,varname,maskin=None):
        # in the sense of volumetric mean
        self.testz()
        self.testvar(varname)        
                        
        avg=ma.empty(self.time.shape)
        exec('loc=self.'+varname+'.copy()')
        print('In avgptial \n dz: %s  \n Field: %s'%( len(self.dz.shape),len(loc.shape)))

        # DIMENSIONAL CASE 
        # 2D VARIABLE
        print(loc.shape)
        if ( loc.shape[1]==1 ):
            print('2D variable')
            # MASKING
            if (maskin is not None):
                print('Found Mask with ' + str(maskin.sum()) + ' masked points.')
                if len(maskin.shape)==3:
                    print('2D mask .. OK')
                    for t in xrange(self.time.shape[0]):
                        loc[t]=ma.masked_where(maskin,loc[t])

            # AVERAGING
            for t in xrange(self.time.shape[0]):
                bi=loc[t]*self.dy*self.dx
                vol=ma.masked_where(bi.mask,self.dy*self.dx*np.ones(bi.shape))
                avg[t] = ma.sum(bi)/ma.sum(vol)

        # 3D VARIABLE
        elif (loc.shape[1]>1):
            if (maskin is not None):
                print('Masking of 3D vars not implemented yet .. Complete G3D_class.py') 
            print("3D variable")
            for t in xrange(self.time.shape[0]):
                bi=loc[t]*self.dz*self.dy*self.dx
                vol=ma.masked_where(bi.mask,self.dz*self.dy*self.dx)
                avg[t] = ma.sum(bi)/ma.sum(vol)
        else:
            print('this dimensional case is not condidered yet')
                
        # return should be 1D (time)    
        return(avg)
        
######################################################################
# UTILITY : test z        
    def testz(self):
        try:
            self.dz
        except :
            print('dz not found -> loading')
            self.instance_z()
            self.instance_zi()
            self.instance_dz()
            
######################################################################
# UTILITY : test variable
            
    def testvar(self,varname,doload=True):
        try:
            exec('self.'+varname)
            isthere=True
        except:
            print('%s not found '%(varname))
            isthere=False
            if doload:
                print ('Loading %s'%(varname) )
                self.gload(varname)
                isthere=True

        return(isthere)
                
            
######################################################################
# UTILITY : test time
            
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time')
            self.dates = [dt.datetime(1858,11,17)+dt.timedelta(days=int(t)) for t in self.time]

######################################################################
# UTILITY : Apply a mask on a variable, handling dimension cases. 

    def maskvar(self,loc, maskin):
    
        if (len(loc.shape)==4) and (len(maskin.shape)==2):
            print('Masking '+ str(len(loc.shape))+ 'D variable with '+str(len(maskin.shape))+'D mask')
            for t in xrange(loc.shape[0]):
                for k in  xrange(loc.shape[1]):
                    loc[t,k,:,:]=ma.masked_where(maskin,loc[t,k,:,:]) # ma.expand_dims(ma.expand_dims(maskin,0),0)
        else:
            print('NOT Masking '+ str(len(loc.shape))+ 'D variable with '+str(len(maskin.shape))+'D mask')
                
        return loc

        
######################################################################
# PROCESS : Horizontally-averaged profile in z coordinate, taking sigma coordinate in consideration

    def avgprofile(self,varname,
                   ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)]),
                   maskin=None
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

        if maskin is not None:
            loc = self.maskvar(loc,maskin)
       
        if (len(self.dz.shape)==3)and(len(loc.shape)==4):            
            gridZU = self.zi[1:]
            gridZD = self.zi[:-1]
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
# PROCESS : Horizontally-averaged profile in sigma coordinate

    def avgprofileSIGMA(self,varname, maskin=None):
        self.testvar(varname)
        self.testtime()
        exec('loc=self.'+varname)
        
        if maskin is not None:
            loc = self.maskvar(loc,maskin)

        avg=ma.average(loc,(2,3))
        return avg

######################################################################
# UTILITY : test spatial coordinates

    def test_coord(self,c1,c2):
        if (type(c1) is int) & (type(c2) is int):
           #considered as indexes                                                                                                                                   
            i=c1
            j=c2
        else:
           #considered coordinates   
            X=c1-self.lon
            i=abs(X).argmin()
            Y=c2-self.lat
            j=abs(Y).argmin()
        return i,j
        
######################################################################    
# UTILITY : provide lon lat coordinates for indexes
    
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
        isvar=self.testvar(varname,False)
        self.testtime()
        i,j=self.test_coord(c1,c2)

        zforplot=self.z[:,j,i]
        if isvar:
            # for some reason the 4d variable is already loaded in memory
            exec('loc=self.'+varname)
            lloc=loc[:,:,j,i]
        else:
            # It's not loaded, but we'll load only what we need
            self.gload(varname,i=i,j=j)
            exec('lloc=self.'+varname+'i'+str(i)+'j'+str(j))

        return lloc,zforplot

##########################################################################
# PROCESS : Vertical Integration

    def vertint(self,varname,zinf=-10000,zsup=2):

        self.testz()
        self.testvar(varname)
        self.testtime()

        exec('loc=self.'+varname)
        vint=ma.empty( (loc.shape[0],loc.shape[2],loc.shape[3]) )

        for t in xrange(self.time.shape[0]):
            vint[t]=ma.sum(loc[t]*self.dz,0)
        
        return vint

#########################################################################                                                                                                                                          
# PROCESS : Vertical Mean
  
    def vertmean(self,varname,zinf=-10000,zsup=2):

        self.testz()
        self.testvar(varname)
        self.testtime()

        exec('loc=self.'+varname)
        vmean=ma.empty( (loc.shape[0],loc.shape[2],loc.shape[3]) )
        print(vint.shape)
        for t in xrange(self.time.shape[0]):
            vmean[t]  = ma.sum ( loc[t]*self.dz , 0)
            vol      = ma.sum (        self.dz , 0)
            vmean[t] = vmean[t]/vol

        return vmean

############################################################################
# VARIABLE : Density

    def instance_DEN(self,i=None,j=None,k=None):

        self.testz()
        tlat=np.tile(self.lat, (self.z.shape[0],self.z.shape[2],1))
        tlat=np.transpose(tlat,(0,2,1))
        tlon=np.tile(self.lon, (self.z.shape[0], self.z.shape[1],1))

        if (i is None) and (j is None) and (k is None):
            p  = gsw.p_from_z(self.z,tlat)
            self.testvar('SAL')
            self.testvar('TEM')
            self.DEN = ma.empty(self.SAL.shape)
            for t in xrange(self.SAL.shape[0]):
                SA = gsw.SA_from_SP(self.SAL[t],p,tlon,tlat)
                self.DEN[t]=gsw.rho(SA,self.TEM[t],p)
            
            if False:
                self.gstore('DEN')

            del self.SAL
            del self.TEM

        elif (i is not None) and (j is not None) and (k is None):
            # need to be computed only for one profile
            p  = gsw.p_from_z(self.z[:,j,i],tlat[:,j,i])
#            self.DEN = ma.empty(self.SAL.shape)
            for t in xrange(self.time.shape[0]):
                self.gload('SAL',i=i,j=j)
                self.gload('TEM',i=i,j=j)
                exec('SALloc=self.SALi'+str(i)+'j'+str(j))
                exec('TEMloc=self.TEMi'+str(i)+'j'+str(j))
                SA = gsw.SA_from_SP(SALloc[t],p,tlon,tlat)
                DENloc=gsw.rho(SA,TEMloc[t],p)
                exec('self.DENi'+str(i)+'j'+str(j)+'=DENloc')

            if False:
                self.gstore('DEN')

############################################################################
# UTILITY : unload to free some memory


############################################################################
# VARIABLE SSS

    def instance_SSS(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.gload('SAL',k=30) 
            self.SSS=self.SALk30
            del self.SALk30
            self.SSS=ma.expand_dims(self.SSS,1)
            
        else:
            print(' Case not ready, please code .. ')

        self.SSSinfo.units='p.s.u.'
        self.SSSinfo.longname='Surface Salinity'
        self.SSSinfo.dims='2D'


############################################################################                                                                                                                                       # VARIABLE CCCn

    def instance_CCCn(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.testvar('DEN')
            self.testvar('TEM')
            self.CC   = ma.masked_where( (self.DEN<1014.5)|(self.TEM > 8.35), self.DEN * (self.TEM-8.35) * 3985 / 1e6 )
            self.CCCn = self.vertint('CC')
            self.CCCn = ma.expand_dims(self.CCCn,1)
        else:
            print(' Case not ready, please code .. ')


############################################################################
# VARIABLE PEA : Potential Energy Anomaly

    def instance_PEA(self, i=None,j=None,k=None):
        if (i is None) and (j is None) :
            self.testvar('TEM')
            self.testvar('SAL')
            self.testvar('DEN')

            self.AVRDEN = np.tile(self.vertmean('DEN'),(1,31,1,1))

            for t in xrange(self.time.shape[0]):
                self.PEAv[t] = 9.81 *self.z*(self.AVRDEN[t]-self.DEN[t])

            self.PEA = self.vertint(PEAv)
            self.PEA= ma.expand_dims(self.PEA,1)

        else:
            print(' Case not ready, please code .. ')


############################################################################                                                                                                             
# VARIABLE : Mixed layer depth  
#
# NOT READY !!! 
#
#        def instance_MLD(t):
#            #  Might probably be ipmroved a lot.                                                                                                                                                                  #             #  Currently, there are two 1D interpolation by i,j
#            Dloc=DEN[t]
#            zloc=z[t]
#            MLDloc=Dloc[30].copy()*0 # for init                                                                                                             
#            for i in xrange(Dloc.shape[1]):
#                for j in xrange(Dloc.shape[2]):
#                    if ma.is_masked(MLDloc[i,j]):
#                        continue
#                    f = interp1d(zloc[:,i,j], Dloc[:,i,j])
#                    d3=f(-3)
#                    f = interp1d(Dloc[:,i,j], zloc[:,i,j])
#            MLDloc[i,j]=f(d3+deltasig)#
#
#            return MLDloc


############################################################################   
# VARIABLE tracer age ... should be made generic for tracer number

    def instance_T1age(self, i=None, j=None, k=None):

        if (i is None) and (j is None) and (k is None):
            self.testvar('T1A')
            self.testvar('T1C')
            self.T1age=ma.masked_where(self.T1C<1e-2, self.T1A/self.T1C)
        else:
            self.gload('T1A',i=i,j=j)
            self.gload('T1C',i=i,j=j)
            exec('T1Aloc=self.T1Ai'+str(i)+'j'+str(j))
            exec('T1Cloc=self.T1Ci'+str(i)+'j'+str(j))
            exec('self.T1agei'+str(i)+'j'+str(j)+'=ma.masked_where(T1Cloc<1e-2, T1Aloc/ T1Cloc)')
