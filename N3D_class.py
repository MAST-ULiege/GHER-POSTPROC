import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt

import G3D_class

class N3D(G3D_class.G3D): 
    '''This is a class for NEMO model outputs exploration. It is based on G3D class but overides specific methods (for z, lon, lat definitions)'''
######################################################################
    def __init__(self,infile,YAML_FILE='local.yml'):
        print(' *******  \n')
        try:
#            YAML_FILE = 'local.yml'
            print("\nLaunching with YAML file: %s" % YAML_FILE)
            # Read yaml configuration file
            with open(YAML_FILE, 'r') as stream:
                config = yaml.load(stream)
        except Exception:
            print("".join(("\n A file called local.yml should be present","'\n")))

        try:
            self.model    = config['MODEL'] 
        except :
            self.model    = 'GHER'
        self.batfile      = config['BATFILE']
        self.figoutputdir = config['PLOTDIR']
        self.resultdir    = config['RESULTDIR']

        if os.path.isfile(infile):
            self.infile=infile
            print(self.infile + ' -> OK')
            self.found=True
        elif os.path.isfile(config['RESULTDIR']+infile):
            self.infile=config['RESULTDIR']+infile
            print(self.infile + ' -> OK')
            self.found=True
        else:
            print(infile+' can not be found in '+config['RESULTDIR'])
            self.found=False
            return

        if not os.path.isdir(self.figoutputdir):
            os.mkdir(self.figoutputdir)

        #The diag filename can be used to store computed diagnostic
        try: 
            self.diagdir    = config['DIAGDIR']
            self.diagfile   = self.diagdir+infile[:-3]+".diag.nc"
        except:
            self.diagfile= self.infile[:-3]+".diag.nc"

        self.instance_bat()
        self.testtime()

        try:
            self.verbose      = config['VERBOSE']
            self.sparemem     = config['SPAREMEM']
    # instantiate the dictionnary with model parameters
            self.initparams(config['PARAMFILE'])
    # needs an update based on namelist (model specific)
            self.updateparams(config['CONFIGFILE'])
        except: 
            print('some available setup info not found in the yaml file, check N3D_class.py: _init_, if needed')


        self.timevarname  = config['TIMEVARNAME'] if ('TIMEVARNAME' in config) else 'time_counter'
        self.timedimname  = config['TIMEDIMNAME'] if ('TIMEDIMNAME' in config) else 'time_counter'
        self.depthvarname = config['DEPTHVARNAME'] if ('DEPTHVARNAME' in config) else 'deptht'
        self.depthdimname = config['DEPTHDIMNAME'] if ('DEPTHDIMNAME' in config) else 'deptht'
        self.latvarname   = config['LATVARNAME'] if ('LATVARNAME' in config) else 'nav_lat'
        self.latdimname   = config['LATDIMNAME'] if ('LATDIMNAME' in config) else 'y'
        self.lonvarname   = config['LONVARNAME'] if ('LONVARNAME' in config) else 'nav_lon'
        self.londimname   = config['LONDIMNAME'] if ('LONDIMNAME' in config) else 'x'




###################################################################### 

######################################################################                                                                                                                                             
# UTILITY : update param dictionnary with run-specific values                                                                                                                                                  
    def updateparams(self,f):
        fp = open(f, 'r')
        line='a'
        while 'namtrc_bamhbi' not in line: # skip all lines before relevant ones
            line = fp.readline()
        pline='='
        while pline:
            pline=fp.readline().split('!')[0]
            if '=' not in pline:
                break
            p,v=pline.split('=')
            exec('v='+v)
            self.paramd[p.strip()]=float(v)


######################################################################
# VARIABLE : BAT
    def instance_bat(self):
        if os.path.isfile(self.batfile):
            print(self.batfile + ' -> OK')
        else:
            print(self.batfile + ' can not be found') 

        self.ksurface=0


        with Dataset(self.batfile,'r') as nc:
            self.kbottom = nc.variables['mbathy'][:]  # 3D
            self.kbottom = self.kbottom[:,None,:,:]   # 4D
            self.z = nc.variables['gdept_0'][:]       # 3D # Wrong for now, we don't consider free surface !! 
            print('gdept_0 shape')
            print(self.z.shape)

            self.bat=ma.empty_like(self.kbottom[:])
            for i in range(self.bat.shape[3]):
                for j in range(self.bat.shape[2]):
                    self.bat[0,0,j,i]=self.z[0,self.kbottom[0,0,j,i],j,i]

            self.bat= ma.masked_where(self.kbottom==0,self.bat)

            self.dx= nc.variables['e1t'][:]
            self.dy= nc.variables['e2t'][:]
            
            self.landmask= nc.variables['tmask'][:]   # 3D 

        self.kbottom = ma.where(self.kbottom!=0, self.kbottom-1,self.kbottom)
######################################################################
# VARIABLE : Z
    def instance_z(self):
        # Dynamic z considering ETA

        with Dataset(self.batfile,'r') as nc:
            self.z   = nc.variables['gdept_0'][:]     # 3D # Wrong for now, we don't consider free surface !!           
            self.lon = nc.variables['glamt'][:][0,0,:]
            self.lat = nc.variables['gphit'][:][0,:,0]
            self.dz  = nc.variables['e3t_0'][:]
######################################################################
# VARIABLE : ZI
    def instance_zi(self):
        # For now let's build a constant dz,
        with Dataset(self.batfile,'r') as nc:
            self.zi  = nc.variables['gdepw'][:]
######################################################################
# UTILITY : test z        
    def testz(self):
        try:
            self.dz
        except :
            print('dz not found -> loading')
            self.instance_z()
#####################################################################
# UTILITY : test variable
    def testvar(self,varname,doload=True,i=None,j=None, k=None, verbose=True):
        try:
            print( 'Checking presence of v=%s, i=%s, j=%s, k=%s'%(varname,i,j,k))
            if (i is None) and (j is None) and (k is None):
                exec('self.'+varname)
            elif (k is not None):
                exec('self.'+varname+'k'+str(k))
            elif (i is not None) and (j is not None):
                print('self.'+varname+'i'+str(i)+'j'+str(j))
                exec('self.'+varname+'i'+str(i)+'j'+str(j))
            isthere=True
        except:
            print('Not Loaded: of v=%s, i=%s, j=%s, k=%s'%(varname,i,j,k))
            isthere=False
            if doload:
                if (any([x is not None for x in [i,j,k]])) and self.testvar(varname, doload=False, verbose=False):
                    print ('Getting v=%s for i :%s, j:%s, k:%s from the complete loaded %s'%(varname, i,j,k,varname) )
                    if (k is not None) and (i is None) and (j is None):
                        if (k=='bottom'):
                            exec('self.'+varname+'kbottom=ma.empty_like(self.'+varname+'[:,self.ksurface])[:,None,:,:]')
                            for i in range(self.bat.shape[2]):  
                                for j in range(self.bat.shape[3]):
                                    if (not ma.is_masked(self.kbottom[0,0,i,j])):
                                        exec('self.'+varname+'kbottom[:,0,i,j]= self.'+varname+'[:,self.kbottom[0,0,i,j],i,j]')
                        elif (k=='surface'):
                            exec('self.'+varname+'ksurface=ma.empty_like(self.'+varname+'[:,self.ksurface])[:,None,:,:]')
                            exec('self.'+varname+'ksurface=self.'+varname+'[:,self.ksurface][:,None,:,:]')
                        else: 
                            exec('self.'+varname+'k'+str(k)+'=self.'+varname+'[:,k])[:,None,:,:]')
                    if (i is not None) and (j is not None):
                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'=self.'+varname+'[:,:,j,i])[:,:,None,None]')
                else:
                    self.gload(varname,i=i,j=j,k=k)
                    isthere=True
                    
                # NEMO NEEDS EXTRA MASKING
                if (i is None) and (j is None):
                    self.testz()
                    print('remasking ' +varname)
                    if (k is not None):
                        varname=varname+'k'+k
                    exec('self.'+varname+'=ma.masked_array(self.'+varname+',mask=False)')
                    # This ensures that all variables comes with 4 dimension, eventually with length=1
                    exec('locshape=self.'+varname+'.shape')
                    if len(locshape)!=4: # a dimension is missing. Most probably it's depth for a 2D variable, so we'll deal with that for the moment.
                        if locshape==(len(self.dates),len(self.lat),len(self.lon)):
                            exec('self.'+varname+'=self.'+varname+'[:,None,:,:]')
                        else:
                            print('missing dim in testvar, for %s of dimensions %s'%(varname,locshape))
                    exec('nt=self.'+varname+'.shape[0]')
                    for t in range(nt):
                    #    if (k is not None):
                    #        exec('self.'+varname+'k'+k+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'k'+k+'.shape[1]]==0, self.'+varname+'k'+k+'[t])')
                    #    else:
                        exec('self.'+varname+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'.shape[1]]==0, self.'+varname+'[t])')

                #elif (k is not None):
                #    print('remasking ' +varname+str(k))
#                else:
#                    print('remasking for N vars is not yet implemented for cases where i or j are not None')
        return(isthere)
######################################################################
# UTILITY : test time
    def testtime(self):
        try:
            self.time
        except:
            print('%s not found -> loading'%('time'))
            self.gload('time_counter') # In NEMO time_counter is the number of seconds since referene time 
            self.time=self.time_counter
            del self.time_counter
            with Dataset(self.infile,'r') as nc:
                t0 = nc.variables['time_counter'].time_origin 
            print('Time Origin : %s'%(t0))    
            self.dates = [dt.datetime.strptime(t0,'%Y-%m-%d %H:%M:%S')+dt.timedelta(seconds=int(t)) for t in self.time]
            self.timevarname='time_counter'

######################################################################
# UTILITY : 
    def instance_SAL(self,i=None,j=None, k=None):
        self.testvar('vosaline',i=i,j=j, k=k)
        if k is not None:
            exec('self.SALk'+str(k)+'=self.vosalinek'+str(k))
            exec('del self.vosalinek'+str(k))
        elif (i is not None) and (j is not None):
            exec('self.SALi'+str(i)+'j'+str(j)+'=self.vosalinei'+str(i)+'j'+str(j))
            exec('del self.vosalinei'+str(i)+'j'+str(j))
        else: 
            self.SAL=self.vosaline
            del self.vosaline

    def instance_TEM(self,i=None,j=None,k=None):
        self.testvar('votemper',i=i,j=j, k=k)
        if k is not None:
            exec('self.TEMk'+str(k)+'=self.votemperk'+str(k))
            exec('del self.votemperk'+str(k))
        elif (i is not None) and (j is not None): 
            exec('self.TEMi'+str(i)+'j'+str(j)+'=self.votemperi'+str(i)+'j'+str(j))
            exec('del self.votemperi'+str(i)+'j'+str(j))
        else:
            self.TEM=self.votemper
            del self.votemper
