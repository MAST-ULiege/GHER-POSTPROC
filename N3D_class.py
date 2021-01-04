import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import gsw
import yaml
import os.path
import datetime as dt
import G3D_class
import glob

######################################################################
def FullLoad(YAML_FILE = 'local.yml', dstring='1d', begdate=None, enddate=None, filetype='ptrc_T'):
    '''
    Used to prepare and manipulate filelist BEFORE instantiation of the class. 
    '''
    import N3D_class

    try:
        print("\n Full Load from YAML file: %s" % YAML_FILE)
        with open(YAML_FILE, 'r') as stream:
                config = yaml.load(stream)
        resultdir=config['RESULTDIR']
        dstring  =config['DSTRING'] 
    except Exception:
        print("".join(("\n Error in FullLoad : A file called local.yml should be present with RESULTDIR (repertory with model ouptuts) and DSTRING (eg. '1m' or '1d') values","'\n")))

    mlist =  [f for f in glob.glob(resultdir+"*"+dstring+"*"+filetype+"*.nc*")]
    mlist.sort()
    mlist=[m.replace(resultdir,'') for m in mlist]

    if ((begdate is not None) and (enddate is not None)) :
#        try : 
#        for mi,mm in enumerate(mlist):
#            Gl  = N3D_class.N3D(mm,YAML_FILE, instancebat=False)
#            if all(Gl.dates<begdate) or all(Gl.dates>enddate):
#                mlist.drop(mi)
#            if mm==mlist[0]:
#                Ga=Gl
#        else:
#            Ga.dates = ma.append(Ga.dates,Gl.dates,0)

            # Assuming NEMO file conventions
        mlist = [ m for m in mlist if not ((N3D_class.N3D(m,YAML_FILE, instancebat=False).dates[-1]<begdate) or (N3D_class.N3D(m,YAML_FILE, instancebat=False).dates[0]>enddate))]

    return(mlist)

######################################################################

class N3D(G3D_class.G3D): 
    '''
    This is a class for NEMO model outputs exploration. It is based on G3D class but overides specific methods (for z, lon, lat definitions)
    '''
    def __init__(self,infile,YAML_FILE='local.yml', instancebat=True):
        #print(' *******  \n')
        try:
            # Read yaml configuration file
            with open(YAML_FILE, 'r') as stream:
                config = yaml.load(stream)
        except Exception:
            print("".join(("\n A file called local.yml should be present","'\n")))

        try:
            self.model    = config['MODEL'] 
        except :
            self.model    = 'NEMO'
        self.batfile      = config['BATFILE']
        self.figoutputdir = config['PLOTDIR']
        self.resultdir    = config['RESULTDIR']

        if os.path.isfile(infile):
            self.infile=infile
#            print(self.infile + ' -> OK')
            self.found=True
        elif os.path.isfile(config['RESULTDIR']+infile):
            self.infile=config['RESULTDIR']+infile
#            print(self.infile + ' -> OK')
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

        try:
            self.verbose      = True
            self.verbose      = config['VERBOSE']
            self.sparemem     = config['SPAREMEM']
    # instantiate the dictionnary with model parameters
            self.initparams(config['PARAMFILE'])
    # needs an update based on namelist (model specific)
            self.updateparams(config['CONFIGFILE'])
        except: 
            if self.verbose: print('some available setup info not found in the yaml file, check N3D_class.py: _init_, if needed')

        if instancebat:
            self.instance_bat()
        self.testtime()

        self.timevarname  = config['TIMEVARNAME'] if ('TIMEVARNAME' in config) else 'time_counter'
        self.timedimname  = config['TIMEDIMNAME'] if ('TIMEDIMNAME' in config) else 'time_counter'
        self.depthvarname = config['DEPTHVARNAME'] if ('DEPTHVARNAME' in config) else 'deptht'
        self.depthdimname = config['DEPTHDIMNAME'] if ('DEPTHDIMNAME' in config) else 'deptht'
        self.latvarname   = config['LATVARNAME'] if ('LATVARNAME' in config) else 'nav_lat'
        self.latdimname   = config['LATDIMNAME'] if ('LATDIMNAME' in config) else 'y'
        self.lonvarname   = config['LONVARNAME'] if ('LONVARNAME' in config) else 'nav_lon'
        self.londimname   = config['LONDIMNAME'] if ('LONDIMNAME' in config) else 'x'

###################################################################### 
    def timeclean(self, begdate=None, enddate=None):
        '''
        UTILITY : Ensure monotonic, non-redundant dates in an instance of the class, 
        as could result from gathering files with redundant times entries (eg. from model run interruption).

        All fields with valid time dimensions are reduced accordingly. 
        '''
        ld=list(self.dates)   
        indx=[ld.index(i) for i in sorted(np.unique(self.dates))]
        
        otl = len(self.dates)
        # list of array attributes of self with valid time dimension 
        bb  = [ a for a in self.__dict__  if isinstance(self.__getattribute__(a),np.ndarray) ]
        bbb = [ b for b in bb if self.__getattribute__(b).shape[0]==otl]

        for b in bbb:
            print('Time-cleaning '+b)
            self.__setattr__(b, self.__getattribute__(b)[indx])

        if ((begdate is not None) and (enddate is not None)):
            indx = [i for i,d in enumerate(self.dates) if (d>=begdate and d<=enddate)]
            for b in bbb:
                self.__setattr__(b, self.__getattribute__(b)[indx])

######################################################################                                                                                                                                             
    def updateparams(self,f):
        '''
        UTILITY : update model parameters dictionnary with run-specific values         
        '''
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
    def instance_bat(self):
        '''
        VARIABLE : Bathymetry
        '''
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
    def instance_z(self):
        '''
        VARIABLE : Depth
        Dynamic z considering ETA.
        !! Does not currently consider Free Surface !! 
        '''
        with Dataset(self.batfile,'r') as nc:
            self.z   = nc.variables['gdept_0'][:]     # 3D # Wrong for now, we don't consider free surface !!           
            self.lon = nc.variables['glamt'][:][0,0,:]
            self.lat = nc.variables['gphit'][:][0,:,0]
            self.dz  = nc.variables['e3t_0'][:]
######################################################################
    def instance_zi(self):
        '''                                                                                                                                                                                                                VARIABLE : Depth at interfaces                                                                                                                                            
        Dynamic z considering ETA.
        !! Does not currently consider Free Surface !!
        '''
        with Dataset(self.batfile,'r') as nc:
            self.zi  = nc.variables['gdepw'][:]
######################################################################
    def testz(self):
        '''
        UTILITY : test if depth is defined, load it if not.
        '''
        try:
            self.dz
        except :
            if self.verbose: print('dz not found -> loading')
            self.instance_z()
#####################################################################
    def testvar(self,varname,doload=True,i=None,j=None, k=None):
        '''
        UTILITY : test presence of a given variable, load it if needed. 
        '''
        if self.verbose: print( 'Checking presence of v=%s, i=%s, j=%s, k=%s'%(varname,i,j,k))
        if (i is None) and (j is None) and (k is None):
            isthere=hasattr(self,varname)
        elif (k is not None):
            isthere=hasattr(self,varname+'k'+str(k))
        elif (i is not None) and (j is not None):
            isthere=hasattr(self,varname+'i'+str(i)+'j'+str(j))
        if doload and not isthere :
            print ('Loading v=%s for i :%s, j:%s, k:%s'%(varname, i,j,k) )
            if (any([x is not None for x in [i,j,k]])) and self.testvar(varname, doload=False):
                if self.verbose: print ('Getting v=%s for i :%s, j:%s, k:%s from the complete loaded %s'%(varname, i,j,k,varname) )
                if (k is not None) and (i is None) and (j is None):
                    if (k=='bottom'):
#                        exec('self.'+varname+'kbottom=ma.empty_like(self.'+varname+'[:,self.ksurface])[:,None,:,:]')
                        #setattr(self, varname+'kbottom',
                        bid=ma.empty_like(getattr(self,varname)[:,self.ksurface])[:,None,:,:] 
                        for i in range(self.bat.shape[2]):  
                            for j in range(self.bat.shape[3]):
                                if (not ma.is_masked(self.kbottom[0,0,i,j])):
                                    bid[:,0,i,j]=getattr(self,varname)[:,self.kbottom[0,0,i,j],i,j]
#                                    exec('self.'+varname+'kbottom[:,0,i,j]= self.'+varname+'[:,self.kbottom[0,0,i,j],i,j]')
                        setattr(self, varname+'kbottom',bid)
                    elif (k=='surface'):
#                        bid=ma.empty_like(getattr(self,varname)[:,self.ksurface])[:,None,:,:])
 #                           exec('self.'+varname+'ksurface=ma.empty_like(self.'+varname+'[:,self.ksurface])[:,None,:,:]')
                        setattr(self, varname+'ksurface',getattr(self,varname)[:,self.ksurface][:,None,:,:])
                        #exec('self.'+varname+'ksurface=self.'+varname+'[:,self.ksurface][:,None,:,:]')
                    else:
                        setattr(self, varname+'k'+str(k),getattr(self,varname)[:,k][:,None,:,:])
#                            exec('self.'+varname+'k'+str(k)+'=self.'+varname+'[:,k])[:,None,:,:]')
                if (i is not None) and (j is not None):
                    setattr(self, varname+'i'+str(i)+'j'+str(j), getattr(self,varname)[:,:,j,i][:,:,None,None])
#                        exec('self.'+varname+'i'+str(i)+'j'+str(j)+'=self.'+varname+'[:,:,j,i])[:,:,None,None]')
            else:
                self.gload(varname,i=i,j=j,k=k)
                isthere=True
                    
                # NEMO NEEDS EXTRA MASKING
            if (i is None) and (j is None):
                self.testz()
                if self.verbose: print('remasking ' +varname)
                if (k is not None):
                    varname=varname+'k'+k
                setattr(self,varname, ma.masked_array(getattr(self,varname),mask=False))
#                exec('self.'+varname+'=ma.masked_array(self.'+varname+',mask=False)')
                # This ensures that all variables comes with 4 dimension, eventually with length=1
                locshape=getattr(self,varname).shape
#                exec('locshape=self.'+varname+'.shape')
                if len(locshape)!=4: # a dimension is missing. Most probably it's depth for a 2D variable, so we'll deal with that for the moment.
                    if locshape==(len(self.dates),len(self.lat),len(self.lon)):
#                        exec('self.'+varname+'=self.'+varname+'[:,None,:,:]')
                        setattr(self,varname, getattr(self,varname)[:,None,:,:])
                    else:
                        print('missing dim in testvar, for %s of dimensions %s'%(varname,locshape))
                nt=getattr(self,varname).shape[0]
                nz=getattr(self,varname).shape[1]
                bid=getattr(self, varname)
                for t in range(nt):
                    #    if (k is not None):
                    #        exec('self.'+varname+'k'+k+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'k'+k+'.shape[1]]==0, self.'+varname+'k'+k+'[t])')
                    #    else:
                    bid[t]=ma.masked_where(self.landmask[0,:nz]==0,bid[t])
#                    exec('self.'+varname+'[t]=ma.masked_where(self.landmask[0,:self.'+varname+'.shape[1]]==0, self.'+varname+'[t])')

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
            if self.verbose: print('%s not found -> loading'%('time'))
            self.gload('time_counter') # In NEMO time_counter is the number of seconds since referene time 
            self.time=self.time_counter
            del self.time_counter
            with Dataset(self.infile,'r') as nc:
                t0 = nc.variables['time_counter'].time_origin 
            if self.verbose: print('Time Origin : %s'%(t0))    
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

    def instance_U(self,i=None,j=None,k=None):
        ''' VARIABLE : Centered U-velocity - [m/s]
        '''
        # Challenge here is to recompute centered velocity, while it's initially defined for another grid, stored in another file. 
        # First attempt would be to temporarilly change "infile" (and other) properties of the instance to access the other file. 
        print(self.infile)
        self.infile=self.infile.replace('grid_T','grid_U')
        print(self.infile)
        self.testvar('vozocrtx', i=i,j=j, k=k)
        if all ([x is None for x in [i,j,k]]):
            U1=self.vozocrtx
            U2=np.concatenate((U1[:,:,:,1:],U1[:,:,:,[-1]]), axis=3)
            self.U=(U1+U2)/2
            del self.vozocrtx
        elif k is not None:
            U1 = getattr(self,'vozocrtxk'+str(k))
#            exec('U1=self.vozocrtxk'+str(k))
            U2=np.concatenate((U1[:,:,:,1:],U1[:,:,:,[-1]]), axis=3)
#            exec('self.Uk'+k+'=(U1+U2)/2')
            setattr(self,'Uk'+k,(U1+U2)/2)
            delattr(self,'vozocrtxk'+k)
#            exec('del self.vozocrtxk'+k)
        elif (i is not None) and (j is not None):
#            exec('U1=self.vozocrtxi'+str(i)+'j'+str(j))
            U1 = getattr(self,'vozocrtxi'+str(i)+'j'+str(j))
            self.testvar('vozocrtx', i=i+1,j=j, k=k)
            U2 = getattr(self,'vozocrtxi'+str(i+1)+'j'+str(j))
#            exec('U2=self.vozocrtxi'+str(i+1)+'j'+str(j))
#            exec('self.Ui'+str(i)+'j'+str(j)+'=(U1+U2)/2')
            setattr(self,'Ui'+str(i)+'j'+str(j),(U1+U2)/2)

        self.infile=self.infile.replace('grid_U','grid_T')

    def instance_V(self,i=None,j=None,k=None):
        ''' VARIABLE : Centered V-velocity - [m/s]
        '''
        # Challenge here is to recompute centered velocity, while it's initially defined for another grid, stored in another file. 
        # First attempt would be to temporarilly change "infile" (and other) properties of the instance to access the other file. 
        print(self.infile)
        self.infile=self.infile.replace('grid_T','grid_V')
        print(self.infile)
        self.testvar('vomecrty', i=i,j=j, k=k)
        if all ([x is None for x in [i,j,k]]):
            V1=self.vomecrty
            V2=np.concatenate((V1[:,:,:,1:],V1[:,:,:,[-1]]), axis=2)
            self.V=(V1+V2)/2
            del self.vomecrty
        elif k is not None:
            V1 = getattr(self,'vomecrtyk'+str(k))
#            exec('V1=self.vomecrtyk'+str(k))
            V2=np.concatenate((V1[:,:,1:,:],V1[:,:,[-1],:]), axis=2)
            setattr(self,'Vk'+k,(V1+V2)/2)
            delattr(self,'vomecrtyk'+k)
#            exec('self.Vk'+k+'=(V1+V2)/2')
#            exec('del self.vomecrtyk'+k)
        elif (i is not None) and (j is not None):
#            exec('V1=self.vomecrtyi'+str(i)+'j'+str(j))
            V1=getattr(self,'vomecrtyi'+str(i)+'j'+str(j))
            self.testvar('vomecrty', i=i,j=j+1, k=k)
            V2=getattr(self,'vomecrtyi'+str(i)+'j'+str(j+1))
#            exec('V2=self.vomecrtyi'+str(i)+'j'+str(j+1))
            setattr(self,'Vi'+str(i)+'j'+str(j),(V1+V2)/2)
#            exec('self.Vi'+str(i)+'j'+str(j)+'=(V1+V2)/2')

        self.infile=self.infile.replace('grid_V','grid_T')

############################################################################
    '''    
    def avgprofile(self,varname,
                   ztab=-1*np.concatenate([np.arange(0,10,2), np.arange(10,40,5),np.arange(50,120,10),np.arange(120,300,50),np.arange(300,1000,200)]),
                   maskin=None
                   ):
        
        #The idea is to get an average profile as a function of time                                                                                                                                               
        #return is 2D : [z,time]                                                                                                                                                                                   
        #sigma-space makes it a bit complicate                                                                                                                                                           
         

        self.testz()
        self.testvar(varname)
        self.testtime()

        avg = ma.empty((len(self.time),ztab.shape[0]-1))

        loc=getattr(self,varname)

        if maskin is not None:
            loc = self.maskvar(loc,maskin)

        if (len(self.dz.shape)==3)and(len(loc.shape)==4):
            gridZU = self.zi[0,1:]
            gridZD = self.zi[0,:-1]
            for k in range(ztab.shape[0]-1):
                if self.verbose: print('%s / %s'%(k+1,ztab.shape[0]-1))
                dzloc= ma.maximum(ma.zeros(self.dz.shape), np.minimum(gridZU, ztab[k])-np.maximum(gridZD, ztab[k+1]))
                for t in range(len(self.time)):
                    vol=ma.masked_where(loc[t].mask,dzloc*self.dy*self.dx)
                    bi=loc[t]*vol
                    #avg[t,k]= ma.sum(bi)/ma.sum(vol)                                                                                                                                                               
                    avg[t,k]= ma.masked_where( (float(bi.count())/float(bi.size))<0.1,ma.sum(bi)/ma.sum(vol))
                    if (t==10):
                        if self.verbose: print('k %s : from %s to %s' %(k,ztab[k],ztab[k+1]))
                        if self.verbose: print('total volume considered for this layer :  %s km^3' %(ma.sum(vol)/1e9))
                        if self.verbose: print('mean val :  %s' %(avg[t,k]))

        zforplot=(ztab[:-1]+ztab[1:])/2
        return avg, zforplot
    '''
