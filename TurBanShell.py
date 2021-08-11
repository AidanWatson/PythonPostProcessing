#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 10:37:40 2021

@author: watsonai
"""

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import Ofpp
from os.path import basename, realpath, exists

def IndexSearch(WN, Val):#returns index of a point (Val) in array (WN)
    i=0#index of current closest fit
    j=0
    for data in WN:
        if(abs(data-Val)<=abs(WN[i]-Val)):
            i=j
        j=j+1
    return i

class openfoam(object):
    """
        Takes as input a directory name under shelldirname, and a timeset from filenum (usually 0 to load variables) using cylinderAdjustFormat
        sets internal shelldirname and filenum to given values
        Main object only targets one file
        Loadslice takes values from vars2l, which is a dictionary of {'variablename':'b'} where b is a boundary field and v is a volume field
        Loadslice requires OF to print Cxyz files for cell centre, and stores a dictionary of the form
        {'Var':'array'}
        Cell centres can be loaded as volume fields, or using mesh.read_cell_centres, then mesh.cell_centres,
        which returns an array ((cx1,cy1,cz1),(cx2...))
        running loadslice takes the timestep as an argument
    
        LoadLog Method parses dsmcLog file, and outputs values into __dict__ with the values as keys, alongside time data
        LoadLog uses parseTime and ParseLog to obtain time arrays and other arrays
        
        MakeAvg creates an array of time averaged data
        Make Avg needs the key to the array to be averaged as an argument, and makes an array in __dict__ with key+Avg as the name
    """
    #file defined somehow
    #filedir='/home/watsonai/OpenFOAM/OpenFOAM-v1706/FaunFiles/CylinderAdjust1/'
    
    def __init__(self,shelldirname=None,filenum=None):#shelldirname:name of FileDirectory filenum=file forces are found in
        # If no rundir specified
        self.mesh=None
        if shelldirname is None: 
            shelldirname = input('Please enter the rundir: ') 
            self.rundir = realpath(shelldirname)
            self.dirname= basename(self.rundir)
        else:
            self.shelldirname=shelldirname
      
          
          # If filenum not given
          
        if filenum is None:
           self.filenum=input("Please enter the file number to load (e.g. 000): ")
        else:
           self.filenum=filenum
        # self.fileDir=self.shelldirname()+self.filenum #sets FileDir to name of target file
           
           
        ## Load run specific parameters using Ofpp
        ## Make sure to have all the run metadata available
        ## e.g. self.lx = length of simulation in x direction
        ##      self.npoints = number of mesh points
        ##      self.mesh = ...
        ##       Other useful quantities
        ##      
    def vars2load(self,v2lu):
        self.vars2l=v2lu
    def TargetFile(self,fname):
        self.fileDir=fname
    
    def loadslice(self,num):#num=> time being loaded
        if self.mesh==None:
            print(self.shelldirname)
            self.mesh=Ofpp.FoamMesh(str(self.shelldirname))
            
        for vrbl in self.vars2l.keys():
            print(vrbl)
            # self.__dict__[vrbl] = Ofpp.parse_boundary_field(fdFile)## USE OFPP TO LOAD VRBL
            targFile=str(self.shelldirname)+'/'+str(num)+'/'+str(vrbl)
            if self.vars2l[vrbl]=='v':
                self.__dict__[vrbl] = Ofpp.parse_internal_field(targFile)## USE OFPP TO LOAD VRBL
            elif self.vars2l[vrbl]=='b':

                self.__dict__[vrbl]=Ofpp.parse_boundary_field(targFile)[b'cylinder'][b'value']
            else:
                print("invalid field type")
                
    def loadLog(self):#parses dsmclog file, reads out time, values of variables in log into dict with keys specified by dictNames
        LogOut=ParseLog(self.shelldirname)
        dictNames=LogOut[0]
        for i in range(len(LogOut)-1):
            self.__dict__[dictNames[i]]=LogOut[i+1]
        print(self.__dict__.keys())
        
        
    def makeAvgScal(self, key=None): #makes array of average of a parameter and puts it in __dict__ with the key "oldkey+Mean"
        if key==None:
            key=input('Key required to make average')
        times=self.__dict__['STime']
        TargArr=self.__dict__[key]
        MeanArr=[]
        dt=times[len(times)-1]-times[len(times)-2]

        for i in range(len(times)):
            if i>=2:
                runtime=times[i-1]

                M1=(1/(runtime+dt))*(runtime*MeanArr[i-1]+dt*TargArr[i])
                MeanArr.append(M1)
            else:
                MeanArr.append(TargArr[i])
        NewKey=key+'Mean'
        self.__dict__[NewKey]=MeanArr
        
    def makeAvgvec(self, key=None): #makes array of average of a vector parameter and puts it in __dict__ with the key "oldkey+Mean"
        if key==None:
            key=input('Key required to make average')
        times=self.__dict__['STime']
        TargArr=self.__dict__[key]
        MAX=[]
        MAY=[]
        MAZ=[]
        dt=times[len(times)-1]-times[len(times)-2]
        TAX=TargArr[0]
        TAY=TargArr[1]
        TAZ=TargArr[2]

        for i in range(len(times)):
            if i>=2:
                runtime=times[i-1]

                MX1=(1/(runtime+dt))*(runtime*MAX[i-1]+dt*TAX[i])
                MY1=(1/(runtime+dt))*(runtime*MAY[i-1]+dt*TAY[i])
                MZ1=(1/(runtime+dt))*(runtime*MAZ[i-1]+dt*TAZ[i])

                MAX.append(MX1)
                MAY.append(MY1)
                MAZ.append(MZ1)


            else:
                MAX.append(TAX[i])
                MAY.append(TAY[i])
                MAZ.append(TAZ[i])


        NewKey=key+'Mean'
        MeanArr=[MAX,MAY,MAZ]
        self.__dict__[NewKey]=MeanArr
            
    
    def makeRAvg(self, key=None,N=1):
        if key==None:
            key=input('Key required to make average')
        TargArr=self.__dict__[key]
        MeanArr=AvgFile(TargArr,N)
        NewKey=key+'Mean'
        self.__dict__[NewKey]=MeanArr        
        print(len(MeanArr))
        
    def makeRAvgVec(self,key=None, N=1):
        if key==None:
            key=input('Key required to make average')
        TargArr=self.__dict__[key]
        MeanArr=[]
        for i in range(len(TargArr)):
            MeanArr.append(AvgFile(TargArr[i],N))
        NewKey=key+'Mean'
        self.__dict__[NewKey]=MeanArr    
            
        


"""
Parses .dat file of integrated forces into array of |time|force in x||fy|fz||mean force in x|mfy|mfz

"""
def loadForce(fname,S=600):#parses .dat file of force into separate arrays for force and time
    try:
        forcefile=str(fname)+'postProcessing/forces1/0/force_0.dat'
        infile = open(forcefile, 'r')
    except:
        print(str(fname)+" does not have a force0 file")#catch for weirdness where some files only generate a forces file (and some generate identical forces files, and only forces0 contains the useful information
        forcefile=str(fname)+'postProcessing/forces1/0/force.dat'
        infile = open(forcefile, 'r')
        print('used force file instead')
        
        # print('checked')
    # infile = open(forcefile, 'r')
    # print('doublechecked')
    skip=4
    times=[]
    fdx=[]
    fdy=[]
    fdz=[]
    fdxMean=[]
    fdyMean=[]
    fdzMean=[]
    i=0
    for line in infile:
        string=line.split('\t')
        if(i>=skip):

            times.append(float(string[0]))#parses time, fd and adds to array
            fdxi=float(string[1].strip("(").strip(")").split(' ')[0])
            fdyi=float(string[1].strip("(").strip(")").split(' ')[1])
            fdzi=float(string[1].strip("(").strip(")").split(' ')[2])
            fdx.append(fdxi)
            fdy.append(fdyi)
            fdz.append(fdzi)
            
            if(len(times)>2):#calculates means

                dt=times[len(times)-1]-times[len(times)-2]
                runtime=times[len(times)-1]
                
                WFdx=(1/(runtime+dt))*(runtime*fdxMean[len(fdxMean)-1]+dt*fdxi)
                fdxMean.append(WFdx)
                WFdy=(1/(runtime+dt))*(runtime*fdyMean[len(fdyMean)-1]+dt*fdyi)
                fdyMean.append(WFdy)
                WFdz=(1/(runtime+dt))*(runtime*fdzMean[len(fdzMean)-1]+dt*fdzi)
                fdzMean.append(WFdz)
            else:
                fdxMean.append(fdxi)
                fdyMean.append(fdyi)
                fdzMean.append(fdzi)
           

        i+=1
    fdxRA=AvgFile(fdx,S)
    fdyRA=AvgFile(fdy,S)
    fdzRA=AvgFile(fdz,S)
    
    return [times,fdx,fdy,fdz,fdxRA,fdyRA,fdzRA,fdxMean,fdyMean,fdzMean]

def loadForcePD(fname,S=600):#parses .dat file of force from pdfoam format into separate arrays for force and time
    try:
        forcefile=str(fname)+'postProcessing/forces1/0/force_0.dat'
        infile = open(forcefile, 'r')
    except:
        print(str(fname)+" does not have a force0 file")#catch for weirdness where some files only generate a forces file (and some generate identical forces files, and only forces0 contains the useful information
        forcefile=str(fname)+'postProcessing/forces1/0/force.dat'
        infile = open(forcefile, 'r')
        print('used force file instead')
        
        # print('checked')
    # infile = open(forcefile, 'r')
    # print('doublechecked')
    skip=4
    times=[]
    fdx=[]
    fdy=[]
    fdz=[]
    fdxA=[]
    fdyA=[]
    fdzA=[]
    fdxMean=[]
    fdyMean=[]
    fdzMean=[]
    i=0
    for line in infile:
        string=line.split('\t')
        if(i>=skip):

            times.append(float(string[0]))#parses time, fd and adds to array
            fdxiA=float(string[1].strip("(").strip(")").split(' ')[0])
            fdyiA=float(string[1].strip("(").strip(")").split(' ')[1])
            fdziA=float(string[1].strip("(").strip(")").split(' ')[2])
            fdxA.append(fdxiA)
            fdyA.append(fdyiA)
            fdzA.append(fdziA)#A inserted to represent unadjusted values of fd, before adjusting to take gradient
            fdxi=fdxiA-fdxA[i-skip-1]
            fdyi=fdyiA-fdyA[i-skip-1]
            fdzi=fdziA-fdzA[i-skip-1]

            fdx.append(fdxi)
            fdy.append(fdyi)
            fdz.append(fdzi)
            
            if(len(times)>2):#calculates means

                dt=times[len(times)-1]-times[len(times)-2]
                runtime=times[len(times)-1]
                
                WFdx=(1/(runtime+dt))*(runtime*fdxMean[len(fdxMean)-1]+dt*fdxi)
                fdxMean.append(WFdx)
                WFdy=(1/(runtime+dt))*(runtime*fdyMean[len(fdyMean)-1]+dt*fdyi)
                fdyMean.append(WFdy)
                WFdz=(1/(runtime+dt))*(runtime*fdzMean[len(fdzMean)-1]+dt*fdzi)
                fdzMean.append(WFdz)
            else:
                fdxMean.append(fdxi)
                fdyMean.append(fdyi)
                fdzMean.append(fdzi)
           

        i+=1
    fdxRA=AvgFile(fdx,S)
    fdyRA=AvgFile(fdy,S)
    fdzRA=AvgFile(fdz,S)
    
    return [times,fdx,fdy,fdz,fdxRA,fdyRA,fdzRA,fdxMean,fdyMean,fdzMean]
"""
Plots force from a selection of cases against a parameter
needs array of parameter values, and of file numbers (assumes files are stored in file fname, with names differing only by the number at the end)
pltGraphsConv and pltGraphSum default to false, and represent plotting the convergence graphs and final plot of parameter vs drag force, respectively.
currently needs file numbers to be in ascending order of parameter, but eventually should sort by parameter
"""
def PFAnalyse(Parameters,CANumber,paramName=None,units=None,PltGraphsConv=False,pltGraphSum=False,fname=None):#parameter for analysis, number of cylinderAdjustFile, bool for plotting convergence plots, bool for plotting summary scatter plot.
    #plots fd0 from array of indices (corresponding to cylinderadjust files)
    bP=Parameters
    fdxFinal=[]
    fdyFinal=[]
    if fname==None:
        fname='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CylinderAdjust'
    if paramName==None:
        paramName="UnknownParameter"
    if units==None:
        units='arb.'
    for i in range(len(CANumber)):
        ind=CANumber[i]
        fdArrays=loadForce(fname+ind+'/')
        time=fdArrays[0]
        if(PltGraphsConv):
            fig1 = plt.figure()#plots fd convergences
            plt.plot(fdArrays[0],fdArrays[1],label='fdx')
            plt.plot(fdArrays[0],fdArrays[4])
            plt.plot(fdArrays[0],fdArrays[7])

            plt.plot(fdArrays[0],fdArrays[2],label='fdy')
            plt.plot(fdArrays[0],fdArrays[5])
            plt.plot(fdArrays[0],fdArrays[8])

            plt.title('Case Number '+ind+' ('+str(bP[i])+')')
            plt.xlabel('Time')
            plt.ylabel('Force')
            plt.legend()
            plt.show()
        
        final=IndexSearch(time,max(time))
        print('fdxFinal='+str(fdArrays[4][final]))
        fdxFinal.append(fdArrays[4][final])
        fdyFinal.append(fdArrays[5][final])
    if(pltGraphSum):
        fig5,(ax1,ax2)=plt.subplots(1,2)
        ax1.plot(bP,fdxFinal,'--o', color='darkGreen',linewidth=0.5,label='fdX')
        ax2.plot(bP,fdyFinal,'--o', color='darkOrange',linewidth=0.5,label='fdY')
        ax1.grid(alpha=0.5)
        ax2.grid(alpha=0.5)
        ax1.set_title('fdX')
        ax2.set_title('fdy')
        ax1.set_ylabel('Force(N)')
        ax2.set_ylabel('Force(N)')
        ax1.set_xlabel(paramName+" ("+units+')')
        ax2.set_xlabel(paramName+" ("+units+")")
    print(fdxFinal)
    print(fdyFinal)
    
    plt.show()

def UFAnalyse(velocities,CANumber,PltGraphsConv=False,pltGraphSum=False,fname=None):#velocities for analysis, number of cylinderAdjustFile, bool for plotting convergence plots, bool for plotting summary scatter plot.
    #plots fd0 from array of indices (corresponding to cylinderadjust files)
    bU=velocities
    fdxFinal=[]
    fdyFinal=[]
    if fname==None:
        fname='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CylinderAdjust'
    for ind in CANumber:
        fdArrays=loadForce(fname+ind+'/')
        time=fdArrays[0]
        if(PltGraphsConv):
            fig1 = plt.figure()#plots fd convergences
            plt.plot(fdArrays[0],fdArrays[1],label='fdx')
            plt.plot(fdArrays[0],fdArrays[4])
            plt.plot(fdArrays[0],fdArrays[2],label='fdy')
            plt.plot(fdArrays[0],fdArrays[5])
            plt.title('CylinderAdjust'+ind)
            plt.xlabel('Time')
            plt.ylabel('Force')
            plt.legend()
            plt.show()
        
        final=IndexSearch(time,max(time))
        print('fdxFinal='+str(fdArrays[4][final]))
        fdxFinal.append(fdArrays[4][final])
        fdyFinal.append(fdArrays[5][final])
    if(pltGraphSum):
        fig5,(ax1,ax2)=plt.subplots(1,2)
        ax1.plot(bU,fdxFinal,'--o', color='darkGreen',linewidth=0.5,label='fdX')
        ax2.plot(bU,fdyFinal,'--o', color='darkOrange',linewidth=0.5,label='fdY')
        ax1.grid(alpha=0.5)
        ax2.grid(alpha=0.5)
        ax1.set_title('fdX')
        ax2.set_title('fdy')
        ax1.set_ylabel('Force(N)')
        ax2.set_ylabel('Force(N)')
        ax1.set_xlabel('Flow speed(km/s)')
        ax2.set_xlabel('Flow speed(km/s)')
    print(fdxFinal)
    print(fdyFinal)
    
    plt.show()
    
"""
takes dsmc.log file and outputs total runtime parseTime, plots with runTimeCheck
parseTime returns total runtime, simulation time
RuntimeCheck takes filename and number (assumes files are numbered as cylinderAdjust files)
also takes NArray and MakeNGraphs, for an array of NEquiv values, and boolean to plot runtime against NEquiv
returns array of final times
"""

def ParseTime(fname=None):
    if fname==None:
        print('no file given to parseTime')
        fname='  '
    try:
        forcefile=str(fname)+'log.pdFoam'
        infile = open(forcefile, 'r')
    except:
        print('no log file in '+fname)
    ExTimeArray=[]#execution time
    ClockTimeArray=[]#clockTime
    SimTimeArray=[]#simulation time

    for line in infile:
        if 'Time =' in line:
            if 'Execution' in line:
                splitline=line.split(' ')
                ExTimeArray.append(float(splitline[2]))
                ClockTimeArray.append(float(splitline[7]))
            else:
                splitline=line.split(' ')
                SimTimeArray.append(float(splitline[2].split('\n')[0] )) 
    returnArray=[ExTimeArray,ClockTimeArray,SimTimeArray]
    return returnArray

    
"""
ParseLog takes filename as input, and outputs array 
|time|particles inserted|numCollisions|numParticles|NumMolecules|SystemMass |px|py|pz|PMag|KEAvg|IEAvg|ETotAv|ExTime
"""
def ParseLog(fname=None, PD=False):
    ForceCheck=True#janky way to distinguish between force and moment data
    if fname==None:
        fname=input('type file for parseLog to check')
    try:
        if(PD==False):
            logfile=str(fname)+'log.dsmcFoam'
        else:
            logfile=str(fname)+'log.pdFoam'
        infile = open(logfile, 'r')
    except:
        print('no log file in '+fname)
    if True:
        STime=[]
        PartIns=[]
        NumColl=[]
        NumP=[]
        NumMol=[]
        SysMass=[]
        px=[]
        py=[]
        pz=[]
        pMag=[]
        KEAvg=[]
        IEAvg=[]
        ETotAvg=[]
        ExTime=[]
        
        FTX=[]
        FTY=[]
        FTZ=[]
        
        FPX=[]
        FPY=[]
        FPZ=[]
        
        FVX=[]
        FVY=[]
        FVZ=[]
        
        MTX=[]
        MTY=[]
        MTZ=[]
        
        MPX=[]
        MPY=[]
        MPZ=[]
        
        MVX=[]
        MVY=[]
        MVZ=[]
        
    for line in infile:
        if 'Time =' in line:
            if 'Execution' in line:
                splitline=line.split(' ')
                ExTime.append(float(splitline[2]))
            else:
                splitline=line.split(' ')
                STime.append(float(splitline[2].split('\n')[0]))

        if 'inserted' in line:#checks for particles inserted
            splitline=line.split('= ')
            PartIns.append(float(splitline[1].split('\n')[0]))
        
        if 'collisions' in line:#checks for collissions
            if '=' in line:
               splitline=line.split('= ')
               NumColl.append(float(splitline[1].split('\n')[0]))
            else:
                NumColl.append(0)
                
        if 'dsmc particles' in line:
            splitline=line.split('= ')
            NumP.append(float(splitline[1].split('\n')[0]))
        
        
        if 'molecules' in line:
            splitline=line.split('= ')
            NumMol.append(float(splitline[1].split('\n')[0]))            
            

        if 'in system' in line:
            splitline=line.split('= ')
            SysMass.append(float(splitline[1].split('\n')[0]))
        
        if 'momentum|' in line:
            splitline=line.split('= ')
            pMag.append(float(splitline[1].split('\n')[0]))

        if 'kinetic' in line:
            splitline=line.split('= ')
            KEAvg.append(float(splitline[1].split('\n')[0]))
 
        if 'internal ' in line:
            splitline=line.split('= ')
            # print(splitline)
            IEAvg.append(float(splitline[1].split('\n')[0]))
 
        if 'total energy' in line:
            splitline=line.split('= ')
            ETotAvg.append(float(splitline[1].split('\n')[0]))

        if 'momentum ' in line:
            splitline=line.split('(')[1].split(')')[0].split(' ')
            px.append(float(splitline[0]))
            py.append(float(splitline[1]))
            pz.append(float(splitline[2]))
            
        if 'of forces' in line:#bad solution please dont judge me, future me
            ForceCheck=True
        if 'of moments' in line:
            ForceCheck=False
            
        if 'Total    :' in line:
            splitline=line.split('(')[1].split(')')[0].split(' ')
            if(ForceCheck):
                FTX.append(float(splitline[0]))
                FTY.append(float(splitline[1]))
                FTZ.append(float(splitline[2]))
            else:
                MTX.append(float(splitline[0]))
                MTY.append(float(splitline[1]))
                MTZ.append(float(splitline[2]))   
                
        if 'Pressure :' in line:
            splitline=line.split('(')[1].split(')')[0].split(' ')
            if(ForceCheck):
                FPX.append(float(splitline[0]))
                FPY.append(float(splitline[1]))
                FPZ.append(float(splitline[2]))
            else:
                MPX.append(float(splitline[0]))
                MPY.append(float(splitline[1]))
                MPZ.append(float(splitline[2])) 
            
        if 'Viscous  :' in line:
            splitline=line.split('(')[1].split(')')[0].split(' ')
            if(ForceCheck):
                FVX.append(float(splitline[0]))
                FVY.append(float(splitline[1]))
                FVZ.append(float(splitline[2]))
            else:
                MVX.append(float(splitline[0]))
                MVY.append(float(splitline[1]))
                MVZ.append(float(splitline[2]))                 
        
        
        
        
        
    DictEntries=['STime', 'PartIns', 'NumColl', 'NumP','NumMol', 'SysMass','px','py','pz', 'pMag', 'KEAvg', 'ETotAvg','FTot','FPress','FVisc','MTot','MPress','MVisc', 'ExTime']
    StuffArray=[DictEntries,STime, PartIns, NumColl, NumP,NumMol, SysMass,px,py,pz, pMag, KEAvg, ETotAvg, [FTX,FTY,FTZ],[FPX,FPY,FPZ],[MVX,MVY,MVZ],[MTX,MTY,MTZ],[MPX,MPY,MPZ],[MVX,MVY,MVZ],ExTime]
    return StuffArray

def CheckConvergence(fname=None,S=100):
    fdArrays=loadForce(fname,S)
    fig1 = plt.figure()#plots fd convergences
    plt.plot(fdArrays[0],fdArrays[1],label='fdx')
    plt.plot(fdArrays[0],fdArrays[4],label='fdx (RA)')
    plt.plot(fdArrays[0],fdArrays[7],label='fdx (M)')
    plt.plot(fdArrays[0],fdArrays[2],label='fdy')
    plt.plot(fdArrays[0],fdArrays[5],label='fdy (RA)')
    plt.plot(fdArrays[0],fdArrays[8],label='fdy (M)')

    plt.title('Convergence of force mean data')
    plt.xlabel('Time')
    plt.ylabel('Force')
    plt.legend()
    plt.show()    
        
def CheckConvergencePD(fname=None,S=100):
    fdArrays=loadForcePD(fname,S)
    fig1 = plt.figure()#plots fd convergences
    plt.plot(fdArrays[0],fdArrays[1],label='fdx')
    plt.plot(fdArrays[0],fdArrays[4],label='fdx (RA)')
    plt.plot(fdArrays[0],fdArrays[7],label='fdx (M)')
    plt.plot(fdArrays[0],fdArrays[2],label='fdy')
    plt.plot(fdArrays[0],fdArrays[5],label='fdy (RA)')
    plt.plot(fdArrays[0],fdArrays[8],label='fdy (M)')

    plt.title('Convergence of force mean data')
    plt.xlabel('Time')
    plt.ylabel('Force')
    plt.legend()
    plt.show()    

def RunTimeCheck(CANumber,fname=None,NArray=None,MakeNGraphs=False):
    if fname==None:
        fname='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CA_VParam/CylinderAdjust'
        print('no filename given')
    fig6=plt.figure()
    FinalETimes=[]
    FinalCTimes=[]
    FinalSTimes=[]
    for i in range(len(CANumber)):#parses time files and finds final time for each
        ind=CANumber[i]
        thisFile=fname+ind+'/'
        TimeArrays=ParseTime(thisFile)
        FinalETimes.append(np.max(TimeArrays[0]))
        FinalCTimes.append(np.max(TimeArrays[1]))     
        FinalSTimes.append(np.max(TimeArrays[2]))
    if MakeNGraphs:
        if NArray==None:
            print('NEquiv Array not included')
        else:
            plt.plot(NArray,FinalETimes,'--o', linewidth=0.5)
            plt.xlabel('NEquivelant (arb.)')
            plt.ylabel('Execution Time (s)')
            plt.title('Runtime for simtime '+str(np.mean(FinalSTimes))+" s")
            plt.show()
    returnArray=[FinalETimes,FinalCTimes,FinalSTimes]
    return returnArray
        
def AvgFile(TargArr,N=1):
    MeanArr=[]
    for i in range(len(TargArr)):
        if(i<N):
            MeanArr.append(np.mean(TargArr[0:i+N]))
                
        elif(i>(len(TargArr)-N)):
            MeanArr.append(np.mean(TargArr[i-N:len(TargArr)-1]))
        else:
            MeanArr.append(np.mean(TargArr[i-N:i+N]))
    return MeanArr
    



'''
Test files for checking methods work

'''
fileNameV1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CA_VParam/CylinderAdjust'
fileNameNE1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CA_NEquivParam/CylinderAdjustNEquiv'
fileNameT1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CA_TParam/CylinderAdjustT'
fileNameTest='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/DemoFile/CA_DemoFile/'

fileNameSX='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CA_SXParam/CylinderAdjustSX'
filenameTimeTest=fileNameV1+'5/'
FileNameLat2='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/SSOWeekendRun/CA_SSOLatDep_L270H500_L'
FileNameLat3='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/ToRunCases/SSOReRunNative/CA_SSOLatDep_L270H500_L'

FileNameLat1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/JakubSSOFiles/CA_SSOLatDep_L270H500_L'
FilePD1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/NewCanonAttempt/cylinderCaseWorking/'
FilePDHO1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/CylinderHOCase/'
FilePDO1='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/CylinderOCase/'
FilePDO2='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/CylinderHOCase2/'
FilePDO3='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/CylinderHCase/'
FilePDO4='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/CylinderHOCase/'

Lats1=[0,30,60,90,300,330]
LatS1=[str(i) for i in Lats1]
speeds=[7.62e+03, 7.62e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.62e+03]#for first half of SSO run
RHN2=[5.626e+10, 2.943e+10, 4.134e+10, 4.089e+10, 4.011e+10, 6.772e+10]
RHO=[4.597e+12, 4.643e+12, 5.894e+12, 5.715e+12, 6.137e+12, 6.203e+12]
RHN=[5.021e+10, 8.415e+10, 1.743e+11, 1.626e+11, 1.343e+11, 9.675e+10]
RHH=[1.378e+11, 1.554e+11, 1.730e+11, 1.563e+11, 1.519e+11, 1.358e+11]
RHe=[1.591e+12, 2.129e+12, 1.771e+12, 1.590e+12, 1.678e+12, 1.314e+12]

RhoTot=[]
for i in range(len(RHN2)):
    RhoTot.append(RHN2[i]+RHO[i]+RHN[i]+RHe[i]+RHH[i])
    # print('Lat: '+str(Lats1[i])+" RhoN: "+str(RhoTot[i]))

def CheckLat(Lats,LatNames,fname=FileNameLat1):
    PFAnalyse(Lats, LatNames, 'Latitude', 'deg.',False, True, fname)
    

# PFAnalyse(Lats2, LatS2,'Latitude','deg',False,True,FileNameLat3)

##SSO Weekend run stats
Lats2=[0,30,60,90,120,150,180,210,240,270,300,330]
LatS2=[str(i) for i in Lats2]
fdx=[3.8907e-08, 4.0202e-08, 4.9561e-08, 4.7816e-08, 5.0953e-08, 5.0904e-08, 4.7978e-08, 2.9922e-08, 1.9036e-08, 1.7705e-08, 2.3929e-08, 3.6904e-08]
fdy=[-2.0527e-09, -2.1628e-09, -2.5812e-09, -2.4662e-09, -2.6309e-09, -2.5991e-09, -2.447e-09, -1.6789e-09, -1.1824e-09, -1.0873e-09, -1.4252e-09, -1.9653e-09]
NEQUIVS=[3.22e+05, 3.52e+05, 4.03e+05, 3.83e+05, 4.07e+05, 3.91e+05, 3.66e+05, 2.93e+05, 2.37e+05, 2.17e+05, 2.68e+05, 3.15e+05]
SPEEDS=[ 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03, 7.63e+03]
RHON2=[5.626e+10, 2.943e+10, 4.134e+10, 4.089e+10, 4.011e+10, 6.772e+10, 1.087e+11, 1.441e+10, 2.827e+09, 2.541e+09, 6.103e+09, 4.961e+10]
RHOO2=[1.160e+09, 4.378e+08, 6.650e+08, 5.885e+08, 6.712e+08, 1.660e+09, 2.417e+09, 1.994e+08, 3.182e+07, 2.535e+07, 8.303e+07, 1.263e+09]
RHON=[5.021e+10, 8.415e+10, 1.743e+11, 1.626e+11, 1.343e+11, 9.675e+10, 7.023e+10, 2.195e+10, 1.365e+10, 1.787e+10, 1.612e+10, 3.233e+10]
RHOO=[4.597e+12, 4.643e+12,5.894e+12, 5.715e+12, 6.137e+12, 6.203e+12, 5.805e+12, 3.317e+12, 1.802e+12, 1.691e+12, 2.438e+12, 4.329e+12]
RHOAr=[6.831e+05, 2.420e+05, 2.756e+05, 3.447e+05, 2.623e+05, 8.987e+05, 1.643e+06, 1.145e+04, 1.882e+04, 3.127e+03, -4.620e+03, 6.891e+05]
RHOHe=[1.591e+12, 2.129e+12, 1.771e+12, 1.590e+12, 1.678e+12, 1.314e+12, 1.153e+12, 2.213e+12, 2.551e+12, 2.282e+12, 2.580e+12, 1.663e+12]
RHOH=[1.378e+11, 1.554e+11, 1.730e+11, 1.563e+11, 1.519e+11, 1.358e+11, 1.783e+11, 2.912e+11, 3.613e+11, 3.453e+11, 3.180e+11, 2.236e+11]
TN=[8.394e+02, 8.064e+02, 8.360e+02, 8.185e+02, 8.322e+02, 8.724e+02, 8.917e+02, 7.542e+02, 6.652e+02, 6.542e+02, 7.050e+02, 8.484e+02]
RhoM=[1.367e-16, 1.411e-16, 1.746e-16, 1.683e-16, 1.794e-16, 1.792e-16, 1.689e-16, 1.045e-16, 6.585e-17, 6.116e-17, 8.309e-17, 1.295e-16]
RhoM2=[i*1e3 for i in RhoM]#RhoM, in kg/m3
FComp=[2*np.abs(i) for i in fdy]
FdTheory=[]
for i in range(len(RhoM2)):
    FdTheory.append((SPEEDS[i]**2)*(RhoM2[i]/2)*1.09*9e-3)
print(FdTheory)

RhoTot=[]
for i in range(len(RHON2)):
    RhoTot.append(RHON2[i]+RHOO[i]+RHON[i]+RHOHe[i]+RHOH[i]+RHOO2[i]+RHOAr[i])
    print('Lat: '+str(Lats2[i])+" RhoN: "+str(RhoTot[i]))


# print(fdx2)
# print(fdy2)
# PFAnalyse(Lats2, LatS2,'Latitude','deg',False,True,FileNameLat3)

def printForceGraphs(IndVar,fx,fy,paramName, units):#prints graphs of fdx/fdy against pre-written variable lists
    fig5,(ax1,ax2)=plt.subplots(1,2)
    IndVar, fx,fy = (list(t) for t in zip(*sorted(zip(IndVar, fx,fy))))
    ax1.plot(IndVar,fx,'-o', color='darkGreen',linewidth=0.5,label='fdX')
    ax2.plot(IndVar,fy,'-o', color='darkOrange',linewidth=0.5,label='fdY')
    ax1.grid(alpha=0.5)
    ax2.grid(alpha=0.5)
    ax1.set_title('Drag force')
    ax2.set_title('Compression forces')
    ax1.set_ylabel('Force(N)')
    ax2.set_ylabel('Force(N)')
    ax1.set_xlabel(paramName+" ("+units+')')
    ax2.set_xlabel(paramName+" ("+units+")")

# printForceGraphs(RhoM2,fdx,FComp,'Mass density','kg/m^3')


# ArrArr1=[NEQUIVS,SPEEDS,RHON2,RHOO2,RHON,RHOO,RHOAr,RHOHe,RHOH,TN,RhoM]
# ArrName1=['NEQUIVS','SPEEDS','RHON2','RHOO2','RHON','RHOO','RHOAr','RHOHe','RHOH','TN','RhoM']

ArrArr1=[RhoTot,TN,RhoM2]
ArrName1=['RhoTot','TN','RhoM']
def plotParams(ArrArr,ArrName, fdx1):
    for i in range(len(ArrArr)):
        # colmap=np.linspace(0,1,len(ArrArr[i]))
        plt.figure(i)
        plt.scatter(ArrArr[i],fdx1,c=Lats2,cmap=RBCmap())
        cbar =plt.colorbar()
        cbar.set_label('Latitude')

        plt.xlabel(str(ArrName[i]))
        plt.title(str(ArrName[i]))
        plt.ylabel('Fdx (N)')
        plt.show()

def RBCmap():#fades from red to blue (in this case)
    cdict = {'red':   [(0.0,  1.0, 1.0),
                   (0.25,  0.5, 0.5),
                   (0.5,  0.0, 0.0),
                    (0.75,  0.5, 0.50),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.25,  0.5, 0.5),
                   (0.5,  1.0, 1.0),
                   (0.75,  0.5, 0.5),
                   (1.0,  0.0, 0.0)]}
    cmapMine=matplotlib.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmapMine

def SSOCmap(): #returns cmap of day-night cycle eventually
     cdict = {'red':   [(0.0,  1.0, 1.0),
                   (0.25,  0.75, 0.75),
                   (0.75,  0.25, 0.250),
                   (1.0,  0.0, 0.0)],

         'green': [(0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.25,  0.25, 0.25),
                   (0.75,  0.75, 0.75),
                   (1.0,  1.0, 1.0)]}
     cmapMine=matplotlib.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
     return cmapMine
   
# plotParams(ArrArr1,ArrName1,fdx)

# linArr=np.linspace(0,10,20)
# plt.scatter(linArr,linArr,c=linArr,cmap=RBCmap())
def theoryForcePlots(IndVar,fx,fC):
    IndVar, fx,fC = (list(t) for t in zip(*sorted(zip(IndVar, fx,fC))))

    plt.figure(5)
    plt.plot(IndVar,fx,'-o',c='red', linewidth=0.5,label='Simulated Drag force')#,cmap='autumn')#'winter')
    plt.plot(IndVar,fC,'-o',c='green', linewidth=0.5,label='Stokes drag force')#Lats2,cmap='winter')
    plt.legend()
    plt.xlabel('Mass density')
    # plt.title('Mass density')
    plt.ylabel('Fdx (N)')
    plt.show()



# theoryForcePlots(RhoM2,fdx,FdTheory)
# CheckLat(Lats1,LatS1)


# SXYN=[1,2,3,4,5,6,7,8,9,10,14,18,30,40,50,60,70,80] #testing squared realtionship for NCElls vs SX/SY
# NCells=[44,152,342,608,950,1368,1862,2432,3078,3800,7448,12312,34200,60800,95000,136000,186000,243200]
# SX2=[i**2 for i in SXYN]
# plt.scatter(SX2,NCells)
# plt.plot(np.linspace(SX2[0],SX2[len(SX2)-1],30),np.linspace(NCells[0],NCells[len(NCells)-1],30),'--')

def ProbeTest(fname,n,times):
    FdArray=[]
    mesh = Ofpp.FoamMesh(fname)
    mesh.read_cell_centres(fname+'1e-05/C')
    PArr=[]
    for entry in times:
        t=str(entry)
        PFile=fname+t+'/p'
        P = Ofpp.parse_internal_field(PFile)
        PArr.append(P)

    meshgridx=[m[0] for m in mesh.cell_centres]
    meshgridy=[m[1] for m in mesh.cell_centres]
    plt.figure(2)
    
    # meshMap=np.linspace(0,1,len(meshgridx))
    plt.scatter(meshgridx,meshgridy)
    for N in n:
        plt.scatter(meshgridx[N],meshgridy[N])
    plt.show()
    print(len(meshgridy))
    
    
    TimeValues=[float(i) for i in times]
    plt.figure(1)
    plt.title("Pressure probes")
    plt.xlabel("Simulation time (s)")
    plt.ylabel("Pressure (Pa)")
    for N in n:
        ValArray=[]
        for i in range(len(PArr)):
            ValArray.append(PArr[i][N])
        plt.scatter(TimeValues,ValArray,label='('+'{:.3e}'.format(meshgridx[N])+", "+'{:.3e}'.format(meshgridy[N])+")")
    plt.legend()       
        

timesArr=['1e-05','2e-05','3e-05','4e-05','5e-05','6e-05','7e-05','8e-05','9e-05','0.0001','0.00011','0.00012','0.00013','0.00014','0.00015','0.00016','0.00017','0.00018','0.00019','0.0002']
# ProbeTest(fileNameTest,[400,500,600,700],timesArr)

def TestLog(fname):#graphs for presentation of log reading tools
    swag=openfoam(fname,'0')
    swag.loadLog()
    print(swag.__dict__.keys)
    plt.figure(3)
    plt.xlabel('Simulation Time (s)')
    plt.ylabel('Force (N)')
    plt.title('Force on Obstacle')
    Fmag=[np.sqrt(swag.FTot[0][i]**2+swag.FTot[1][i]**2+swag.FTot[2][i]**2) for i in range(len(swag.FTot[0]))]
    Pmag=[np.sqrt(swag.FPress[0][i]**2+swag.FPress[1][i]**2+swag.FPress[2][i]**2) for i in range(len(swag.FPress[0]))]
    Vmag=[np.sqrt(swag.FVisc[0][i]**2+swag.FVisc[1][i]**2+swag.FVisc[2][i]**2) for i in range(len(swag.FVisc[0]))]

    plt.plot(swag.STime[1:],Fmag,label='total Force')
    plt.plot(swag.STime[1:],Pmag,label='Pressure Contribution')
    plt.plot(swag.STime[1:],Vmag,label='Viscous Contribution')
    plt.legend()
    plt.show()
    # plt.plot(swag.STime,swag.KEAvg)
    # plt.plot(swag.STime,swag.ETotAvg)

def TheoryForceRatio(IndVar,fx,fC):
    ForceRatio=[]
    plt.figure(6)
    for i in range(len(fx)):
        ForceRatio.append((fx[i]/fC[i]))
    plt.scatter(IndVar,ForceRatio,c='red',label='Force Ratio')#,cmap='autumn')#'winter')
    # plt.legend()
    plt.xlabel('Mass density (kg/m^3)')
    plt.ylabel('Fx/FStokes')
    plt.title("Ratio of simulated force to stokes drag Force")
    plt.show()


# theoryForcePlots(RhoM2,fdx,FdTheory)
# TheoryForceRatio(RhoM2,fdx,FdTheory)

def IonComparisonPlot():
    RhoIon=[]
    RhoNeut=[]
    CDRatio=[]
    for i in range(len(RHOH)):
        # print(len(RHOH))
        # print(i)
        RhoIon.append(RHOH[i]+RHOHe[i]+RHOO[i]+RHON[i]+RHOAr[i])
        RhoNeut.append(RHON2[i]+RHOO2[i])
        CDRatio.append((RHOH[i]+RHOHe[i]+RHOO[i]+RHON[i]+RHOAr[i])/(RHON2[i]+RHOO2[i]))
    print(RhoIon)
    print(Lats2)
    plt.figure(7)

    plt.plot(Lats2,CDRatio)
    plt.legend()
    plt.show()
    
def presentationPlots():
    fdxNano=[]
    for i in range(len(fdx)):
        fdxNano.append(fdx[i]*10e9)
    plt.figure(1)
    plt.plot(Lats2,fdxNano,'--o', linewidth=0.5)
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Force (nN)")
    plt.title("Drag force throughout SSO")
    plt.show()
    plt.figure(2)
    plt.plot(Lats2,RhoM2)
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Mass density (kg/m^3)")
    plt.title("Mass density through SSO")
    plt.show()
    
    
    
# ForcesHO=loadForce(FilePDHO1)
# ForcesO=loadForce(FilePDO1)
# CheckConvergencePD(FilePDHO1)
CheckConvergencePD(FilePDO4)
RunTimeCheck(FilePDO1)

# CheckConvergence(FilePDO1)
# ForceCheck=[]
# ForceCheck.append(ForcesHO[1][0])
# for i in range(len(ForcesHO[1])-1):
#     ForceCheck.append(ForcesHO[1][i+1]-ForcesHO[1][i])
# plt.plot(ForcesHO[0],ForceCheck)    



# presentationPlots()
# TestLog(fileNameTest)
# IonComparisonPlot()
# swag=openfoam(fileNameT1+'1/','0')
# swag.loadLog()
# swag.makeRAvg('pMag',5)
# swag.makeRAvgVec('FTot',30)


# print(swag.__dict__.keys())
# plt.plot(swag.__dict__['STime'],swag.__dict__['FTot'][0])
# plt.plot(swag.__dict__['STime'],swag.__dict__['FTotMean'][0])

# PFAnalyse([1,2,3,4,5,6,7,8,9],['1','2','3','4','5','6','7','8','9'],'SX','arb',False,True,fileNameSX)
# RunTimeCheck(['1','2','3','4','5','6','7','8','9'],fileNameSX,[1,2,3,4,5,6,7,8,9],True)
# CheckConvergence(fileNameSX+'4/',800)

# UFAnalyse([0.5, 1.0,1.736, 2.0, 3.0, 4.0, 5.0],['5','6','11','7','8','9','10'],False,True,fileName)
# ParseTime(filenameTimeTest)
# RunTimeCheck(LatS1,FileNameLat1,RhoTot,True)

# PFAnalyse([0.5, 1.0,1.736, 2.0, 3.0, 4.0, 5.0],['5','6','11','7','8','9','10'],'Velocity','km/s',False,True,fileNameV1)

# PFAnalyse([0.5, 1.0,1.736, 2.0, 3.0, 4.0, 5.0],['5','6','11','7','8','9','10'],'Velocity','km/s',False,True,fileNameV1)
# PFAnalyse([5e3,1e4,5e4,8e4,1e5,5e5],['6','5','1','2','3','4'],"NEquiv",'arb',False,True,fileNameNE1)
# checkArrs=ParseLog(fileNameT1+'1/')
# print(len(checkArrs[1]))
# for i in range(len(checkArrs)):
#     print(str(len(checkArrs[i]))+'   '+str(checkArrs[i][0:10]))

# # fig2 =plt.figure()
# # plt.plot(checkArrs[0],checkArrs[12])



# plt.xlabel('Simulation Time (seconds)')
# plt.show()
# print(checkArrs[1][0:10])




# TimeCheck = ParseTime('/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/CylinderAdjustSymmetryAttempt/')
# fig2=plt.figure()
# plt.plot(TimeCheck[0],TimeCheck[2][:-1])
# plt.show()

# PFAnalyse([0,30,60,90,300,330],['0','30','60','90','300','330'],"Latitude",'arb',False,True,FileNameLat1)

# PFAnalyse([10,1000],['1','2'],"Temp.",'K',False,True,fileNameT1)

