#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 11:24:01 2021

@author: aidan
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import Ofpp
from os.path import basename, realpath, exists
from scipy import interpolate




REarth=6371000
G=6.673e-11
MEarth=5.972e24

def findInd(Param, LatArr,LatGoal):
    targInd=0
    for i in range(len(LatArr)):
        if(np.abs(LatArr[targInd]-LatGoal)>=np.abs(LatArr[i]-LatGoal)):
            targInd=i
    return targInd


class OrbitParam(object):
    """
    Takes orbital parameters as input and outputs atmospheric parameters
    currently designed for Phi,Theta, Ecc=0, and varying RE
    RE: height above surface in metres (input in km)
    Ecc:eccentricity
    Phi:angle away from being in plane with daynight boundary (check word for this)
    Theta:Angle in plane
    (need to check real words for these, but theta is the angle specifying progress through orbit, phi is the angle of the orbit all orbits pass over poles)
    
    
    """
    
    
    
    def __init__(self,RE=0,Theta=0,Ecc=0,Phi=0):
        if RE==0:
            RE=input('please specify the height above the surface of earth in m:')
        else:
            self.RE=RE*1000
        self.Theta=Theta
        self.Ecc=Ecc
        self.Phi=Phi
        # CalcVel()
        
    
    
    def CalcVel(self, Height=None):
        """
        calculates velocity relative to the earth (currently perpendicular velocity)
        currently uses basic orbital velocity formula
        """
        if Height==None:
            Rtot=self.RE+REarth
        else:
            Rtot=(Height*1000)+REarth
        VelPerp=np.sqrt(G*MEarth)*1/np.sqrt(Rtot)
        self.VelPerp=VelPerp
        # print('perpendicular velocity: '+str(VelPerp)+' (m/s)')
        return VelPerp
    
    def CalcVeLat(self, Lat):
        """
        calculates velocity relative to the earths atmosphere for SSO
        """
        Rtot=self.RE+REarth
        VelPerp=(np.sqrt(G*MEarth)*1/np.sqrt(Rtot))
        VeLat=(2*np.pi*Rtot/86400*np.cos(np.radians(Lat)))
        # print(VelPerp)
        # print(VeLat)
        VeRel=(np.sqrt(VelPerp**2+VeLat**2))
        self.VelPerp=VelPerp
        # print('perpendicular velocity: '+str(VelPerp)+' (m/s)')
        return VeRel
    
    def LoadDataMSIS(self,filename):
        self.__dict__['MSISFname']=filename
        Dat=readMSISFile(filename)
        dataNames=['height','RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
        for i in range(len(Dat)):
            rawName=dataNames[i]
            self.__dict__[rawName]=Dat[i]
    
    
    def getParam(self, param=None,TargHeight=None):#interpolates parameter from array at given height
        if param==None:
            param=input('enter param to be output:')
        height=self.__dict__['height']

        targ=self.__dict__[param]
        tck = interpolate.splrep(height, targ, s=1)
        if TargHeight==None:
            PH=interpolate.splev((self.RE/1000),tck,der=0)
        else:
            PH=interpolate.splev((TargHeight),tck,der=0)

        print('interpolated '+param+' succesfully')
        return PH



    def CalcRhoTot(self):#returns total number density of particles
        RhoNames=['RhoO',"RhoN2","RhoO2","RhoHe","RhoAr","RhoH","RhoN"]
        RhoTot=[]
        for i in range(len(self.__dict__['RhoO'])):
            RhoTot.append(0)
            for nam in RhoNames:
                RhoTot[i]=RhoTot[i]+self.__dict__[nam][i]*1e6
        self.__dict__['RhoTot']=RhoTot
        return RhoTot
    
    def OutputParams(self):
        Vel=self.CalcVel()
        RhO=self.getParam('RhoO')*1e6        
        RhN2=self.getParam('RhoN2')*1e6        
        RhO2=self.getParam('RhoO2')*1e6        
        RhM=self.getParam('RhoMass')*1e3        
        Tn=self.getParam('Tn')        
        RhHe=self.getParam('RhoHe')*1e6
        RhAr=self.getParam('RhoAr')*1e6        
        RhH=self.getParam('RhoH')*1e6        
        RhN=self.getParam('RhoN')*1e6       
        print(" ")
        print('   __________________Ionosphere Parameters_____________________')
        print('At a height of '+str(self.__dict__['RE']/1000)+' km, parameters as extracted from MSIS file '+self.__dict__['MSISFname']+" ")
        print("Mass density: "+str(RhM)+' '+'{:.4e}'.format(RhM)+" (kg/m^3) with number density: "+'{:.3e}'.format(np.sum([RhO,RhN2,RhO2,RhM,RhHe,RhAr,RhH,RhN]))+" (1/m^3)")
        print('Perpendicular velocity: '+'{:.3e}'.format(Vel)+' (m/s) (calculated for spherical orbit only)')
        print('Neutral Temperature: '+'{:.3e}'.format(Tn))
        print('_________ Density of Ions ________')
        print('O:'+'{:.3e}'.format(RhO)+'      N2:'+'{:.3e}'.format(RhN2)+'       O2:'+'{:.3e}'.format(RhO2))
        print('He: '+'{:.3e}'.format(RhHe)+'  Ar:'+'{:.3e}'.format(RhAr)+'  H:'+'{:.3e}'.format(RhH)+'   N:'+'{:.3e}'.format(RhN))
        print('Ion densities in 1/m^3')
        print('mass interpolation currently broken, as is ion number density below 90km for H,N,O ')
    
    def WriteFile(self,fname, Harray): #prints file of rows of variables for bash script to run
        writeSwag=open(fname, "w+")
        VariablesArray=['RhoO', 'RhoN2', 'RhoO2',  'Tn', 'RhoHe', 'RhoAr', 'RhoH', 'RhoN']
        VarString='Heights: '
        for Hval in Harray:
            VarString=VarString+str(Hval)+' '
        writeSwag.write(VarString+' \n')
        print(VarString)
    
        for VarKey in VariablesArray:
            VarString=VarKey+' '
            for Hval in Harray:
                print(VarKey)
                print(Hval)
                VarString=VarString+' '+'{:.3e}'.format(self.getParam(VarKey,Hval)*1e6)
            print(VarString)
            writeSwag.write(VarString+' \n')
        VarString=' RhoNum: '
        for Hval in Harray:
            VarString=VarString+' '+'{:.2e}'.format(self.FindN(Hval,5e3))
        writeSwag.write(VarString+' \n')

    def WriteFileLat(self,fname, LArray): #prints file of rows of variables for bash script to run
        writeSwag=open(fname, "w+")
        VariablesArray=['RhoO', 'RhoN2', 'RhoO2',   'RhoHe', 'RhoAr', 'RhoH', 'RhoN']
        VarString='latitudes: '
        for Hval in LArray:
            VarString=VarString+str((Hval))+' '
        writeSwag.write(VarString+' \n')
        print(VarString)
        VarString='Tn: '
        for Hval in LArray:
            VarString=VarString+' '+'{:.3e}'.format(self.getParam('Tn',Hval))
        writeSwag.write(VarString+' \n')
        
        for VarKey in VariablesArray:
            VarString=VarKey+' '
            for Hval in LArray:
                print(VarKey)
                print(Hval)
                VarString=VarString+' '+'{:.3e}'.format(self.getParam(VarKey,Hval)*1e6)
            print(VarString)
            writeSwag.write(VarString+' \n')
        VarString=' NEquivs: '
        for Hval in LArray:
            VarString=VarString+' '+'{:.2e}'.format(self.FindN(Hval,5e3))
        writeSwag.write(VarString+' \n')
        VarString='RelSpeeds: '
        for LVal in LArray:
            VarString=VarString+' '+'{:.2e}'.format(self.CalcVeLat(LVal))
        writeSwag.write(VarString+' \n')

            
            
            
    
    def FindN(self,Height, NGoal=None):#calculates total number density per cubic metre from component ions, adjusts nequiv to suit
        if NGoal==None:
            NGoal=5e3
        RhO=self.getParam('RhoO',Height)*1e6        
        RhN2=self.getParam('RhoN2',Height)*1e6        
        RhO2=self.getParam('RhoO2',Height)*1e6        
        RhHe=self.getParam('RhoHe',Height)*1e6
        RhAr=self.getParam('RhoAr',Height)*1e6        
        RhH=self.getParam('RhoH',Height)*1e6        
        RhN=self.getParam('RhoN',Height)*1e6  
        totN=np.sum([RhO,RhN2,RhO2,RhHe,RhAr,RhH,RhN])
        
        NEquiv=NGoal*totN/1e11
        return NEquiv
            
        
    




def readMSISFile(filename=None):#reads a height dependent MSIS file
    """
    Reads from MSIS files from https://ccmc.gsfc.nasa.gov/cgi-bin/modelweb/models/vitmo_model.cgi
    
    Parameters
    ----------
    filename : string, optional

    Returns
    -------
    DatArray : TYPE
        contains list of height(or other independent variable), density of O,N2,O2,Mass, temperature, density of He,Ar,H,N.
    density in units of number per cm^3, except mass density (g/cm^3), height in km, and T in K

    """
    if filename==None:
        filename=input('input file to read data from: ')
    try:
        logfile=str(filename)
        infile = open(logfile, 'r')
    except:
        print('no MSIS file in '+filename)  
    if True:
        indVar=[]
        RhoO=[]
        RhoN2=[]
        RhoO2=[]
        RhoMass=[]
        Tn=[]
        RhoHe=[]
        RhoAr=[]
        RhoH=[]
        RhoN=[]
    DatArray=[indVar,RhoO,RhoN2,RhoO2,RhoMass,Tn,RhoHe,RhoAr,RhoH,RhoN]
    
    for line in infile:
        SL1=line.strip().split('  ')#getting around weirdness of the file printing
        #fudging a missing space in the line
        splitline=[SL1[0],SL1[1],SL1[2],SL1[3].split(' ')[0],SL1[3].split(' ')[1],SL1[4],SL1[5],SL1[6],SL1[7],SL1[8]]

        for i in range(len(DatArray)):
            DatArray[i].append(float(splitline[i]))
        
    return DatArray
    

def HInterp(Height=None, TargDat=None,heightFac=None,S=0,d=0):#takes height and target data as input, returns array of new height and interpolated data
   #HeightFac is the fineness factor, how many new points are included between each data point, S is the smoothing (usually 0)
    if(Height==None):
        print("No height data given to HInterp")
    if TargDat==None:
        print("no data to interpolate")

    
    tck = interpolate.splrep(Height, TargDat, s=S)#interpolation
    
    Hnew=np.linspace(Height[0],Height[len(Height)-1],len(Height)*heightFac)#new array to extrapolate to
    ynew = interpolate.splev(Hnew, tck, der=d)
        
    return [Hnew,ynew]
    

    
    
def CombineLatFiles(File1,File2, FileFinalName):
    """
    Parameters
    ----------
    File1 : Lat dep MSISE file, starting from South pole at -90 degrees.
    File2 : Lat dep MSISE file, be ending at north pole at 360 degrees.
    FileFinalName : File to be written to.
    Returns
    -------
    prints a File to FileFinalName, containing the latitude dependent data from files 1 and 2 
    Takes two files running from -90 to 90, and returns one, starting at 0 (north pole), with the data from File1, and continuing around to 360 degrees
    """
    try:
        # print(File1)
        logfile=str(File1)
        infile = open(logfile, 'r')
    except:
        print('no MSIS file in '+File1)  
    LatArr=[]
    DatArr=[]
    L2Arr=[]
    DatArr2=[]
    for line in infile:
        SL1=line.strip().split('  ')
        LatArr.append(SL1[0])
        # print(SL1)
        splitline=[SL1[1],SL1[2],SL1[3].split(' ')[0],SL1[3].split(' ')[1],SL1[4],SL1[5],SL1[6],SL1[7],SL1[8]]
        DatArr.append(splitline)
    infile.close()
    
    
    try:
        logfile=str(File2)
        infile = open(logfile, 'r')
    except:
        print('no MSIS file in '+File2)  
    for line in infile:
        SL1=line.strip().split('  ')
        L2Arr.append(SL1[0])
        splitline=[SL1[1],SL1[2],SL1[3].split(' ')[0],SL1[3].split(' ')[1],SL1[4],SL1[5],SL1[6],SL1[7],SL1[8]]
        DatArr2.append(splitline)
    infile.close()    
    
    
    writeSwag=open(FileFinalName, "w+")
    LatFake=np.linspace(0,180,181)
    for i in range(len(LatFake)):
        splitline=DatArr[len(LatFake)-i-1]
        for j in [0,1,3,4,5,6,7,8]:
            splitline[j]=splitline[j]+' '#FUDGEFUDGEFUDGE #adding the right amount of spaces to mimic the MSISE files
        NewString='  '+str(LatFake[i])+"  "
        for dat in splitline:
            NewString=NewString+dat+' '
        NewString=NewString+'\n'
        writeSwag.write(NewString)
        
    LatFake2=np.linspace(181,360,(len(DatArr2)-1))
    for i in range(len(LatFake2)):
        
        splitline=DatArr2[i+1]
        for j in [0,1,3,4,5,6,7,8]:
            splitline[j]=splitline[j]+' '#FUDGEFUDGEFUDGE #adding the right amount of spaces to mimic the MSISE files
        NewString='  '+str(LatFake2[i])+"  "
        for dat in splitline:
            NewString=NewString+dat+' '
        NewString=NewString+'\n'
        writeSwag.write(NewString)    
    writeSwag.close()
        
     
        
     

    
    







def test1():
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'
    swaggy=OrbitParam(300,0,0,0)
    swaggy.CalcVel()
    print(swaggy.__dict__)
    FDat1=readMSISFile(fname)
    densIndex=[1,2,3,6,7,8,9]#index of the density files
    dataNames=['height (km)','r/$ \rho $ O',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
    datType=['-','--','-','-','-','-','-','-','--','--']
    
    
    fig1=plt.figure()
    for ind in densIndex:
        plt.plot(FDat1[0],FDat1[ind],datType[ind],linewidth=0.8,label=dataNames[ind])
    plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.title("Species density")
    plt.xlabel('Height (km)')
    plt.legend()
    plt.show()
    
    fig2=plt.figure()
    plt.plot(FDat1[0],FDat1[5])
    plt.yscale('log')
    plt.ylabel('T(K)')
    plt.xlabel('Height (km)')
    plt.show()

    fig3=plt.figure()
    plt.plot(FDat1[0],FDat1[4])
    plt.yscale('log')
    plt.ylabel('Mass density(g/cm^3)')
    plt.xlabel('Height (km)')
    plt.show()
    


def test2():#testing interpolation for various values of interpolation parameters
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFile'
    swaggy=OrbitParam(300,0,0,0)
    swaggy.CalcVel()
    print(swaggy.__dict__)
    swaggy.LoadDataMSIS(fname)
    print(swaggy.__dict__.keys())
    height=swaggy.__dict__['height']
    TargName='Tn'
    Targ=swaggy.__dict__[TargName]
    
   
    fig4=plt.figure()
    # HInterp(height,Tn,HFac,0)
    plt.plot(height,Targ,label='Raw Data')
    for Sc in [1,40]:
        for S in [0,1,2,10]:
            intDat=HInterp(height,Targ,Sc,S,0)
            plt.plot(intDat[0], intDat[1],label='Sc,S:'+str(Sc)+','+str(S),linewidth=0.4)
        # for S in [0,1]
    plt.xlim([20,60])
    plt.ylim([150,300])
    plt.xlabel('height') 
    plt.ylabel(TargName)
    # plt.yscale('log')
    plt.legend()
    plt.title('Interpolation test')
   
def test3():#plot of longitude in degrees
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'
    swaggy=OrbitParam(300,0,0,0)
    swaggy.CalcVel()
    print(swaggy.__dict__)
    FDat1=readMSISFile(fname)
    densIndex=[1,2,3,6,7,8,9]#index of the density files
    
    IndVarName='height (m)'
    
    dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
    latt=[(a) for a in FDat1[0]]
    IndVar=latt
    fig1=plt.figure()
    for ind in densIndex:
        plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])
    plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.xlabel(IndVarName)
    plt.title('Ion species by '+IndVarName)
    plt.legend()
    plt.show()
    
    intT=HInterp(IndVar,FDat1[5],1)
    fig2=plt.figure()
    plt.plot(IndVar,FDat1[5])
    plt.plot(intT[0],intT[1])
    plt.yscale('log')
    plt.ylabel('T(K)')
    plt.xlabel(IndVarName)
    plt.title("Temperature by "+IndVarName)
    plt.show()
    
    intM=HInterp(IndVar,FDat1[5],1)

    fig3=plt.figure()
    plt.plot(IndVar,FDat1[4])
    plt.plot(intM[0],intM[1])
    plt.yscale('log')
    plt.ylabel('Mass density(g/cm^3)')
    plt.title('Mass density by '+IndVarName)
    plt.xlabel(IndVarName)
    plt.show()

def test4():#combines two opposing latitude files along same axis
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDep1'
    fname2='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDep2'

    swaggy=OrbitParam(300,0,0,0)
    swaggy.CalcVel()
    print(swaggy.__dict__)
    FDat2=readMSISFile(fname)
    FDat3=readMSISFile(fname2)
    FDat1=[]
    
    for i in range(len(FDat2)):
        CombEntry=[]
        for j in range(len(FDat2[i])):
            CombEntry.append(FDat2[i][j])
        for j in range(len(FDat3[i])):
            CombEntry.append(FDat3[i][len(FDat3[i])-j-1])    
        FDat1.append(CombEntry)

    densIndex=[1,2,3,6,7,8,9]#index of the density files
    
    IndVarName='latitude(deg)'
    
    dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
    latt=[(a+90) for a in FDat1[0]]
    IndVar=latt
    fig1=plt.figure()
    for ind in densIndex:
        plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])
    plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.xlabel(IndVarName)
    plt.title('Ion species by '+IndVarName)
    plt.legend()
    plt.show()
    
    fig2=plt.figure()
    plt.plot(IndVar,FDat1[5])
    # plt.yscale('log')
    plt.ylabel('T(K)')
    plt.xlabel(IndVarName)
    plt.title("Temperature by "+IndVarName)
    plt.show()

    fig3=plt.figure()
    plt.plot(IndVar,FDat1[4])
    # plt.yscale('log')
    plt.ylabel('Mass density(g/cm^3)')
    plt.title('Mass density by '+IndVarName)
    plt.xlabel(IndVarName)
    plt.show()

def HTest(height):#tests output of parameters
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'
    swaggy=OrbitParam(height,0,0,0)
    swaggy.LoadDataMSIS(fname)
    print(swaggy.__dict__.keys())
    swaggy.CalcVel()
    swaggy.OutputParams()

def HTest2(heights):#outputs velocity for array of input heights
    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'
    swaggy=OrbitParam(400,0,0,0)
    swaggy.LoadDataMSIS(fname)
    print(swaggy.__dict__.keys())
    for h in heights:
        print('height: '+str(h)+' V:'+str(swaggy.CalcVel(h)))
    swaggy.CalcVel()
    swaggy.OutputParams()

def test6():#checking interpolations of values for low numbers

    fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDepCombi_Long270H500'
    swaggy=OrbitParam(300,0,0,0)
    swaggy.CalcVel()
    print(swaggy.__dict__)
    FDat1=readMSISFile(fname)
    densIndex=[1,2,3,6,7,8,9]#index of the density files
    
    IndVarName='height (m)'
    
    dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
    latt=[(a) for a in FDat1[0]]
    IndVar=latt
    fig1=plt.figure()
    for ind in densIndex:
        intIon=HInterp(IndVar,FDat1[ind],2)
        plt.plot(intIon[0],intIon[1],'--')
        plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])
    plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.xlabel(IndVarName)
    plt.title('Ion species by '+IndVarName)
    plt.legend()
    plt.show()
    
    intT=HInterp(IndVar,FDat1[5],2)
    fig2=plt.figure()
    plt.plot(IndVar,FDat1[5])
    plt.plot(intT[0],intT[1])
    plt.yscale('log')
    plt.ylabel('T(K)')
    plt.xlabel(IndVarName)
    plt.title("Temperature by "+IndVarName)
    plt.show()
    
    intM=HInterp(IndVar,FDat1[5],1)
    fig3=plt.figure()
    plt.plot(IndVar,FDat1[4])
    plt.plot(intM[0],intM[1])
    plt.yscale('log')
    plt.ylabel('Mass density(g/cm^3)')
    plt.title('Mass density by '+IndVarName)
    plt.xlabel(IndVarName)
    plt.show()
    
def FilePrint(HArr,fileIn=None,fileOut=None):#tests output of parameters
    if fileIn==None:
        
        fname='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'
    else:
        fname=fileIn
    swaggy=OrbitParam(HArr[0],0,0,0)
    swaggy.LoadDataMSIS(fname)
    print(swaggy.__dict__.keys())
    swaggy.CalcVel()
    swaggy.OutputParams()
    if fileOut==None:
        swaggy.WriteFile('/home/aidan/Documents/IonSpeciesFiles/FileWriteTest', HArr)
    else:
        swaggy.WriteFile(fileOut, HArr)       


def SSOOrbitTest(fname,Height,Lats):#takes input of file, height, and latitudes, and writes to a file the appropriate parameters
    swaggy=OrbitParam(Height,0,0,0)
    swaggy.LoadDataMSIS(fname)
    swaggy.CalcRhoTot()
    # plt.plot(swaggy.height,swaggy.RhoTot)
    # plt.title("Number Density SSO")
    # plt.xlabel('')
    swaggy.WriteFileLat('/home/aidan/Documents/IonSpeciesFiles/ParamFileLatDepSSO2',Lats)  
    # for L in Lats:
    #     print("Latitude: "+str(L)+"  Tn:"+str(swaggy.getParam("Tn",L)))

def LatCombiTest(FLat1,FLat2,FLat3,RName=' ',Bool=False):#combines two lat dependent files.

       CombineLatFiles(FLat1,FLat2,FLat3)
       # readMSISFile(FLat3)
       
       fname=FLat3
       swaggy=OrbitParam(500,0,0,0)
       swaggy.CalcVel()
       print(swaggy.__dict__)
       FDat1=readMSISFile(fname)
       densIndex=[1,2,3,6,7,8,9]#index of the density files
        
       IndVarName='Latitude (degrees)'
       
       dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
       
       if Bool:
           FDat=[]
           l1=np.linspace(0,180,181)
           l2=np.linspace(180,1,180)
           
           for i in range(len(l1)):
               FDat.append(l1[i])
           for i in range(len(l2)):
               FDat.append(l2[i])
               
           latt=[(a) for a in FDat]
       else:
           latt=[(a) for a in FDat1[0]]

       IndVar=latt
       fig1=plt.figure()
       for ind in densIndex:
           plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])
    
       plt.yscale('log')
       plt.ylabel('Density (1/cm^3)')
       plt.xlabel(IndVarName)
       plt.title('Ion species by '+IndVarName+' '+RName)
       plt.legend()
       plt.show()
        
       intT=HInterp(IndVar,FDat1[5],1)
       fig2=plt.figure()
       plt.plot(IndVar,FDat1[5])
       plt.plot(intT[0],intT[1])
       plt.yscale('log')
       plt.ylabel('T(K)')
       plt.xlabel(IndVarName)
       plt.title("Temperature by "+IndVarName+' '+RName)
       plt.show()
        
       intM=HInterp(IndVar,FDat1[5],1)
    
       fig3=plt.figure()
       plt.plot(IndVar,FDat1[4])
       # plt.plot(intM[0],intM[1])
       plt.yscale('log')
       plt.ylabel('Mass density(g/cm^3)')
       plt.title('Mass density by '+IndVarName+' '+RName)
       plt.xlabel(IndVarName)
       plt.show()

def LatCompare(FLat1,FLat2,RName=' ', Bool=False):#takes two lat dependent files, circling an area, and plots the difference.
    FDat2=readMSISFile(FLat1)
    FDat3=readMSISFile(FLat2)
    if Bool:
        FDat=[]
        l1=np.linspace(0,180,181)
        l2=np.linspace(180,1,180)
           
        for i in range(len(l1)):
            FDat.append(l1[i])
        for i in range(len(l2)):
            FDat.append(l2[i])
               
        latt=[(a) for a in FDat]
    else:
        latt=[(a) for a in FDat2[0]]    
    FDat1=[]
    for i in range(len(FDat2)):
        tempF=[]
        for j in range(len(FDat2[i])):
            tempF.append(FDat2[i][j]-FDat3[i][j])
        FDat1.append(tempF)
    
    
    densIndex=[1,2,3,6,7,8,9]#index of the density files

    IndVarName='Latitude (degrees)'   
    dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]

    IndVar=latt
    fig1=plt.figure()
    for ind in densIndex:
        plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])

    # plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.xlabel(IndVarName)
    plt.title('Ion species Difference by '+IndVarName+' '+RName)
    plt.legend()
    plt.show()
        
    intT=HInterp(IndVar,FDat1[5],1)
    fig2=plt.figure()
    plt.plot(IndVar,FDat1[5])
    plt.plot(intT[0],intT[1])
    # plt.yscale('log')
    plt.ylabel('T(K)')
    plt.xlabel(IndVarName)
    plt.title("Temperature Difference by "+IndVarName+' '+RName)
    plt.show()
        
    intM=HInterp(IndVar,FDat1[5],1)
    
    fig3=plt.figure()
    plt.plot(IndVar,FDat1[4])
       # plt.plot(intM[0],intM[1])
    # plt.yscale('log')
    plt.ylabel('Mass density(g/cm^3)')
    plt.title('Mass density Difference by '+IndVarName+' '+RName)
    plt.xlabel(IndVarName)
    plt.show()        
        
            
def test7(fname,Height,Lats): #tests output of mass density file
    swaggy=OrbitParam(Height,0,0,0)
    swaggy.LoadDataMSIS(fname)
    swaggy.CalcRhoTot()
    MassInd=[]
    for lat in Lats:
        MassInd.append(swaggy.__dict__['RhoMass'][findInd(swaggy.__dict__['RhoMass'],swaggy.height,lat)])
    print(MassInd)
    plt.plot(swaggy.height,swaggy.RhoMass)
    
def LatPlot(FLat1,Indices=[0]):#plots parameters by latitude from file
    FDat1=readMSISFile(FLat1)
    latt=[(a) for a in FDat1[0]]
    densIndex=[1,2,3,6,7,8,9]#index of the density files

    IndVarName='Angle(degrees)'#'Latitude (degrees)'   
    dataNames=[IndVarName,'RhoO',"RhoN2","RhoO2","RhoMass","Tn","RhoHe","RhoAr","RhoH","RhoN"]
    IonRatio=[]
    for i in range(len(FDat1[1])):
        IonRatio.append((FDat1[1][i]+FDat1[6][i]+FDat1[7][i]+FDat1[8][i]+FDat1[9][i])/(FDat1[2][i]+FDat1[3][i]))
    IndVar=latt
    # fig1=plt.figure()
    fig1, axs1=plt.subplots()
    axs1.set_xticks([0, 30,60,90,120,150,180,210,240,270,300,330,360], minor=False)  
    axs1.xaxis.grid(True, which='major',linestyle='--',alpha=0.5)
    for ind in densIndex:
        plt.plot(IndVar,FDat1[ind],linewidth=0.8,label=dataNames[ind])

    plt.yscale('log')
    plt.ylabel('Density (1/cm^3)')
    plt.xlabel(IndVarName)
    plt.title('Ion species by Latitude')
    plt.legend()
    plt.show()
    
    fig4,axs=plt.subplots()
    plt.plot(IndVar,IonRatio)
    # plt.yscale('log')
    axs.set_xticks([0, 30,60,90,120,150,180,210,240,270,300,330,360], minor=False)  
    axs.xaxis.grid(True, which='major',linestyle='--',alpha=0.5)
    plt.ylabel('Ni/Nn')
    plt.xlabel(IndVarName)
    plt.title("Relative Ion Density by Latitude")
    plt.show()   
    
    intT=HInterp(IndVar,FDat1[5],1)
    # fig2=plt.figure()
    fig2,axs=plt.subplots()
    plt.plot(IndVar,FDat1[5])
    plt.plot(intT[0],intT[1])
    # plt.yscale('log')
    axs.set_xticks([0, 30,60,90,120,150,180,210,240,270,300,330,360], minor=False)  
    axs.xaxis.grid(True, which='major',linestyle='--',alpha=0.5)
    plt.ylabel('T(K)')
    plt.xlabel(IndVarName)
    plt.title("Temperature by Latitude")
    plt.show()
        
    intM=HInterp(IndVar,FDat1[5],1)
    
    # fig3=plt.figure()
    fig3,ax=plt.subplots()
    plt.plot(IndVar,FDat1[4])
    ax.set_xticks([0, 30,60,90,120,150,180,210,240,270,300,330,360], minor=False)  
    ax.xaxis.grid(True, which='major',linestyle='--',alpha=0.5)
    plt.ylabel('Mass density(g/cm^3)')
    plt.title('Mass density Through SSO')
    plt.xlabel(IndVarName)
    plt.show()        
        
          
# test1()
# test2()
# test3()
# test4()

# HTest(400)
# test6()
# FilePrint([300,320,340,360,380,400])
fname1='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDep_Long270_H500'
FLat1='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDep_Long90_H500'
FLat2='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDep_Long270_H500'
FLat3='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesLatDepCombi_Long270H500'

FLat4='/home/aidan/Documents/IonSpeciesFiles/ISLatDep_L0H500'
FLat5='/home/aidan/Documents/IonSpeciesFiles/ISLatDep_L180H500'
FLat6='/home/aidan/Documents/IonSpeciesFiles/ISLatDep_Combi_L0H500'
FH1='/home/aidan/Documents/IonSpeciesFiles/IonSpeciesFileFine'

# HTest2([100,500])
# SSOOrbitTest(FLat3,500,[270,240,210,180,120,150])
# LatCombiTest(FLat1,FLat2,FLat3, 'SSO')
# LatCombiTest(FLat4,FLat5,FLat6, "0-180")
# LatCompare(FLat3,FLat6)
# test6()
# test7(FLat3,500,[0,30,60,90,120,150,180,210,240,270,300,330])
# SSOOrbitTest(FLat3,500,[0,30,60,90,120,150,180,210,240,270,300,330])
LatPlot(FLat3,'Angle')

# test1()














            