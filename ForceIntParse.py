#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:19:49 2021

Methods for loading force data from a force.dat sheet for python
Plotting methods stored in ForceDatPlot python file




@author: aidan
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




     
def AvgFile(TargArr,N=1):
    """
    Takes running average, smoothing dataset (smoothing properly, not like openfoams method)

    Parameters
    ----------
    TargArr : Array to be averaged
    N : Number of points on either side to include in average

    Returns
    -------
    MeanArr : running mean (smoothed dataset)

    """
    MeanArr=[]
    for i in range(len(TargArr)):
        if(i<N):
            MeanArr.append(np.mean(TargArr[0:i+N]))
                
        elif(i>(len(TargArr)-N)):
            MeanArr.append(np.mean(TargArr[i-N:len(TargArr)-1]))
        else:
            MeanArr.append(np.mean(TargArr[i-N:i+N]))
    return MeanArr


def loadForcePD(fname,S=600): #parses .dat file of force from pdfoam format into separate arrays for force and time
    """
    Takes a file of .dat values from pdfoam 
    Returns array of times, xyz force values, xyz Running average, and OF method of xyz means
    """
    
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
        fdArrays=loadForcePD(fname+ind+'/')
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
        fdArrays=loadForcePD(fname+ind+'/')
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

def StdDevCheck(time,fdArray,S=100,ST=0.001):
    """
    Parameters
    ----------
    time : Array of time values
    fdArray : Array of force values to be analysed (only raw values, not mean)

    Returns
    -------
    Array of standard deviation over time for a range of S points
    Value of standard deviation for the time after time ST
    """
    SDvArr=[]
    # print(len(time))
    for i in range(len(fdArray)):
        if i<=S:
            fd2=fdArray[0:i+S]
            SDvArr.append(np.std(fd2))
        elif i>=(len(fdArray)-S):
            fd2=fdArray[i-S:]
            SDvArr.append(np.std(fd2))
        else:
            fd2=fdArray[i-S:i+S]
            SDvArr.append(np.std(fd2))
    
    
    TInd=IndexSearch(time,ST)
    STDEV=np.std(fdArray[TInd:])
    return ([SDvArr,STDEV])
