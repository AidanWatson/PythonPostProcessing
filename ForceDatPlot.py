#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 11:31:06 2021

@author: aidan
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import Ofpp
from os.path import basename, realpath, exists
import ForceIntParse as FIP

def ConvergencePlots(fname,S1=200,S2=200):
    '''
    Parameters
    ----------
    fname : File name within 3rdtrypdrun
    S1 : number of points to make mean from
        DESCRIPTION. The default is 200.
    S2 : Number of points to check for standard deviation
        DESCRIPTION. The default is 200.
    Returns
    -------
    Plots forces, means, and standard deviation for a given file
    '''
    FDArr=ReadForceDat(fname,S1)
    CheckConvergencePD(None, FDArr,S1)
    STDPlot(FDArr,S2,0.001)
    
    
FPath='/home/aidan/OpenFOAM/OpenFOAM-v1706/FaunFiles/3rdTryPDRun/'
def IndexSearch(WN, Val):#returns index of a point (Val) in array (WN)
    i=0#index of current closest fit
    j=0
    for data in WN:
        if(abs(data-Val)<=abs(WN[i]-Val)):
            i=j
        j=j+1
    return i

def ReadForceDat(fname=None,S=100):
    """
    Parameters
    ----------
    fname : File Name within 3rdTryPDRun directory
    S: number of poitns to smooth by
    Returns
    -------
    LoadForceArray from FIP file
    """
    fname1=FPath+fname+'/'
    fdArrays=FIP.loadForcePD(fname1,S)
    return fdArrays

def CheckConvergencePD(fname=None,FDat=None,S=100,plotY=True):
    """
    Method for plotting the convergence of force data with time

    Parameters
    ----------
    fname : File name within 3rdTryPDRun directory
    FDat : Alternative way to call: inputs data as an array if already loaded
    S : number of points on either side to smooth with
    plotY: Determines if the y component of force should be plotted alongside x

    Returns
    -------
    plots the x and y components of the force data by timestep, as well as the 
    Runnning average and mean of this data

    """
    if fname!=None:
        fname1=FPath+fname+'/'
        fdArrays=FIP.loadForcePD(fname1,S)
    elif FDat!=None:
        fdArrays=FDat
    else: 
        print("needs a source of force information")
        
    fig1 = plt.figure()#plots fd convergences
    plt.plot(fdArrays[0],fdArrays[1],label='fdx')
    plt.plot(fdArrays[0],fdArrays[4],label='fdx (RA)')
    plt.plot(fdArrays[0],fdArrays[7],label='fdx (M)')
    if plotY:
        plt.plot(fdArrays[0],fdArrays[2],label='fdy')
        plt.plot(fdArrays[0],fdArrays[5],label='fdy (RA)')
        plt.plot(fdArrays[0],fdArrays[8],label='fdy (M)')

    plt.title('Convergence of force mean data')
    plt.xlabel('Time')
    plt.ylabel('Force')
    plt.legend()
    plt.show()  



def STDPlot(fdArr,S=100,ST=0.001):
    '''
    Parameters
    ----------
    fdArr : standard array from loadforce file
    S : Area to smooth, optional
    ST : Time after which standard deviation is measured, also optional.

    Returns
    -------
    Prints standard deviation for last timestep
    plots standard deviation over time

    '''
    STDs=FIP.StdDevCheck(fdArr[0],fdArr[1],S, ST)
    STDArr=STDs[0]
    print("Standard deviation for last "+str(ST)+"s: "+"{:.4e}".format(STDs[1]))
    fig1=plt.figure()
    time=fdArr[0]
    plt.plot(time,STDArr,label='standard deviation of FdX')
    plt.xlabel('Time (s)')
    plt.ylabel("Standard Deviation")
    plt.legend()
    plt.show()
    
        
def CheckConvSTD(fdArr,S=[100,200],ST=0.001): 
    """
    plots standard deviation curves on the same axis (gimicky, not for normal use)
    """
    fig2=plt.figure()
    time=fdArr[0]
    for s in S:
       STDs=FIP.StdDevCheck(fdArr[0],fdArr[1],s, ST)
       plt.plot(time,STDs[0],label='S='+str(s))
    plt.xlabel('Time (s)')
    plt.ylabel("Standard Deviation")
    plt.legend()
    plt.show()

