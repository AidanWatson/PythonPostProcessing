#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:15:01 2021

@author: aidan
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import Ofpp
from os.path import basename, realpath, exists
import ForceIntParse as FIP
import ForceDatPlot as FDP
from scipy.fft import fft, fftfreq

FArr=FDP.ReadForceDat('CylinderHOCase3',300)
# FDP.STDPlot(FArr,100)
STDArr=FIP.StdDevCheck(FArr[0],FArr[1],100)
# fig1=plt.figure(1)
# plt.plot(FArr[0],STDArr[0])
# plt.show()
# FDP.CheckConvSTD(FArr,[100,200,500,1000,3000,5000])

FDP.ConvergencePlots('CylinderHOCase3',300,300)
# F=4.87353826587844e-10
# print("s: "+"{:.2e}".format(F))

# STD=FIP.StdDevCheck(FArr[0],FArr[1],100)
# print(len(STD[0]))
# print(STD[1])
# print(type(FArr[0]))
# print(type(FArr[1][1:5]))



def FourierBS(time=FArr[0],Arr=STDArr[0]):#attempts to deconstruct signal (not working at present)
    N=len(time)#number of timesteps
    dT=time[1]-time[0]#sample spacing
    Arr1=np.array(Arr)
    yf=2/N*np.abs(fft(Arr1))#[:N//2]) #outputs frequency coefficients (only second half give positive frequencies)
    xf=fftfreq(N, dT)#[:N//2] #outputs the frequencies being analysed
    print(len(yf))
    print(yf)
    fig1=plt.figure(2)
    plt.plot(xf,yf)
    plt.show()
    reconstructBS(time,yf,xf)
    
def reconstructBS(time,Coeff,freqs):#attempts to reconstruct the signal
    N=len(time)#number of timesteps
    dT=time[1]-time[0]#sample spacing
    signal=np.zeros(N)
    for j in range(0,1500):#len(Coeff)):
        if j%100==0:
            print(j)
        for i in range(len(signal)):
            signal[i]=signal[i]+Coeff[j]*np.sin(2*np.pi*time[i]*freqs[j]/N)
    fig3=plt.figure(3)
    plt.plot(time,signal)
    plt.show()

        
FourierBS() 

