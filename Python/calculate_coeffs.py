import matplotlib.pyplot as plt
import numpy as np
import ParallelFilterDesign 
import time

Fs = 44100
fname = "Engl"

data = np.genfromtxt(fname+'.csv', delimiter=';')
w = np.zeros(len(data))
tf = np.zeros(len(data))

for i in range(len(data)):
    w[i] = data[i][0]
    tf[i] = data[i][1]

H = ParallelFilterDesign.MinPhaseN(tf,w)

#Dual log
nPoles1 = 9
nPoles2 = 11
pLog1 = np.logspace(np.log10(10),np.log10(500),nPoles1)
pLog2 = np.logspace(np.log10(780),np.log10(20000),nPoles2)
pLog = np.concatenate((pLog1,pLog2))
p = ParallelFilterDesign.freqpoles(pLog,Fs)

(Bm, Am, FIR) =  ParallelFilterDesign.parFiltDesFr(w,H,p,1)

np.savetxt(fname+"_pyt_Bm_dl.csv",Bm,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_Am_dl.csv",Am,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_FIR_dl.csv",FIR, delimiter=";", fmt="%.15f")


#Single warp
nPoles = 40
lambda1 = 0.92
p = ParallelFilterDesign.warpPolesFr(w,1,H,nPoles,lambda1,5)

(Bm, Am, FIR) =  ParallelFilterDesign.parFiltDesFr(w,H,p,1)

np.savetxt(fname+"_pyt_Bm_w.csv",Bm,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_Am_w.csv",Am,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_FIR_w.csv",FIR, delimiter=";", fmt="%.15f")


#Dual warp
CrossFreq = 500.0
WindowCrossoverLength = 50
lambda1 = 0.986
lambda2 = 0.65
nPoles1 = 18
nPoles2 = 22

C = int( np.argmax(w>2*np.pi*CrossFreq/Fs))


# a.shape[0] = rows, a.shape[1] = cols

p = ParallelFilterDesign.DualWarpPolesFr(w,1,tf,C,WindowCrossoverLength,lambda1,lambda2,nPoles1,nPoles2,5)

(Bm, Am, FIR) =  ParallelFilterDesign.parFiltDesFr(w,H,p,1)

np.savetxt(fname+"_pyt_Bm_dw.csv",Bm,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_Am_dw.csv",Am,delimiter=";", fmt="%.15f")
np.savetxt(fname+"_pyt_FIR_dw.csv",FIR, delimiter=";", fmt="%.15f")
