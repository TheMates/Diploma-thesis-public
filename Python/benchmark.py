import matplotlib.pyplot as plt
import numpy as np
import ParallelFilterDesign 
import time

fname = "Engl"


CrossFreq = 500.0
WindowCrossoverLength = 50
lambda1 = 0.986
lambda2 = 0.65
nPoles1 = 18
nPoles2 = 22
Fs = 44100

data = np.genfromtxt(fname+'.csv', delimiter=';')
w = np.zeros(len(data))
tf = np.zeros(len(data))

for i in range(len(data)):
    w[i] = data[i][0]
    tf[i] = data[i][1]


C = int( np.argmax(w>2*np.pi*CrossFreq/Fs))


# a.shape[0] = rows, a.shape[1] = cols
H = ParallelFilterDesign.MinPhaseN(tf,w)

ITERS = 100
start = time.time()
for k in range(ITERS):

    #poslat tf
    p = ParallelFilterDesign.DualWarpPolesFr(w,1,tf,C,WindowCrossoverLength,lambda1,lambda2,nPoles1,nPoles2,5)

    (Bm, Am, FIR) =  ParallelFilterDesign.parFiltDesFr(w,H,p,1)

end = time.time()

t = end-start

print("Whole time: %.3f ms \nOne iter time: %.3f ms" % (t*1000, t*1000/ITERS))


