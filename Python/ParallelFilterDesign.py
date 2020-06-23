import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline as CubSpline
from scipy.signal import hilbert
from scipy.linalg import solve as scipySolve
from numpy import pi as pi


def MinPhaseN(Magn, Wspec, npoints = 1024):
    #Magn - double vector, Wspec - double vector
    LinW = np.linspace(0,np.pi,npoints+1)
    lambda1 = -0.9
    lambda2 = lambda1*lambda1
    wLinW = np.arctan2( (1-lambda2)*np.sin(LinW),(1+lambda2)* np.cos(LinW) -2*lambda1)

    #inter = interp1d(Wspec,Magn,kind='cubic', fill_value='extrapolate')     #chová se jinak než matlab 
    inter1 = CubSpline(Wspec,Magn,bc_type='not-a-knot', extrapolate=True)

    #LinMagn = np.abs(inter(wLinW))
    LinMagn = np.abs(inter1(wLinW))

    N = len(LinMagn)-1
    
    PerLinMagn = np.concatenate((LinMagn[:-1],LinMagn[:0:-1] ))        
    PerLinPhase = -np.imag(hilbert(np.log(PerLinMagn)))             #tady už mám taky jiný výsledky než matlab, takže jsou to opět funkce spline a hilbert
    LinPhase = PerLinPhase[:N+1]

    #phSpl = interp1d(wLinW,LinPhase,kind='cubic', fill_value='extrapolate')
    phSpl = CubSpline(wLinW,LinPhase,bc_type='not-a-knot', extrapolate=True)
    Phase = phSpl(Wspec)

    H = Magn * np.exp(Phase*1j)

    return H


def DualWarpPolesFr(w, X, Y, CrossFr, CrossFade, lambda1, lambda2, nPoles1, nPoles2, ITER):
    #Y musí být abs() takže tam nemusím posílat to minphase verzi, je to zbytečná práce zde, ale minphase verzi budu potřebovat v parfiltdes
    if X == 1:
        X = np.ones(len(Y))

    HL = len(X)
    C = CrossFr
    WL = 2*round(CrossFade/2)
    window = np.hanning(2*WL+1)
    window = window[:WL]
    
    whf = np.concatenate((np.zeros(int(C-WL/2)) , window, np.ones(int(HL-C-WL/2))))
    wlf = 1-whf

    X1 = X*wlf + X[C]*whf
    X2 = X*whf + X[C]*wlf
    Y1 = Y*wlf + Y[C]*whf
    Y2 = Y*whf + Y[C]*wlf

    Y1 = MinPhaseN(Y1,w)
    Y2 = MinPhaseN(Y2,w)

    pwarp1 = warpPolesFr(w,X1,Y1,nPoles1,lambda1,ITER)
    pwarp2 = warpPolesFr(w,X2,Y2,nPoles2,lambda2,ITER)

    return np.concatenate((pwarp1,pwarp2))


def warpPolesFr(w,X,Y,PNUM,lambda1,ITER):
    #X real, Y complex
    lambda2 = lambda1*lambda1
    
    w_wp = np.arctan2( (1-lambda2)*np.sin(w),(1+lambda2)* np.cos(w) -2*lambda1)

    Awp = lsidFr(w_wp,X,Y,PNUM,PNUM,ITER)
    pwp = np.roots(Awp)
    p = (pwp+lambda1)/(1+lambda1*pwp)
    return p

def lsidFr(w,X,Y,nB, nA, ITER):
    L = len(Y)
    Z = np.exp(w*-1j)
    PARL = nA+nB+1

    M = np.zeros((PARL,L),complex)

    for k in range(0,nB+1):
        M[k,:] = X*Z**k
    for k in range(0,nA):
        M[k+nB+1,:] = -Y*(Z**(k+1))

    # M = np.transpose(np.conjugate(M)) # udělá jen conjugate, pak je třeba udělat 
    Am = (M.dot(np.transpose(np.conjugate(M)))).real
    b = (np.conjugate(M).dot(Y)).real
    par = np.linalg.solve(Am ,b)
    #par = scipySolve(Am,b, overwrite_a = True, overwrite_b = True, check_finite = False)

    #B = par[0:nB+1]
    A = np.concatenate(([1] ,par[nB+1:]))

    MW = np.zeros((PARL,L),complex)
    YW = np.zeros(L,complex)        #optimalizace? beztak že ne
    for k in range(ITER):
        YW.fill(0)
        AWt = abs(np.polyval(np.concatenate(([1] ,A)),Z))**-1
        for i in range(PARL):
            MW[i,:] = M[i,:]*AWt
        YW = Y*AWt

        Am = (MW.dot(np.transpose(np.conjugate(MW)))).real
        b = (np.conjugate(MW).dot(YW)).real
        par = np.linalg.solve(Am ,b)
        #B = par[0:nB+1]
        A = np.concatenate(([1] ,par[nB+1:]))

    return A

def parFiltDesFr(W, H, p, NFIR):
    # W - double vector, H - complex vector - minphase version, p - complex vector - poles, NFIR - int scalar
    for pole in p:
        if abs(pole)>1. :
            pole = 1./np.conjugate(pole)            #flip poles inside unit circle
        if abs(pole)>0.9995 and np.imag(pole) == 0:
            pole = 0.995                            #if its DC, give it less gain, it tends to feedback vigorously
    
    p = cplxpair(p)
    nPoles = len(p)
    nPolesEven = int(2*np.floor(nPoles/2.))
    ODD = 0
    if nPoles>nPolesEven:
        ODD = 1
    
    Z = np.exp(W*-1j)   
    Z2 = Z**2
    M = np.zeros((len(H),nPolesEven+NFIR),complex)

    for k in range(0,nPolesEven,2):
        A = np.poly(p[k:k+2])
        resp = 1./(A[0] + A[1]*Z + A[2]*Z2)
        M[:,k] = resp
        M[:,k+1] = Z*resp

    if ODD:
        A = np.poly(p[-1])
        resp = 1./(A[0] + A[1]*Z)
        M[:,nPoles-1] = resp
    
    for k in range(NFIR):
        M[:,nPoles+k] = Z**(k)
    
    A = (np.transpose(np.conjugate(M)).dot(M)).real
    b = (np.transpose(np.conjugate(M)).dot(H)).real
    par = np.linalg.solve(A ,b)
    #par = scipySolve(A,b, overwrite_a = True, overwrite_b = True, check_finite = False)

    #Constructing Bm, Am, matrices
    Am = np.zeros((3, int(np.ceil(nPoles/2))))      #poles 
    Bm = np.zeros((2, int(np.ceil(nPoles/2))))      #zeros

    for k in range(int(nPolesEven/2)):
        Am[:,k] = np.poly(p[2*k:2*k+2])
        Bm[:,k] = par[2*k:2*k+2]

    if ODD:
        Am[:,-1] = np.concatenate((np.poly(p[nPoles]),[0]))
        Bm[:,-1] = np.concatenate((par[nPoles],[0]))

    FIR = 0
    if NFIR>0:
        FIR = par[-NFIR:]

    return (Bm,Am,FIR)

def cplxpair(p):
    p = np.sort_complex(p)
    reals = p[np.imag(p)==0]
    p = np.delete(p,np.where(np.imag(p)==0))
    return np.concatenate((p,reals))

def freqpoles(fr,Fs= 44100,Q = -1):
    wp = 2*pi*fr/Fs
    if Q <0:
        wp = np.sort(wp)
        pnum = len(wp)
        dwp = np.zeros(pnum)
        for k in range(1,pnum-1):
            dwp[k] = (wp[k+1]-wp[k-1])/2.
        dwp[0] = wp[1] - wp[0]
        dwp[pnum-1] = wp[pnum-1]-wp[pnum-2]
    if Q>0:
        dwp = wp/Q
    p = np.exp(-dwp/2.)*np.exp(1j*wp)
    p = np.concatenate((p, np.conj(p)))

    return p

