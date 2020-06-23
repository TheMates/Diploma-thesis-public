function [Bm, Am, FIR] = ParallelFilterDesignMex(H,W, Fs,lambda1,lambda2,nPoles1,nPoles2,crossFreq,crossLength,NFIR,useNAKspline)
%function [Bm, Am, FIR] = ParallelFilterDesignMex(H,W, Fs,lambda1,lambda2,nPoles1,nPoles2,crossFreq,crossLength,NFIR,useNAKspline)
% 
% Does the dirty work in C++
% H - magnitude of frequency response
% W - angular frequencies
% Fs - sample rate
% lambda1 - first warping parameter
% lambda2 - second warping parameter
% nPoles1 - number of poles in first band
% nPoles2 - number of poles in second band
% crossFreq - crossing frequency of two bands
% crossLength - number of frequency bins where bands cross
% NFIR - number of FIR coefficieints
% useNAKspline - bool, use not-a-knot or natural spline
% 
% Returns:
% Bm - numerator coefficients
% Am - denominator coefficients
% FIR - FIR part coefficients
% 
% Notes:
% number of iterations of lsidfr = 5, NFIR = 5
%
%
% First minphasen() will be called upon TF, then dualwarppolesfr() will be
% executed and then parfiltdesfr().




