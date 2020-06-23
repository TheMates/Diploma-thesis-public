function [Bm,Am,FIR] = filterDesignMatlab(W,target, freqStart, freqStop, crossFreq, crossLength,nPoles1, nPoles2, lambda1, lambda2, NFIR,ITER,Fs,natSpline)

startW = 2*pi*freqStart/Fs;
start = find(W>startW,1);
stopW = 2*pi*freqStop/Fs;
stop = find(W>stopW,1);
if isempty(stop)
    stop = length(W);
end

TF = bankCore.minphasen(target(start:stop),W(start:stop),2^10,[],natSpline); 

C=find(W>2*pi*crossFreq/Fs,1); 

warning('off')
pdualwarp=bankCore.dualwarppolesfr(W(start:stop),1,TF,C-start,crossLength,lambda1,lambda2,nPoles1,nPoles2,ITER,[],natSpline);
[Bm,Am,FIR]=bankCore.parfiltdesfr(W(start:stop),TF,pdualwarp,NFIR);
warning('on')

end