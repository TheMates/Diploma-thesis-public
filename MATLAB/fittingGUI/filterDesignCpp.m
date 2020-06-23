function [Bm,Am,FIR] = filterDesignCpp(W,target, freqStart, freqStop, crossFreq, crossLength,nPoles1, nPoles2, lambda1, lambda2, NFIR,ITER,Fs,useNAKspline)

startW = 2*pi*freqStart/Fs;
start = find(W>startW,1);
stopW = 2*pi*freqStop/Fs;
stop = find(W>stopW,1);
if isempty(stop)
    stop = length(W);
end


[Bm,Am,FIR] = ParallelFilterDesignMex(abs(target(start:stop)),W(start:stop),Fs,lambda1,lambda2,nPoles1,nPoles2,crossFreq,crossLength,NFIR,logical(useNAKspline));


end
