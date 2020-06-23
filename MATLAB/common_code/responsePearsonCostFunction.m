function r = responsePearsonCostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)

lambda1 = params(1,1);
lambda2 = params(1,2);
crossFreq = params(1,3);
crossLength = params(1,4);
nPoles1 = round(params(1,5));
nPoles2 = nPoles-nPoles1;

ITER = 5;
C=find(w>crossFreq,1); 

warning('off')
pdualwarp=bankCore.dualwarppolesfr(w,1,target,C,crossLength,lambda1,lambda2,nPoles1,nPoles2,ITER);
[Bm,Am,FIR]=bankCore.parfiltdesfr(w,target,pdualwarp,NFIR);
warning('on')

filterResp = bankCore.parfiltfresp(Bm,Am,FIR,w);

if freqWeightFlag
    Weight = rescale(-hearingThreshold(Fs*w/(2*pi)));           %freq weight
    if levelWeightFlag.Flag
        Weight = Weight.*expanderWeight(db(target),levelWeightFlag.Threshold,levelWeightFlag.Ratio,levelWeightFlag.KneeWidth);       %freq + level weight
    end
    r = 1 - wcorrcoef(db(target),db(filterResp),Weight);
else if levelWeightFlag.Flag
    Weight = expanderWeight(db(target),levelWeightFlag.Threshold,levelWeightFlag.Ratio,levelWeightFlag.KneeWidth);                   %level weight
    r = 1 - wcorrcoef(db(target),db(filterResp),Weight);
else
r = 1 - wcorrcoef(db(target),db(filterResp),1);             %no weight
end


end

