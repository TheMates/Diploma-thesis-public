function MSE = fracOctMSECostFunction(target, w, Fs, nPoles, NFIR, params,freqWeightFlag, levelWeightFlag)

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
    ht = rescale(-hearingThreshold(Fs*w/(2*pi)));           %freq weight
    Weight = fracOctAnalysis(w,ht,Fs,'oct',false);
    Weight = Weight(:,2);
    if levelWeightFlag.Flag
        expw = expanderWeight(db(target),levelWeightFlag.Threshold,levelWeightFlag.Ratio,levelWeightFlag.KneeWidth);                   %level weight
        expw = fracOctAnalysis(w,expw,Fs,'oct',false);
        expw = expw(:,2);
        Weight = Weight.*expw;
    end
    MSE = fracOctMSE(db(target),db(filterResp),'oct',w,Fs,Weight); 
else if levelWeightFlag.Flag
    expw = expanderWeight(db(target),levelWeightFlag.Threshold,levelWeightFlag.Ratio,levelWeightFlag.KneeWidth);                   %level weight
    Weight = fracOctAnalysis(w,expw,Fs,'oct',false);
    Weight = Weight(:,2);
    MSE = fracOctMSE(db(target),db(filterResp),'oct',w,Fs,Weight); 
else
  MSE = fracOctMSE(db(target),db(filterResp),'oct',w,Fs);  
    end

end

