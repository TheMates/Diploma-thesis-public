%% Analysis
% I have 4 different responses - Engl, Numark, Marshall, Fabia
% I have 6 fitness - responseMSE, fracOctMSE, fracThridMSE, criticalMSE, pearson, spectdev
% Make 100 designs
% Save results to file

clc
clearvars

mydir  = pwd;idcs   = strfind(mydir,filesep);root_dir = mydir(1:idcs(end)-1); %one directory up
addpath([root_dir filesep 'common_code'])

NPSO = 100;
NPOLES = 48;

sampleRates  = [44100 44100 48000 48000];
files = {'Engl','Marshall','Numark','Fabia'};


for device = 1:length(sampleRates)
    
    Fs = sampleRates(device);
    path = files{device};
    %% PROBLEM PARAMETERS
    problem.lambda1.Val = 0.979;
    problem.lambda1.MinVal = 0.8;
    problem.lambda1.MaxVal = 0.994;
    
    problem.lambda2.Val = 0.65;
    problem.lambda2.MinVal = 0.5;
    problem.lambda2.MaxVal = 0.9;
    
    problem.crossFreq.Val = [];
    problem.crossFreq.MinVal = 300;     
    problem.crossFreq.MaxVal = 3000;
    
    problem.crossLength.Val = [];
    problem.crossLength.MinVal = 0;
    problem.crossLength.MaxVal = 100;
    
    problem.nPoles = NPOLES;
    problem.Fs = Fs;
    problem.NFIR = 1;
    
    problem.nPolesLow.Val = NPOLES/2;
    problem.nPolesLow.MinVal = NPOLES/2 - NPOLES/4;
    problem.nPolesLow.MaxVal = NPOLES/2 + NPOLES/4;
    
    problem.frequencyWeight = true;
    
    problem.levelWeight.Flag = true;
    problem.levelWeight.Threshold = -15;
    problem.levelWeight.Ratio = 3;
    problem.levelWeight.KneeWidth = 5;
    
    %% PSO PARAMETERS
    kappa = 1;                      
    phi1 = 2.05;
    phi2 = 2.05;
    phi = phi1 + phi2;              
    chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
    
    params.MaxIt = 35;             
    params.nPop = 15;              
    
    params.w = chi;             %inertia coefficient, standardní PSO,
    params.wdamp = 0.96;        %damping ration of inertia coefficient
    params.c1 = chi*phi1;       %personal acceleration coefficient
    params.c2 = chi*phi2;       %social acceleration coefficient
    params.c3 = 0;
    
    params.particleOut = 'mirror';
    %% INPUTS
    %% responseMSE
    problem.CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
        responseMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    
    problems = repmat(problem,NPSO,1);
    for i = 1:NPSO
        problems(i).W = 2*pi*getfield(load([path '.mat']),'fr')/Fs;
        problems(i).Target = bankCore.minphasen(getfield(load([path '.mat']),'H'),problems(i).W,2^15);
    end
    
    
    %% PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    respMSEres = saveResults(out,problems);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fracOct
    for i = 1:NPSO
        problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
            fracOctMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    end
    % PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    fracOctMSEres = saveResults(out,problems);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fracThird
    for i = 1:NPSO
        problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
            fracThirdMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    end
    % PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    fracThirdMSEres = saveResults(out,problems);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% critical
    for i = 1:NPSO
        problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
            criticalBandMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    end
    % PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    criticalMSEres = saveResults(out,problems);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Pearson
    for i = 1:NPSO
        problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
            responsePearsonCostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    end
    % PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    pearsonMSEres = saveResults(out,problems);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SpectDev
    for i = 1:NPSO
        problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
            responseSpectDeviationCostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
    end
    % PSO
    tic
    parfor i = 1:NPSO
        out(i) = parfiltPSO(problems(i), params);
    end
    toc
    
    spectDevMSEres = saveResults(out,problems);
    
    %% Save all
    
    w = problems(1).W;
    H = problems(1).Target;
    
    save([path 'Results.mat'],'w','H','Fs','respMSEres','fracOctMSEres','fracThirdMSEres','criticalMSEres','pearsonMSEres','spectDevMSEres');
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%
function result = saveResults(out,problems)
NFIR = 1;

for i = 1:length(out)
    filterParams = out(i).BestSolution;
    % filter repsonse
    w = problems(i).W;
    lambda1 = filterParams.Position(1,1);
    lambda2 = filterParams.Position(1,2);
    crossFreq = filterParams.Position(1,3);
    crossLength = round(filterParams.Position(1,4));
    nPoles1 = round(filterParams.Position(1,5));
    nPoles2 = problems(i).nPoles - nPoles1;
    
    ITER = 5;
    C=find(w>crossFreq,1);
    TF = problems(i).Target;
    
    warning('off')
    pdualwarp=bankCore.dualwarppolesfr(w,1,TF,C,crossLength,lambda1,lambda2,nPoles1,nPoles2,ITER);
    [Bm,Am,FIR]=bankCore.parfiltdesfr(w,TF,pdualwarp,NFIR);
    warning('on')
    
    H = bankCore.parfiltfresp(Bm,Am,FIR,w);
    
    result.H.(['H' num2str(i)]) = H;
    result.Cost.(['Cost' num2str(i)]) = out(i).BestSolution.Cost;
end


end
