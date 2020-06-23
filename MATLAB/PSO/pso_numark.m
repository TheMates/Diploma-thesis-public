% 17x10 img size
% Filter design for Numark - PSO
% compute with all fitness, make 4 designs, calculate average cost
% 
% save filter response and cost to file

clc
clear

mydir  = pwd;idcs   = strfind(mydir,filesep);root_dir = mydir(1:idcs(end)-1); %one directory up
addpath([root_dir filesep 'common_code'])

Fs  = 48000;
file = 'Numark.mat';

NPOLES = 48;

repair = false;

%% PROBLEM PARAMETERS
problem.lambda1.Val = 0.979;
problem.lambda1.MinVal = 0.8;
problem.lambda1.MaxVal = 0.994;

problem.lambda2.Val = 0.65;
problem.lambda2.MinVal = 0.5;
problem.lambda2.MaxVal = 0.9;

problem.crossFreq.Val = [];
problem.crossFreq.MinVal = 500;
problem.crossFreq.MaxVal = 4000;

problem.crossLength.Val = [];
problem.crossLength.MinVal = 0;
problem.crossLength.MaxVal = 200;

problem.nPoles = NPOLES;
problem.Fs = Fs;
problem.NFIR = 1;

problem.nPolesLow.Val = NPOLES/2;
problem.nPolesLow.MinVal = NPOLES/2 - NPOLES/4;
problem.nPolesLow.MaxVal = NPOLES/2 + NPOLES/4;

problem.frequencyWeight = true;

problem.levelWeight.Flag = false;
problem.levelWeight.Threshold = -15;
problem.levelWeight.Ratio = 4;
problem.levelWeight.KneeWidth = 5;

%% PSO PARAMETERS
kappa = 1;                      % 0..1
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;              %greater than 4
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

params.MaxIt = 35;             
params.nPop = 15;              
%params.nPop = floor(10+2*sqrt(5));     %SPSO 2006

%coeffs
params.w = chi;             %inertia coefficient
params.wdamp = 0.96;        %damping ration of inertia coefficient, vždycky po iteraci
params.c1 = chi*phi1;       %personal acceleration coefficient
params.c2 = chi*phi2;       %social acceleration coefficient
params.c3 = 0;              %mean of dimension coefficient - from Fast Coverging PSO

params.particleOUT = 'mirror';
%% PROBLEM INPUT DATA
%% responseMSE
problem.CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    responseMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);

%4 designs
problems = repmat(problem,4,1);
for i = 1:4
    problems(i).W = 2*pi*getfield(load(file),'fr')/Fs;
    problems(i).Target = bankCore.minphasen(getfield(load(file),'H'),problems(i).W,2^15);
end


%% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

respMSEres = saveResults(out,problems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fracOct
for i = 1:4
    problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    fracOctMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
end
% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

fracOctMSEres = saveResults(out,problems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fracThird
for i = 1:4
    problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    fracThirdMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
end
% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

fracThirdMSEres = saveResults(out,problems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% critical
for i = 1:4
    problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    criticalBandMSECostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
end
% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

criticalMSEres = saveResults(out,problems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pearson
for i = 1:4
    problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    responsePearsonCostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
end
% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

pearsonMSEres = saveResults(out,problems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SpectDev
for i = 1:4
    problems(i).CostFunction = @(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag)...
    responseSpectDeviationCostFunction(target, w, Fs, nPoles, NFIR, params, freqWeightFlag, levelWeightFlag);
end
% PSO
tic
parfor i = 1:4
    out(i) = parfiltPSO(problems(i), params);
end
toc

spectDevMSEres = saveResults(out,problems);

%% Save all

w = problems(1).W;
H = problems(1).Target;

save('NumarkPSOResultsFreqW.mat','w','H','Fs','respMSEres','fracOctMSEres','fracThirdMSEres','criticalMSEres','pearsonMSEres','spectDevMSEres');


%% %%%%%%%%%%%%%%%%%%%%%%% 
function result = saveResults(out,problems)
    NFIR = 1;
    
    for i = 1:4
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




