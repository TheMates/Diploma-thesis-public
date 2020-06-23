function [MSE] = criticalBandMSE(H1,H2,w,Fs,weight)
% H1 - first frequency response in dB
% H2 - second frequency response in dB
% w - frequency vector, angular or absolute
% Fs - optional sampling frequency
% weight - vector of weights - weighted MSE

if nargin<5 || isempty(Fs)
    Fs = 44100;
end

if~exist('weight','var')
    weight = 1;
end

if w(end)<=pi+0.001 
    w = Fs*w/(2*pi);
end

cb1 = criticalBandAnalysis(w,H1,Fs,false,false);
cb2 = criticalBandAnalysis(w,H2,Fs,false,false);

MSE = weightedResponseMSE(cb1(:,2),cb2(:,2),weight);

end

