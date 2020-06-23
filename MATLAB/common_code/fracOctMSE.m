function [MSE] = fracOctMSE(H1,H2,method,w,Fs,weight)
% H1 - first frequency response in dB
% H2 - second frequency response in dB
% method - 'third' or 'oct'
% w - frequency vector, angular or absolute
% Fs - optional sampling frequency
% weight - vector of weights - weighted MSE. Must be same length as number
% of octave analysis

if nargin<5 || isempty(Fs)
    Fs = 44100;
end

if~exist('weight','var')
    weight = 1;
end

if w(end)<=pi+0.001 
    w = Fs*w/(2*pi);
end

frac1 = fracOctAnalysis(w,H1,Fs,method,false,false);
frac2 = fracOctAnalysis(w,H2,Fs,method,false,false);

MSE = weightedResponseMSE(frac1(:,2),frac2(:,2),weight);

end

