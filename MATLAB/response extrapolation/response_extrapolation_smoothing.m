% Demonstrates methods for extrapolation of frequency response and
% smoothing

clc
clearvars

mydir  = pwd;idcs   = strfind(mydir,filesep);root_dir = mydir(1:idcs(end)-1); %one directory up
addpath([root_dir filesep 'common_code'])

load('iPadOrig.mat');

Fs = 88200;

semilogx(fr,db(H),'LineWidth',1.5) 
hold on

w = fr2w(fr,Fs);

[frsm,Hsm] = smoothTF(w,H,Fs,6);    %also oversamples the response to set density (100 bins per octave)

semilogx(frsm,db(Hsm))

[frLow,HLow] =  extrapolateFreqRespLow(frsm,db(Hsm),Fs,10,2,false); %returns in dB scale
[frHigh,HHigh] = extrapolateFreqRespHigh(frsm,db(Hsm),Fs,Fs/2,1,false);

semilogx(frLow,HLow)
semilogx(frHigh,HHigh)
hold off
grid on

xlabel('f [Hz] \rightarrow')
ylabel('Magnitude [dB] \rightarrow')
legend('original','smoothed','low freq extrapolation','high freq extrapolation');
