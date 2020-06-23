function [fr, expH] = extrapolateFreqRespHigh(fr, H, Fs, freqStop, nOctavesFit, useLPF)
% function [fr, H] = extrapolateFreqRespHigh(fr, H, Fs, freqStop, nOctavesFit,useHPF)
% 
% function that extends given frequency response H (in dB) to higher frequencies.
% The last N octaves of characteristics (in dB and log fr) will be fitted with line
% and with this trend the frequency response will be extended
%         
% 
% fr - frequencies vector Hz
% H - frequency response dB
% Fs - sampling frequency
% freqStop - last frequency in extended freq vector
% nOctavesFit - number of octaves, that will be fitted
% useHPF - bool, uses butterworth LPF 2. order at 24000 Hz
% 

if nargin<6
    useLPF = false;
end

%% fitting of last N octaves
frIn = fr;

f1_i = find(frIn>2^(-nOctavesFit)*fr(end),1);   %index of start fitting
f2 = fr(end);         


fr_fit = frIn(f1_i:end);

myfit = fittype('a + b*log(x)',...          % we want constant slope down in log freq, so the fit is logarithmic
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});

P=fit(fr_fit,H(f1_i:end),myfit,'StartPoint',[1 1]);
%we dont want rising slope at high freq, so if P.b >0, we set it to 0
if P.b>0
    P.a = mean(H(f1_i:end));
    P.b = 0;
end
fitted = P.a + P.b*log(fr_fit);                 % slope line at defined frequencies

%% extrapolation
nPerOct = round(length(fr_fit)/nOctavesFit);
octaves = log2(freqStop/f2);
points = round(octaves*nPerOct);

highf = logspace(log10(f2),log10(freqStop),points)';       %extrapolated freq vector
highH = P.a + P.b*log(highf(2:end));                  %excluding the first, which is the same as in frIn
%fitted line + last nOctaves of H
smoothf = [ frIn(f1_i:end) ;highf(2:end)];             
smoothH = [H(f1_i:end); highH];

[logscale,smoothmagn]=smoothTF(smoothf,10.^(smoothH/20), Fs, 6);    %smooth the transition from fit to freq resp - third octave

%different fr points, interpolate to original smoothf
smoothH = interp1(logscale,db(smoothmagn),smoothf,'spline');    

%% cross fade transition area
crossWeight = linspace(1,0,length(fr_fit));
highHOrig = 10.^(H(f1_i:end)/20).*crossWeight';        
highNew = 10.^(smoothH(1:length(fr_fit))/20).*(flip(crossWeight))';

smoothH(1:length(fr_fit)) = db(highNew+highHOrig);

expH = [H(1:end-length(fr_fit)); smoothH];     %now the whole reponse

%% OUT

fr = [frIn; highf(2:end)];

%butter 2. order, f = 24000 Hz
if useLPF
    cutoff = 24000;
    w = fr*2*pi/Fs;
    [Bh,Ah]=butter(2,cutoff/(Fs/2),'low'); 
    LPF=freqz(Bh,Ah,w);
    expH = db(10.^(expH/20).*abs(LPF) );
end

%% plot
if nargout <1
    semilogx(fr,expH,'LineStyle','--')
    hold on
    semilogx(fr_fit,fitted,'g--','LineWidth',1)
    semilogx(frIn,H)
  
    hold off
    grid on
   legend('extrapolated','fit','original');
end

end

