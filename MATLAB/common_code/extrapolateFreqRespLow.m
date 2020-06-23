function [fr, expH] = extrapolateFreqRespLow(fr, H, Fs, freqStart, nOctavesFit, useHPF)
% function [fr, H] = extrapolateFreqRespLow(fr, H, Fs, freqStart, nOctavesFit,useHPF)
% 
% function that extends given frequency response H (in dB) to lower frequencies.
% The first N octaves of characteristics (in dB and log fr) will be fitted with line
% and with this trend the frequency response will be extended
%         
% 
% fr - frequencies vector Hz
% H - frequency response dB
% Fs - sampling frequency
% freqStart - first frequency in extended freq vector
% nOctavesFit - number of octaves, that will be fitted
% useHPF - bool, uses butterworth HPF 4. order at 40 Hz
% 

if nargin<6
    useHPF = false;
end

%% fitting of first N octaves
frIn = fr;
f1 = frIn(1);
f2_i = find(frIn>2^nOctavesFit*f1,1);         %index of stop fitting

fr_fit = frIn(1:f2_i);

myfit = fittype('a + b*log(x)',...          % we want constant slope down in log freq, so the fit is logarithmic
'dependent',{'y'},'independent',{'x'},...
'coefficients',{'a','b'});

P=fit(fr_fit,H(1:f2_i),myfit,'StartPoint',[1 1]);
%we dont want rising slope towards DC, so if P.b <0, we set it to 0
if P.b<0
    P.a = mean(H(1:f2_i));
    P.b = 0;
end
fitted = P.a + P.b*log(fr_fit);                 % slope line at defined frequencies

%% extrapolation
nPerOct = round(length(fr_fit)/2);
octaves = log2(f1/freqStart);
points = round(octaves*nPerOct);

lowf = logspace(log10(freqStart),log10(f1),points)';       %extrapolated freq vector
lowH = P.a + P.b*log(lowf(1:end-1));                %excluding the last, which is the same as in fr
%fitted line + first nOctaves of H
smoothf = [lowf(1:end-1); frIn(1:f2_i)];             
smoothH = [lowH; H(1:f2_i)];

[logscale,smoothmagn]=smoothTF(smoothf,10.^(smoothH/20), Fs, 3);    %smooth the transition from fit to freq resp

%different fr points, interpolate to original smoothf
smoothH = interp1(logscale,db(smoothmagn),smoothf,'spline');    

%% cross fade transition area
crossWeight = linspace(0,1,f2_i);
lowHOrig = 10.^(H(1:f2_i)/20).*crossWeight';        
lowNew = 10.^(smoothH(end-f2_i+1:end)/20).*(flip(crossWeight))';

smoothH(end-f2_i+1:end) = db(lowNew+abs(lowHOrig));

expH = [smoothH; H(f2_i+1:end)];     %now the whole reponse

%% OUT

fr = [lowf(1:end-1);frIn];

%butter 4. order, f = 40 Hz
if useHPF
    cutoff = 40;
    w = fr*2*pi/Fs;
    [Bh,Ah]=butter(4,cutoff/(Fs/2),'high'); 
    HPF=freqz(Bh,Ah,w);
    expH = db(10.^(expH/20).*abs(HPF) );
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

