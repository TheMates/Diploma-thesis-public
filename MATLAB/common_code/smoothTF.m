function [logscale,smoothmagn]=smoothTF(W,H, Fs, fract, avg)

%TFPLOTS - Smoothed transfer fucntion plotting
%   [FREQ,MAGN]=TFPLOTS(W,H, Fs, FRACT, AVG, WINDOW)
%   Fractional-octave smoothing of given frequency response defined on angular frequencies W.
%   First frequency response is transformed to linear scale by spline, then the data is colleced
%   into logaritmically spaced bins and the average is computed for
%   each bin (100/octave). Then this is smoothed (convolved) by a hanning window, where
%   FRACT defines the fractional-octave smoothing (default is 3, meaning third-octave).
%   The length of the smoothing hanning window is the double compared to the distance
%   defined by FRACT.
%   The sampling frequency is set by FS (default is 44.1 kHz) and the plotting color is set by the COLOR variable
%   (default is 'b').
%
%   If the AVG variable is set to 'power' then the power is averaged
%   in the logaritmic bins and during smoothing (this is the default - 
%   on the contrary to the TFPLOT function, where 'comp' is the default),
%   if it is 'abs' then the absolute value, and if to 'comp',
%   it averages the complex transfer function.
%
%
%   If the output arguments FREQ and MAGN is not asked for, then it plots
%   semilogx(FREQ,20*log10(MAGN)); 


octbin=100; % this is how many frequency bins are computed per octave


if nargin<5,
    avg='power';
end;

if nargin<4,
    fract=6;
end;

if nargin<3,
    Fs=44100;
end;


FFTSIZE=2^18;
H=H(:);


if(W(end)>pi)
    W = 2*pi*W/Fs;
end


LinW=pi*[0:(FFTSIZE-1)]/(FFTSIZE-1);

H=abs(spline(W,H,LinW)); % we take the absolute value since the spline might produce negative values even from positive input


logfact=2^(1/octbin);
LOGN=floor(log(Fs/2)/log(logfact));
logscale=logfact.^[0:LOGN]; %logarithmic scale from 1 Hz to Fs/2

% remove data, that are not in the freq vector
logscale(logscale< Fs*W(1)/(2*pi))= [];
logscale(logscale> Fs*W(end)/(2*pi)) = [];
LOGN = length(logscale)-1;

magn=(abs(H));
compamp=H;


%creating 100th octave resolution log. spaced data from the lin. spaced FFT data
clear logmagn;
fstep=Fs/FFTSIZE/2;
for k=0:LOGN,
   start=round(logscale(k+1)/sqrt(logfact)/fstep);
   start=max(start,1);
   start=min(start,FFTSIZE);
   stop=round(logscale(k+1)*sqrt(logfact)/fstep);
   stop=max(stop,1);
   stop=min(stop,FFTSIZE);

   if strcmpi(avg,'comp') | strcmpi(avg,'complex'),   %averaging the complex transfer function
       logmagn(k+1)=mean(compamp(start:stop)); 
   end;
   if strcmpi(avg,'abs'), %averaging absolute value
      logmagn(k+1)=mean(abs(compamp(start:stop))); 
   end;
   if strcmpi(avg,'power') | strcmpi(avg,'pow'), %averaging power
       logmagn(k+1)=sqrt(mean(abs(compamp(start:stop)).^2)); 
   end;

    
    
end;

%creating hanning window
HL=2*round(octbin/fract); %fractional octave smoothing
hh=hanning(HL);

L=length(logmagn);
logmagn(L+1:L+HL)=0;

%Smoothing the log. spaced data by convonvling with the hanning window

   if strcmpi(avg,'comp') | strcmpi(avg,'complex'),   %averaging the complex transfer function
       tmp=fftfilt(hh,logmagn); 
       smoothmagn=tmp(HL/2+1:HL/2+L)/sum(hh);
   end;
   if strcmpi(avg,'abs'), %averaging absolute value
       tmp=fftfilt(hh,logmagn);
       smoothmagn=tmp(HL/2+1:HL/2+L)/sum(hh);
   end;
   if strcmpi(avg,'power') | strcmpi(avg,'pow'), %averaging power
       tmp=fftfilt(hh,logmagn.^2); 
       smoothmagn=sqrt(tmp(HL/2+1:HL/2+L)/sum(hh));
   end;

   logscale = logscale(:);
   smoothmagn = smoothmagn(:);
if nargout<1,
    semilogx(logscale,20*log10(abs(smoothmagn)));
end;

