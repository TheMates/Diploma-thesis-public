function [bandSPL,H] = criticalBandAnalysis(w,H,Fs,apply_A_weight,plot_figure,color)
% Computes energy of given frequency response in fractional-octave bands with A-Weighting.
% 
% w - frequency vector (absolute or angular)
% H - magnitude frequency response in dB
% method - 'oct' for octave bands, 'third' for third-octave bands
        
if nargin <6
    color = 'r';
end
if nargin <5
    plot_figure = false;
end

if nargin <4
    apply_A_weight = true;
end


if w(end)<=pi+0.001 
    w = Fs*w/(2*pi);
end


bands = criticalBandsZwicker(w);
nbands = size(bands,1);

if apply_A_weight
    %[OBSOLETE]
    %Norm data to max value of 94 - dBSPL
    H = H./max(H) + 93;
    
    % filter with A weight
    
    % A-weighting filter coefficients
    c1 = 12194.217^2;
    c2 = 20.598997^2;
    c3 = 107.65265^2;
    c4 = 737.86223^2;
    
    % evaluate the A-weighting filter in the frequency domain
    f = w.^2;
    num = c1*(f.^2);
    den = (f+c2) .* sqrt((f+c3).*(f+c4)) .* (f+c1);
    A = db(1.2589*num./den);
    
    H = H+A;
    
    %remove data below threshold
    idx = find(H<0);
    H(idx) = [];
    w(idx) = [];
end

bandSPL = [bands(:,2) zeros(nbands,1)];

for band = 1:nbands 
        
    wlow = max(bands(band,1),w(1));
    whigh =min(bands(band,3),w(end));
    
    if any(intersect(find(wlow<=w), find(w<whigh)))
        idx = intersect(find(wlow<=w), find(w<whigh));
        bandSPL(band,2) = mean(H(idx));
    end
end

if nargout<1
    semilogx(w,H,'LineWidth',1);
    hold on
    bandF = [];
    bandLvl = [];
    for i=1:nbands
        bandF = [bandF; bands(i,1); bands(i,3); bands(i,3)];
        bandLvl = [bandLvl; bandSPL(i,2); bandSPL(i,2);bandSPL(min(nbands,i+1),2)];
    end
    bandF(1,1) = w(1);
    semilogx(bandF,bandLvl,'color',color);
    hold off
    grid on
end

if plot_figure
        hold on
        bandF = [];
        bandLvl = [];
    for i=1:nbands
        bandF = [bandF; bands(i,1); bands(i,3); bands(i,3)];
        bandLvl = [bandLvl; bandSPL(i,2); bandSPL(i,2);bandSPL(min(nbands,i+1),2)];
    end
    bandF(1,1) = w(1);
    semilogx(bandF,bandLvl,'color',color);
end





