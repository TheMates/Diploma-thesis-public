function threshold = hearingThreshold(f,scaled)
% Computes hearing threshold on specified frequencies in Hz
% Ernst Terhardt - Calculating virtual pitch
% scaled - true - scales hearing threshold so max freq will be 1
if nargin<2
    scaled = true;
end
lowEnd = -0.8;      %-0.8 is standard, but usable in between -0.7..-1.1
highEnd = 4;      %4 is standard, but usable in between 3.75..4
threshold = 3.64*(f./1000).^(lowEnd) - 6.5*exp(-0.6*(f./1000-3.3).^2 ) + (10^-3)*(f./1000).^highEnd;

%20 kHz = weight 0
if scaled
    fmax = 20000;
    cal = 3.64*(fmax./1000).^(-0.8) - 6.5*exp(-0.6*(fmax./1000-3.3).^2 ) + (10^-3)*(fmax./1000).^4;
    threshold = threshold/cal;
    threshold(threshold>1)=1;
end

end

