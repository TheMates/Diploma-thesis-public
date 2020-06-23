function bands = fracOctaveStandard(fr,frac)
% returns bands in certain frac oct resolution suitable for given frequency
% vector, centered around 1 kHz
%
% fr - in Hz
% frac - 3 = third, 6 = sixth, 1 = oct, etc..

base = 1000;    %center 1 kHz

startExp = round(frac*log2(fr(1)/base));       % first center band 
stopExp = round(frac*log2(fr(end)/base));      % last center band

exps = (startExp:stopExp)';

center = 1000*2.^(exps/frac);
bot = center*2.^(-1/(2*frac));
top = center*2.^(1/(2*frac));

bands = [bot center top];
end