function b = criticalBandsZwicker(fr)
% function b = criticalBandsZwicker(fr)
% Generates set of critical bands by Zwicker definition.
% 
% Transforms frequency to bark scale
% 13*atan(0.76.*fr/1000)+3.5*atan((fr/7500).^2)
% 
% 
% returns b - columns fb, fu, fc (bottom, center, up freq)
% 

% Standard
% center = [ 50 150 250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400 4000 4800 5800 7000 8500 10500 13500]';
% bot =    [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500 12000]';
% top =    [100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500 12000 15500]';
% 
% b = [bot center top];

bfr = 13*atan(0.76.*fr/1000)+3.5*atan((fr/7500).^2);
chidx = find(diff(fix(bfr))==1);   
bot = zeros(length(chidx)+1,1);
top = zeros(length(chidx)+1,1);

bot(1) = fr(1);
top(1) = fr(chidx(1)+1);

for i = 2:length(chidx)
    bot(i) = fr(chidx(i-1)+1);
    top(i) = fr(chidx(i)+1);
end
bot(end) = top(end-1);
top(end) = fr(end);

center = mean([bot top],2);
b = [bot center top];
end

