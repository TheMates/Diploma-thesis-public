function r = wcorrcoef(x,y,w)
% weighted normalized correlation (Pearson) coefficient
% https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Weighted_correlation_coefficient

if length(w) == 1
    w = w*ones(length(x),1);
end

r = wcov(x,y,w)./sqrt(wcov(x,x,w).*wcov(y,y,w));
end

function m = wmean(x,w)
% weighted mean
m = sum(w.*x)./sum(w);
end

function c = wcov(x,y,w)
% Weighted covariance
xwm = wmean(x,w);
ywm = wmean(y,w);

c = sum(w.*(x-xwm).*(y-ywm))./sum(w);
end
