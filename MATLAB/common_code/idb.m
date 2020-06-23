function y = idb(x)
% inverse decibel transform 20*log10(x) <-> 10.^(x/20)
y = 10.^(x/20);
end

