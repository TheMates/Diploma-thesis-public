function MSE = responseMSE(H1,H2)
% H1 - first frequency response in dB
% H2 - second frequency response in dB

MSE = immse(H1,H2);
end

