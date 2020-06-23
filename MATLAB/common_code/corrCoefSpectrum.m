function r = corrCoefSpectrum(H1,H2)
r = sum(H1.*H2)/sqrt(sum(H1.^2)*sum(H2.^2) );
end

