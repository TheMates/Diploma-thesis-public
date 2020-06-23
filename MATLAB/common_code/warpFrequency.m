function w_wp = warpFrequency(w,lambda)

lambda2 = lambda^2; 
% len = sqrt(((1+lambda2).*cos(w)-2*lambda).^2 + ((1-lambda2).*sin(w)).^2);   %important for dewarping
w_wp=atan2((1-lambda2).*sin(w),(1+lambda2).*cos(w)-2*lambda);

end

