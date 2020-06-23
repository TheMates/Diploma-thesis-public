function wd = dewarpFrequency(w,len,lambda)
  
lambda2 = lambda^2;
dx = len.*cos(w);
wd =  acos((2.* lambda + dx)./(lambda2 + 1));
end