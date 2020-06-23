function sigma = spectrumDeviation(H1,H2,weight)
% S. A. Ryder, "Methods for comparing frequency response analysis measurements," 
% Conference Record of the the 2002 IEEE International Symposium on Electrical Insulation 
% (Cat. No.02CH37316), Boston, MA, USA, 2002, pp. 187-190.

if nargin <3
    m = mean([H1 H2],2);
    sigma = sum( sqrt( ((H1-m)./m).^2 + ((H2-m)./m).^2  ));  
    return
end

if length(weight) == 1
   weight = weight*ones(length(H1),1);
end
m = mean([H1 H2],2);
sigma = sum( sqrt( weight.*((H1-m)./m).^2 + weight.*((H2-m)./m).^2  ));

end

