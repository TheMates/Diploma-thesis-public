function MSE = weightedResponseMSE(H1,H2,Weight)

if Weight == 1
    Weight = ones(length(H1),1);
end

MSE = sum(Weight.*(H2-H1).^2)/length(H1);
end

