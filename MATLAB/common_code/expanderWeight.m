function [weight,y] = expanderWeight(x,threshold,ratio,kneeWidth)

calib = max(x)+1; 
x = x-calib;        %avoid zero crossing

med = median(x);
threshold = med + threshold;

i_1 = find(x<(threshold-kneeWidth/2));
i_2 = find(x>=(threshold-kneeWidth/2) & x<=(threshold+kneeWidth/2));
i_3 = find(x>(threshold+kneeWidth/2));

x_1 = x(i_1);
x_2 = x(i_2);
x_3 = x(i_3);

y_1 = threshold+(x_1-threshold)*ratio;
y_2 = x_2 + ((1-ratio)*(x_2-threshold-kneeWidth/2).^2  )/(2*kneeWidth);
y_3 = x_3;

y(i_1) = y_1;
y(i_2) = y_2;
y(i_3) = y_3;

weight = x(:)./y(:);
y = y+ calib;
end

