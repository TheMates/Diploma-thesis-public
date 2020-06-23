function [poles] = drawPlot(Bm,Am,FIR,W,Fs,target)

filterH=bankCore.parfiltfresp(Bm,Am,FIR,W); % freq. response of the filter

%% plotting
fr = Fs*W/(2*pi);

DBSTEP=20; % curves are offset by DBSTEP dB
semilogx(fr,db(target),'g-.'); 
hold on;

%matlab data
semilogx(fr,db(filterH),'r');

poles = [];
for i=1:size(Am,2)
    p = roots(Am(:,i));
    poles = [poles ;p(:)];
end

fp=Fs*angle(poles')/(2*pi);
plot(fp,-DBSTEP/4,'kx');
hold off;
grid on
warning('off')
axis([10 Fs/2 -50 30]);
warning('on')


SIZE=10;
set(gca,'FontName','Times','Fontsize',SIZE);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');   
title(strcat('Frequency response modeling with a ', num2str(length(poles)), 'th order parallel filter'));

legend('target','filter','poles')


end