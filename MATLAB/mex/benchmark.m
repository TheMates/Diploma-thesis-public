mydir  = pwd;idcs   = strfind(mydir,filesep);root_dir = mydir(1:idcs(end)-1); %one directory up
addpath([root_dir filesep 'common_code'])

load([root_dir '\devices data\Engl.mat'])

Fs = 44100;

w = 2*pi*fr/Fs;

%% benchmark test

%average 57 ms :(
Fs = 44100;
N =200;
tic
for i=1:N
    [Bm, Am, FIR] = ParallelFilterDesignMex( H,w,Fs,0.986,0.65,18,22,500,50,1,false);
end
time = toc;
disp(['average time C++ :' num2str(time/N) 's']);

Am1 = Am;
%%
% 1.4x faster
% average 41 ms
tic
for i=1:N
    minph = bankCore.minphasen(H,w,2^14);
    C=find(w>2*pi*500/Fs,1);
    p = bankCore.dualwarppolesfr(w,1,H,C,50,0.986,0.65,18,22,5);
    [Bm, Am, FIR] = bankCore.parfiltdesfr(w,minph,p,1);
end
time = toc;
disp(['average time MATLAB :' num2str(time/N) 's']);
