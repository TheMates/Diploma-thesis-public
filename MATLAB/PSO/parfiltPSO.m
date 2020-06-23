function out = parfiltPSO(problem,params)
% function out = parfiltPSO(problem,params)
%
% Particle swarm optimization for input parameters of parallel filter design with dual warping
% 
% problem.Target    - target frequency repsonse - THE MINPHASE VERSION
% problem.W         - angular frequencies of specification
% problem.Fs        - sampling frequency
% 
% problem.lambda1       - lambda for low frequencies, struct of Val, MinVal, MaxVal
% problem.lambda2       - lambda for high frequencies, struct of Val, MinVal, MaxVal
% problem.crossFreq     - cross frequency of dual warping, struct of Val, MinVal, MaxVal
% problem.crossLength   - number of samples where two bands can overlap, struct of Val, MinVal, MaxVal
% problem.nPoles        - total number of poles 
% problem.nPolesLow     - number of poles at low frequency band, number of poles for high frequency band is nPoles-nPolesLow,struct of Val, MinVal, MaxVal
% 
% params        - parameters of PSO
% params.MaxIt  - maximal number of iterations
% params.nPop   - number of particles
% params.w      - intertia coefficient
% params.wdamp  - damping of intertia coefficient
% params.c1     - personal acceleration coefficient
% params.c2     - social acceleration coefficient
% params.c3     - from FCPSO mean of all particles in the dimension
% params.particleOut - type of 

CostFunction = problem.CostFunction;      

%% Parameters of PSO

VarSize = [1 5];        %lambda1, lambda2, crossFreq, crossLength, nPolesLow
MaxIt = params.MaxIt;           
nPop = params.nPop;              

%koeficienty
w = params.w;           %inertia coefficient, standardní PSO, nìkde asi definovaný
wdamp = params.wdamp;   %damping ration of inertia coefficient, vždycky po iteraci
c1 = params.c1;         %personal acceleration coefficient
c2 = params.c2;         %social acceleration coefficient
c3 = params.c3;         %koeficient z prùmìru dimenzí - fast converging pso

% crossFreq to angular
problem.crossFreq.MinVal = 2*pi*problem.crossFreq.MinVal/problem.Fs;
problem.crossFreq.MaxVal = 2*pi*problem.crossFreq.MaxVal/problem.Fs;

% VELOCITY LIMITS
velocityLimit = 0.2;

VelocityLimits.lambda1 = velocityLimit*(problem.lambda1.MaxVal-problem.lambda1.MinVal);
VelocityLimits.lambda2 = velocityLimit*(problem.lambda2.MaxVal-problem.lambda2.MinVal);
VelocityLimits.crossFreq = velocityLimit*(problem.crossFreq.MaxVal-problem.crossFreq.MinVal);
VelocityLimits.crossLength = velocityLimit*(problem.crossLength.MaxVal-problem.crossLength.MinVal);
VelocityLimits.nPolesLow = velocityLimit*(problem.nPolesLow.MaxVal-problem.nPolesLow.MinVal);

if isfield(params,'particleOut')
    particleOut = params.particleOut;
else
    particleOut = 'mirror';
end


%% INITIALIZATION
empty_particle.Position = [];       %pozice
empty_particle.Velocity = [];       %rychlost
empty_particle.Cost = [];           %závisí na aktuálních koeficientech, evaluace CostFunction v této pozici
empty_particle.Best.Position = [];  %Best je taky struct
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop,1);      %nPop øádkù a 1 sloupec

GlobalBest.Cost = inf;                          %nejhorší hodnota

%Initialize population members
for i=1:nPop
    %generate random solution
    
    positions(1,1) = unifrnd(problem.lambda1.Val-0.015, problem.lambda1.Val+0.015);     % lambda1 - initial (0.977), disperse between -0.016 +0.016 - variance of results so far
    positions(1,2) = unifrnd(problem.lambda2.MinVal, problem.lambda2.MaxVal);           % lambda2 - initial (0.65), disperse between 0.5 .. 0.9
    positions(1,3) = unifrnd(problem.crossFreq.MinVal,problem.crossFreq.MaxVal);        % crossFreq - 500-3000, now float, CostFunction will round
    positions(1,4) = unifrnd(problem.crossLength.MinVal,problem.crossLength.MaxVal);    % crossLength - 0-200, now float, CostFunction will round
    positions(1,5) = unifrnd(problem.nPolesLow.Val-5, problem.nPolesLow.Val+5);         % nPolesLow - initial (24), disperse between -5 +5
    
    particle(i).Position = positions;
    
    % Initialize velocity
    velocities(1,1) = unifrnd(-VelocityLimits.lambda1,VelocityLimits.lambda1);
    velocities(1,2) = unifrnd(-VelocityLimits.lambda2,VelocityLimits.lambda2);
    velocities(1,3) = unifrnd(-VelocityLimits.crossFreq,VelocityLimits.crossFreq);
    velocities(1,4) = unifrnd(-VelocityLimits.crossLength,VelocityLimits.crossLength);
    velocities(1,5) = unifrnd(-VelocityLimits.nPolesLow,VelocityLimits.nPolesLow);
    
%     particle(i).Velocity = zeros(VarSize);  %na zaèátu bude 0
    particle(i).Velocity = velocities;  %random velocity       
    
    
    % Evaluation
    
    particle(i).Cost = CostFunction(problem.Target,...
                                    problem.W,...
                                    problem.Fs,...
                                    problem.nPoles,...
                                    problem.NFIR,...
                                    particle(i).Position,...    %nPolesLow
                                    problem.frequencyWeight,...
                                    problem.levelWeight);         
    
    
    % Update personal best, zatím máme jen jeden výsledek, tak to tam dám
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end

BestCosts = zeros(MaxIt,1);     %tady bude best cost v každé iteraci




%% Main Loop of PSO
for it = 1:MaxIt
    
    for i = 1:nPop          %bohužel nepùjde použít parfor, protože se v loopu updatuje global best a ta se používá v další particle
        % Update velocity
        particle(i).Velocity = w*particle(i).Velocity...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position)...
              + c3*rand(VarSize).*(mean(reshape([particle.Position],[VarSize(2),nPop]),2)' - particle(i).Position);
        
        % Apply velocity limits
        particle(i).Velocity(1,1) = max(particle(i).Velocity(1,1), -VelocityLimits.lambda1);     
        particle(i).Velocity(1,1) = min(particle(i).Velocity(1,1),  VelocityLimits.lambda1);
        particle(i).Velocity(1,2) = max(particle(i).Velocity(1,2), -VelocityLimits.lambda2);     
        particle(i).Velocity(1,2) = min(particle(i).Velocity(1,2),  VelocityLimits.lambda2);
        particle(i).Velocity(1,3) = max(particle(i).Velocity(1,3), -VelocityLimits.crossFreq);     
        particle(i).Velocity(1,3) = min(particle(i).Velocity(1,3),  VelocityLimits.crossFreq);
        particle(i).Velocity(1,4) = max(particle(i).Velocity(1,4), -VelocityLimits.crossLength);     
        particle(i).Velocity(1,4) = min(particle(i).Velocity(1,4),  VelocityLimits.crossLength);
        particle(i).Velocity(1,5) = max(particle(i).Velocity(1,5), -VelocityLimits.nPolesLow);     
        particle(i).Velocity(1,5) = min(particle(i).Velocity(1,5),  VelocityLimits.nPolesLow);
   
        
        % Update position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply bounds limits
        switch particleOut
            case 'nearest'  %Nearest neighbor
                particle(i).Position(1,1) = max(particle(i).Position(1,1), problem.lambda1.MinVal);
                particle(i).Position(1,1) = min(particle(i).Position(1,1), problem.lambda1.MaxVal);
                particle(i).Position(1,2) = max(particle(i).Position(1,2), problem.lambda2.MinVal);
                particle(i).Position(1,2) = min(particle(i).Position(1,2), problem.lambda2.MaxVal);
                particle(i).Position(1,3) = max(particle(i).Position(1,3), problem.crossFreq.MinVal);
                particle(i).Position(1,3) = min(particle(i).Position(1,3), problem.crossFreq.MaxVal);
                particle(i).Position(1,4) = max(particle(i).Position(1,4), problem.crossLength.MinVal);
                particle(i).Position(1,4) = min(particle(i).Position(1,4), problem.crossLength.MaxVal);
                particle(i).Position(1,5) = max(particle(i).Position(1,5), problem.nPolesLow.MinVal);
                particle(i).Position(1,5) = min(particle(i).Position(1,5), problem.nPolesLow.MaxVal);
                
            case 'periodic' % Periodic extension
                if particle(i).Position(1,1) >problem.lambda1.MaxVal || particle(i).Position(1,1) <problem.lambda1.MinVal
                    kolik = floor((problem.lambda1.MinVal-particle(i).Position(1,1))/(problem.lambda1.MaxVal-problem.lambda1.MinVal));
                    particle(i).Position(1,1) = particle(i).Position(1,1) + (kolik +1)*(problem.lambda1.MaxVal-problem.lambda1.MinVal);
                    particle(i).Velocity(1,1) = 0;
                end
                if particle(i).Position(1,2) >problem.lambda2.MaxVal || particle(i).Position(1,2) <problem.lambda2.MinVal
                    kolik = floor((problem.lambda2.MinVal-particle(i).Position(1,2))/(problem.lambda2.MaxVal-problem.lambda2.MinVal));
                    particle(i).Position(1,2) = particle(i).Position(1,2) + (kolik +1)*(problem.lambda2.MaxVal-problem.lambda2.MinVal);
                    particle(i).Velocity(1,2) = 0;
                end
                if particle(i).Position(1,3) >problem.crossFreq.MaxVal || particle(i).Position(1,3) <problem.crossFreq.MinVal
                    kolik = floor((problem.crossFreq.MinVal-particle(i).Position(1,3))/(problem.crossFreq.MaxVal-problem.crossFreq.MinVal));
                    particle(i).Position(1,3) = particle(i).Position(1,3) + (kolik +1)*(problem.crossFreq.MaxVal-problem.crossFreq.MinVal);
                    particle(i).Velocity(1,3) = 0;
                end
                if particle(i).Position(1,4) >problem.crossLength.MaxVal || particle(i).Position(1,4) <problem.crossLength.MinVal
                    kolik = floor((problem.crossLength.MinVal-particle(i).Position(1,4))/(problem.crossLength.MaxVal-problem.crossLength.MinVal));
                    particle(i).Position(1,4) = particle(i).Position(1,4) + (kolik +1)*(problem.crossLength.MaxVal-problem.crossLength.MinVal);
                    particle(i).Velocity(1,4) = 0;
                end
                if particle(i).Position(1,5) >problem.nPolesLow.MaxVal || particle(i).Position(1,5) <problem.nPolesLow.MinVal
                    kolik = floor((problem.nPolesLow.MinVal-particle(i).Position(1,5))/(problem.nPolesLow.MaxVal-problem.nPolesLow.MinVal));
                    particle(i).Position(1,5) = particle(i).Position(1,5) + (kolik +1)*(problem.nPolesLow.MaxVal-problem.nPolesLow.MinVal);
                    particle(i).Velocity(1,5) = 0;
                end
                
            case 'mirror'
                if particle(i).Position(1,1) >problem.lambda1.MaxVal || particle(i).Position(1,1) <problem.lambda1.MinVal
                    delka = (problem.lambda1.MaxVal-problem.lambda1.MinVal);
                    kolik = floor((problem.lambda1.MinVal-particle(i).Position(1,1))/delka);
                    particle(i).Position(1,1) = (problem.lambda1.MinVal -kolik*delka) - particle(i).Position(1,1);
                    if mod(kolik,2)
                        particle(i).Position(1,1) = -particle(i).Position(1,1)+ problem.lambda1.MaxVal;
                    else
                        particle(i).Position(1,1) = particle(i).Position(1,1)+ problem.lambda1.MinVal;
                    end
                    particle(i).Velocity(1,1) = 0;
                end
                if particle(i).Position(1,2) >problem.lambda2.MaxVal || particle(i).Position(1,2) <problem.lambda2.MinVal
                    delka = (problem.lambda2.MaxVal-problem.lambda2.MinVal);
                    kolik = floor((problem.lambda2.MinVal-particle(i).Position(1,2))/delka);
                    particle(i).Position(1,2) = (problem.lambda2.MinVal -kolik*delka) - particle(i).Position(1,2);
                    if mod(kolik,2)
                        particle(i).Position(1,2) = -particle(i).Position(1,2)+ problem.lambda2.MaxVal;
                    else
                        particle(i).Position(1,2) = particle(i).Position(1,2)+ problem.lambda2.MinVal;
                    end
                    particle(i).Velocity(1,2) = 0;
                end
                if particle(i).Position(1,3) >problem.crossFreq.MaxVal || particle(i).Position(1,3) <problem.crossFreq.MinVal
                    delka = (problem.crossFreq.MaxVal-problem.crossFreq.MinVal);
                    kolik = floor((problem.crossFreq.MinVal-particle(i).Position(1,3))/delka);
                    particle(i).Position(1,3) = (problem.crossFreq.MinVal -kolik*delka) - particle(i).Position(1,3);
                    if mod(kolik,2)
                        particle(i).Position(1,3) = -particle(i).Position(1,3)+ problem.crossFreq.MaxVal;
                    else
                        particle(i).Position(1,3) = particle(i).Position(1,3)+ problem.crossFreq.MinVal;
                    end
                    particle(i).Velocity(1,3) = 0;
                end
                if particle(i).Position(1,4) >problem.crossLength.MaxVal || particle(i).Position(1,4) <problem.crossLength.MinVal
                    delka = (problem.crossLength.MaxVal-problem.crossLength.MinVal);
                    kolik = floor((problem.crossLength.MinVal-particle(i).Position(1,4))/delka);
                    particle(i).Position(1,4) = (problem.crossLength.MinVal -kolik*delka) - particle(i).Position(1,4);
                    if mod(kolik,2)
                        particle(i).Position(1,4) = -particle(i).Position(1,4)+ problem.crossLength.MaxVal;
                    else
                        particle(i).Position(1,4) = particle(i).Position(1,4)+ problem.crossLength.MinVal;
                    end
                    particle(i).Velocity(1,4) = 0;
                end
                if particle(i).Position(1,5) >problem.nPolesLow.MaxVal || particle(i).Position(1,5) <problem.nPolesLow.MinVal
                    delka = (problem.nPolesLow.MaxVal-problem.nPolesLow.MinVal);
                    kolik = floor((problem.nPolesLow.MinVal-particle(i).Position(1,5))/delka);
                    particle(i).Position(1,5) = (problem.nPolesLow.MinVal -kolik*delka) - particle(i).Position(1,5);
                    if mod(kolik,2)
                        particle(i).Position(1,5) = -particle(i).Position(1,5)+ problem.nPolesLow.MaxVal;
                    else
                        particle(i).Position(1,5) = particle(i).Position(1,5)+ problem.nPolesLow.MinVal;
                    end
                    particle(i).Velocity(1,5) = 0;
                end
            otherwise
        end
        
         % Evaluation
    particle(i).Cost = CostFunction(problem.Target,...
                                    problem.W,...
                                    problem.Fs,...
                                    problem.nPoles,...
                                    problem.NFIR,...
                                    particle(i).Position,...    %nPolesLow
                                    problem.frequencyWeight,...
                                    problem.levelWeight);         

        if particle(i).Cost < particle(i).Best.Cost                      
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;

            %jestliže má lepší Cost než její pøedchozí hodnota, tak je
            %možnost, že je aji novým nejlepším
            if particle(i).Best.Cost < GlobalBest.Cost
%                 if(particle(i).Best.Cost <0)
%                     a = 0;
%                 end                
                GlobalBest = particle(i).Best;
            end
        end
        
    end
    
    % Store the best cost of iteration
    BestCosts(it) = GlobalBest.Cost;
    
    %Display Iteration information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    
    %Damping Inertia Coefficient
    w = w*wdamp;    %utlum, to je dùležitý, jinak bude pøevažovat složka pùvodní rychlosti nad složkami k minimùm
end

disp(['Best Cost = ' num2str(BestCosts(end))]);

out.pop = particle;
out.BestSolution = GlobalBest;
out.BestCosts = BestCosts;


end

