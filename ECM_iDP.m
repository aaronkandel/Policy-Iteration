% Dynamic Programming via Policy Iteration
% Aaron Kandel

clc
clear
close all
%% Define Transition Tables:

% Define Params:
t_max = 400; % [s] Max time
dt = 7.5; % [s] Timestep
maxit = t_max/dt; % iterations
z_targ = 0.8; % [-] Target SOC
z0 = 0.2; % [-] Initial SOC

% Battery Parameters:
C_Batt = [2.3*3600];% [F]
C1 = [2500];% [F]
R1 = [0.01];% [Ohms]
R2 = [0.01];% [Ohms]
% Parameter Structure:
p.C_Batt = C_Batt;
p.C1 = C1;
p.R1 = R1;
p.R2 = R2;
p.dt = dt;


% Discretize states and inputs:
socD = 0.2:0.00125:(z_targ); % [-] State-of-charge discretization
vrcD = -0.02:0.005:0.2; % [V] Capacitor Voltage discretization
Imin = 0;% [A] Minimum Input Current
Imax = 46;% [A] Maximum Input Current
Iit = 1;% [A] Coarseness of Input Discretization
currentD = Imin:Iit:Imax; % [A] Input Current Discretization

% Generate Vector of all possible states:
xvec = zeros(2,length(socD)*length(vrcD)); 
k = 1;
for i = 1:length(socD)
    for j = 1:length(vrcD)
        xvec(:,k) = [socD(i);vrcD(j)];
        k = k + 1;        
    end % END FOR
end % END FOR


% Battery Model:
% Define Equivalent Circuit Model Linear Dynamics:
A = [1, 0;
     0, (1-dt/(R1*C1))];
B = [dt/C_Batt;dt/C1];
% Define Equivalent Circuit Model Nonlinear Output Equation:
% Load VOC Data (Experimental Battery OCV-SOC Curve):
VOC_data = csvread('Voc.dat',1,0);
soc = VOC_data(:,1);
voc = VOC_data(:,2);
VsimMax = 3.6; % [V] Voltage Constraint


% Assemble table of deterministic state transitions:
% With simple modifications, this code can handle probabilistic state
% transitions from a Markov chain model.  However, for this
% application study the underlying system is deterministic, so the
% transition probabilities are either zero or one.
dsoc = (socD(2) - socD(1));
dvrc = (vrcD(2) - vrcD(1));%./min(x1mesh);
for k = 1:length(currentD)
    clear Vsim VOC
    nxvec = A*xvec + B*currentD(k); % Vector of next states
    % Nearest Neighbor Interpolation:
    nxvSOC = round(nxvec(1,:)./dsoc) .* dsoc;
    nxvVRC = round(nxvec(2,:)./dvrc) .* dvrc;
    nxv = [nxvSOC;nxvVRC];
    % Compute indices of next states:
    for i = 1:length(nxv)
        [~,nearIndSOC] = min((nxvec(1,i) - socD).^2);
        [~,nearIndVRC] = min((nxvec(2,i) - vrcD).^2);
        nextInd(i,k) = length(vrcD)*(nearIndSOC-1) + (nearIndVRC);     
        
        % Compute voltage output for constraint satisfaction:
        VOC(1,i) = voc(soc == round(nxvec(1,i),3));
        Vsim(1,i) = VOC(1,i) + nxvec(2,i) + currentD(k).*R2;
    end % END FOR
    
    % Transition cost (SOC reference tracking + constraint penalty):
    tCost(:,k) = (  (nxv(1,:) - z_targ).^2 + 1000*(Vsim > VsimMax)  )'; 

end % END FOR
%% Policy Iteration:

nInputs = length(currentD);
nStates = length(nxv);
discountFactor = 0.998;
optimalPolicy = PolicyIter(nextInd,tCost,nInputs,nStates,discountFactor);
%% Simulate Final Solution:

% Initialize states:
SOC = z0;
VOC = voc(soc == round(SOC,3));%
vrc = 0;
Vsim = VOC;
action = 0;

% Simulate Greedy Policies:
for i = 1:round(maxit)            
    % Find nearest states:
    [socNear,nearIndSOC] = min((socD - SOC(i)).^2);
    [vrcNear,nearIndVRC] = min((vrcD - vrc(i)).^2);
    % Compute state index:
    curState = length(vrcD)*(nearIndSOC-1) + (nearIndVRC);%  
    
    % Use state index with optimal policy lookup table to get optimal input:
    action(i,1) = currentD(optimalPolicy(curState));

    % Compute State Transitions:
    [SOC(i+1),vrc(i+1),VOC(i),Vsim(i)] = env(...
    SOC(i),vrc(i),action(i),voc,soc,p);

end % END FOR

% PLOT RESULTS:
t2 = dt*(0:(maxit+1)); % assign time vector for plotting

figure(1) 
clf
subplot(1,3,1)
hold on
plot(t2(1:length(action)-1),action(1:end-1))%,'Linewidth',2)
xlim([0 t_max])
grid on
xlabel('Time [s]')
ylabel('Current [A]')
subplot(1,3,2)
hold on
plot(t2(1:length(SOC)),SOC)
plot([0,t2(end)],[z_targ, z_targ],'--k')
grid on
xlabel('Time [s]')
ylabel('SOC [-]')
legend('Optimal','Target') % 'Eng.'
ylim([0 1])
xlim([0 t_max])
subplot(1,3,3)
hold on
plot(t2(1:length(Vsim)),Vsim)
plot([0 t_max],[VsimMax,VsimMax])%t,V_S)
ylim([3 4])
xlim([0 t_max])
grid on
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
legend('Optimal','Constraint') %
%% Functions:

function [SOCn,vrcn,VOCn,Vsimn] = env(SOC,vrc,u,voc,soc,p)
% Compute State Transitions:
SOCn = SOC + (p.dt/p.C_Batt)*u;
VOCn = voc(soc == round(SOCn,3));
vrcn = vrc-(p.dt/(p.R1*p.C1))*vrc + p.dt/p.C1*u;
Vsimn = VOCn + vrcn + u.*p.R2; %       
end



