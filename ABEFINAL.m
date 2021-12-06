clear, close all, format compact, clc
% ABE4649 FINAL GROUP PROJECT 
% 
% The purpose of this project is to model the 3 dimensional "eigenforest"
% Model created by our group. 
% The equations detail change in tree population (dx)
%       change in parasitic beetle population (db, nondim dGamma) 
%       and change in human lumbar mill policy (dh) 
%********** PARAMETERS
%   time: integer (> 0 ) 
%       used to develop iterations for vectors and numerical integration
%   dt : float (> 0) 
%       used to develop time vector for numerical integration
%   gamma: equation ( hHf/rTree)
%       serves as nondimensional unit gamma, reps 
%   hotel: equation (gamma*rTree / hf)
%       serves as nondimensional unit H
%   alfa: equation (Bd/rTree) 
%       serves as nondimensional unit alfa
%   beetleEq: equation 
%       serves as nondimensional unit B
%   beetleDamage: float ((0,1])
%       indicates damage to tree by beetle (1/(beetle*time)) 
%   P: float ((0,+infty)) 
%       indicates profitability per tree harvested ($/tree) 
%   beetleMortality: float ([0,1]) 
%       indicates natural loss of beetles (1/time) 
%   tree: ([0,K])
%       indicates numbers of trees (tree) 
%********* END PARAMETERS

%SETUP PARAMETERS

time = 1500;                % establish time
dt = .001;                  % establish time interval
tVector = [1:dt:time]';      % create vector 1-time @ dt increment

xVector = zeros(size(tVector));
gammaVector = zeros(size(tVector));

beetleMortality = .2;       % natural beetle death rate 1/time
rFumigation = .1;           % fumigation rate 1/time
rBeetle = .3;               % intrinsic beetle 1/time
rHuman = .1;               % intrinsic policy 1/ttime
rTree = .08;                 % intrinisc tree growth 1/time
beetleDamage = .12;         % beetle damage to tree 1/(beetle * tree)
tree = 5000;                   % Trees tree
c = .0024;                     % cost of fumigation $/beetle
q = 8;                      % num. fumigations dimless
f = 1;                      %
P = 4;                      % profitability $/tree
h = .23;                     % harvest effor 1/$
K = 3000;                      % tree carry capacity tree
z = .11;                    % habitability for beetle 1/tree
H = 12;                     % policy effort $ 

% BEGIN PRECOMPUTATIONS 
beetleEq = (rFumigation*H)/((z*rBeetle*tree) - (beetleMortality))
theta = (rHuman*P*h*K)/(rTree)                  %
phi = (c*q*rHuman*h*f)/(rTree*beetleDamage)     %
alfa = (beetleEq*beetleDamage)/(rTree)          %
gamma = (h*H*f)/(rTree)                         %

xVector(1) = 900;                               % init x condition
gammaVector(1) = 55.2;                         % init gamma cond. 

% BEGIN EULER 
for t = 1:(length(tVector)-1)               
    % for each incremental unit time
    %    iterate through the provided nondim eqn to populate xVector
    %    and gammaVector 
    a = alfa;                               
    th = theta;         
    p = phi;
    x = xVector(t);
    g = gammaVector(t);
    xVector(t+1) = xVector(t) + dt*(x*(1-x) - a*x - gamma*x);      %dTree
    gammaVector(t+1) = gammaVector(t) + dt*((th*gamma*x) - (p*x)); %dPolicy
end 

%PLOTTING SYSTEM 

figure(101)                           %initialize plot
%ISOCLINE STUFF 
plot([0 0],[-1 1.2],'b'); hold on
xx = linspace(0,1,101);
plot(xx,1-xx,'b');
% y-isoclines (red)
plot([-1 1.2],[0 0],'r');
plot([1 1]*gamma,[-1 1.2],'r');
%xlim([-0.02 1.1]), ylim([-0.02 1.1])  %establish axis boundaries  
set(gca,'fontsize',16)
xlabel('Trees'), ylabel('Policy')              %establish axis title
grid on

plot(xVector,gammaVector,'k');        % plot points

% THIS IS THE STARTING PONT 
plot(xVector(1),gammaVector(1),'ok','markersize',8);

% THIS IS THE END POINT
plot(xVector(end),gammaVector(end),'pk','markersize',10,'markerfacecolor','k');

%-EOF- 
