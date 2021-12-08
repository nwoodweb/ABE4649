clear, close all, format compact
% ABE4649 FINAL GROUP PROJECT 
% STEVEN COLL, KELSEY VOUGHT, NATHAN WOOD
% FALL 2021 
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

time = 2000;                % establish time 
dt = .001;                  % establish time interval
tVector = [1:dt:time]';      % create vector 1-time @ dt increment

xVector = zeros(size(tVector));     % initialize x vector
gammaVector = zeros(size(tVector)); % initialize y vector 

beetleMortality = .15;       % natural beetle death rate 1/time
rFumigation = 0.3;           % fumigation rate 1/time
rBeetle = .35;               % intrinsic beetle 1/time
rHuman = .3;               % intrinsic policy 1/ttime
rTree = .09;                 % intrinisc tree growth 1/time
beetleDamage =  .24;         % beetle damage to tree 1/(beetle * tree)
tree = 900;                   % Trees tree curse the eigentrees
c = .4;                     % cost of fumigation $/beetle
q = 12;                      % num. fumigations dimless
f = 1;                      
P = 55;                      % profitability $/tree
h = .21;                     % harvest effort 1/$
K = 1000;                      % tree carry capacity tree
z = 200;                    % habitability for beetle 1/tree
H = 6;                     % policy effort $ 
gVal = 300;                                         % beetles per dollar

% BEGIN PREClOMPUTATIONS 

beetleEq = (h*rHuman*P*K)/(rTree)               % nondim unit beta 
theta = (rHuman*P*h*K)/(rTree)                  % nondim unit theta
phi = (c*q*rHuman*h*f)/(rTree*beetleDamage)     % nondim unit phi 
alfa = (rFumigation*beetleDamage * gVal)/(h*f*z*rBeetle*K)          % nondim unit alfa
gamma = ((c*q*rFumigation*rHuman)/(z*rBeetle*K*rTree))                      % nondim unit gamma
y = (H*h*f)/(rTree)
s = (beetleMortality)/(z * rBeetle * K) 

% INITIAL CONDITIONS FOR EULER 
xVector(1) = tree;                               % init x condition
gammaVector(1) = H;                         % init gamma cond. 

% BEGIN EULER 
for t = 1:(length(tVector )-1)               
    % for each incremental unit time
    %    iterate through the provided nondim eqn to populate xVector
    %    and gammaVector 
    %    y(n+1) = y(n) + delta*(f(n+1) : n > 0 
    a = alfa;                               
    th = theta;         
    p = phi;
    b = beetleEq;
    x = xVector(t);
    g = gammaVector(t);
    xVector(t+1) = xVector(t) + dt*(x*(1-x) - ((a*x*y)/(x-s)) - y*x);      %dTree
    gammaVector(t+1) = gammaVector(t) + dt*(b*x*y) - ((gamma*y)/(x-f)); %dPolicy
end % END FOR LOOP

%PLOTTING SYSTEM 

figure(101)                           %initialize plot
%ISOCLINE STUFF 
% TO COMPUTE ISOCLINE 1-x-((alfa*y)/(x-delta)) -y = 0 
% WE NEED TO COMPUTE EULER
eulerClineVectorOne = zeros(size(tVector));
eulerClineVectorOne(1) = 0;                 % initial eulerCline value 
% BEGEIN EULER FOR ISOCLINE 
for t = 1:(length(tVector )-1)               
    % for each iteration in t 
    %   iterate through provided isocline eqn to ppulate
    %   eulerClineVectorOne 1 - x - (*y)/(x-delta)) - y = 0 
    a = alfa; 
    x = xVector(t); 
    eulerClineVectorOne(t+1) = eulerClineVectorOne(t) + dt*(((1-x)*(x-s))/(x-s+alfa)); 
end % end isocline euler 

plot(xVector,eulerClineVectorOne,'b'); hold on

%WE NEED TO PLOT ROOTS OF OTHER ISOCLINE
%    beetleEq*xRoot^2 - beetleEq*gamma*xRoot - gamma = 0 
%    accomplish with quadratic eqn
%    -b +- sqrt(b^2 - 4ac) / 2a 
b = beetleEq
xRootOne = ((b* gamma) + sqrt(((b*gamma)^2)+(4*b*(1*gamma)))) / b
xRootTwo = ((b * gamma) - sqrt(((b*gamma)^2)+(4*b*(1*gamma)))) / b


set(gca,'fontsize',14)
title('Beetlemania')
xlabel('Trees'), ylabel('Profit')              %establish axis title
grid on

plot(xVector,gammaVector,'k');                 % plot points
% THIS IS THE STARTING PONT 
plot(xVector(1),gammaVector(1),'ok','markersize',8);
% THIS IS THE END POINT
plot(xVector(end),gammaVector(end),'o','markersize',10,'markerfacecolor','k');
plot(0,xRootOne,'r');
plot(xRootOne,'r');
%---EOF--- 
