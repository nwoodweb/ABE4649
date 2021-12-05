clear, close all, format compact, clc
% ABE4649 FINAL GROUP PROJECT 
% 
% The purpose of this project is to model the 3 dimensional "eigenforest"
% Model created by our group. 
% The equations detail change in tree population (dx)
%       change in parasitic beetle population (db) 
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
%********* END PARAMETERS

%SETUP PARAMETERS

time = 1500;                % establish time
dt = .001;                  % establish time interval
tVector = [1:dt:time]';      % create vector 1-time @ dt increment

xVector = zeros(size(tVector));
gammaVector = zeros(size(tVector));

beetleMortality = .2;        % natural beetle death rate 
rFumigation = .1;            % setup fumigation rate 
rBeetle = .3;
rHuman = .12;     
rTree = .5;
beetleDamage = .12;

tree = 1;
c = 24;
q = 8;
f = 1;
P = 4;
h = 23; 
K = 5; 
z = 11;
H = .12; 

beetleEq = (rFumigation*H)/((z*rBeetle*tree) - (beetleMortality))

theta = (rHuman*P*h*K)/(rTree)
phi = (c*q*rHuman*h*f)/(rTree*beetleDamage)
alfa = (beetleEq*beetleDamage)/(rTree)
gamma = (h*H*f)/(rTree)

xVector(1) = 1.2;           % init x condition
gammaVector(1) = 10.2;
length(tVector)

% BEGIN EULER 
for t = 1:(length(tVector)-1)
    a = alfa;
    th = theta;
    p = phi;
    x = xVector(t);
    g = gammaVector(t);
    xVector(t+1) = xVector(t) + dt*(x*(1-x) - a*x - gamma*x);
    gammaVector(t+1) = gammaVector(t) + dt*((th*gamma*x) - (p*x));
end 

%BEGIN PLOT 
figure(101)                           %initialize plot
%ISOCLINE STUFF 
plot([0 0],[-1 1.2],'b'); hold on
xx = linspace(0,1,101);
plot(xx,1-xx,'b');
% y-isoclines (red)
plot([-1 1.2],[0 0],'r');
plot([1 1]*gamma,[-1 1.2],'r');
%xlim([-0.02 1.1]), ylim([-0.02 1.1])  %establish axis bound        
set(gca,'fontsize',16)
xlabel('x'), ylabel('y')              %establish axis title
grid on
% Phase portrait avnd endpoints (black)
plot(xVector,gammaVector,'k');        % plot points
plot(xVector(1),gammaVector(1),'ok','markersize',8);
plot(xVector(end),gammaVector(end),'pk','markersize',10,'markerfacecolor','k');

%-EOF- 