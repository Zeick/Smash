% SMASH-RGE project
% (C) Timo K?rkk?inen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
% Masses, VEVs and Yukawa couplings at energy scale MZ
mt = 172.44;        % Top quark mass
mb = 4.18;          % Bottom quark mass
mh = 125.09;        % Higgs mass
mtau = 1.777;       % Tau mass
mn = 1e14;          % Heavy neutrino mass
v = 246;            % SM Higgs VEV
vS = 1e14;          % Scalar singlet VEV
yt0 = mt*sqrt(2)/v;  % Top quark Yukawa coupling
yb0 = mb*sqrt(2)/v;  % Bottom quark Yukawa coupling
ytau0 = mtau*sqrt(2)/v;  % Tau Yukawa coupling
yn0 = mn*sqrt(2)/vS; % Majorana neutrino Yukawa coupling
yf0 = 1;
yq0 = 1;
% Scalar potential parameters
lambdaH0 = mh^2/(2*v^2); % SM Higgs self-coupling at MZ
lambdaS0 = 0;            % Scalar singlet self-coupling
lambdaHS0 = 0;           % Scalar singlet-doublet coupling
muH0 = v*sqrt(lambdaH0);  % SM Higgs mu parameter
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter

% Coupling constants et cetera
g10 = 0.357;         % U(1) gauge coupling at MZ
g20 = 0.652;         % SU(2) gauge coupling at MZ
g30 = 1.221;         % SU(3) gauge coupling at MZ
q = -1/3;           % Extra quark-like field charge in e-units
nf = 6;             % Number of flavors
ng = nf/2;

% Constants from Phys. Rev. D 46.3945
b1 = -4/3*ng-1/10;
b2 = 22/3-4/3*ng-1/6;
b3 = 11-4/3*ng;
b = diag([0 136/3 102]) - ng*[[19/5 1/5 11/30];[3/5 49/3 3/2];[44/15 4 76/3]] - [[9/50 3/10 0];[9/10 13/6 0];[0 0 0]];
C = [[17/10 1/2 3/2];[3/2 3/2 1/2];[2 2 0]];

% Initial values
%x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 yq0 lambdaS0 lambdaHS0 yn0];
tspan = [0 20];
[t, x] = ode45('rgeq',tspan,x0);
% plot(t,x);