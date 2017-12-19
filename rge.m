% SMASH-RGE project
% (C) Timo K?rkk?inen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
timeElapsed = cputime;
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
yf0 = 0.05;           % Dirac neutrino Yukawa coupling
yq0 = 0.1;           % New quark Yukawa coupling
% Scalar potential parameters
lambdaH0 = mh^2/v^2; % SM Higgs self-coupling at MZ
lambdaS0 = 1e-11;            % Scalar singlet self-coupling
lambdaHS0 = 1e-6;           % Scalar singlet-doublet coupling
muH0 = mh;
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
tspan = [1.96 20]; % log10(MZ/GeV) = 1.96
%opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
[t, x] = ode45('rgeq',tspan,x0); % ODE-function, solution span, initial values, error tolerance
fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
%g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
%yf = x(:,7);  lambdaH = x(:,8);           muH0 = sqrt(x(:,9));        yq = x(:,10);
%lambdaS = x(:,11);  lambdaHS = x(:,12);   muS0 = sqrt(x(:,13));       yn = x(:,14);

 g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
 yf = x(:,7);  lambdaH = x(:,8);           yq = x(:,9);
 lambdaS = x(:,10);  lambdaHS = x(:,11);   yn = x(:,12);

figure;
subplot(2,2,1);
plot(t,abs(g1)); hold on; plot(t,abs(g2)); plot(t,abs(g3)); 
ylim([0 2.5]);
legend('{\fontsize{15}g_1}','{\fontsize{15}g_2}','{\fontsize{15}g_3}','Location','NorthEast');
subplot(2,2,2);
plot(t,abs(yt)); hold on; plot(t, abs(yn));
legend('{\fontsize{15}y_t}','{\fontsize{15}y_N}','Location','NorthEast');
subplot(2,2,3);
plot(t, abs(yb)); hold on; plot(t, abs(ytau)); plot(t, abs(yf)); plot(t, abs(yq));
legend('{\fontsize{15}y_b}','{\fontsize{15}y_\tau}','{\fontsize{15}y_F}','{\fontsize{15}y_Q}','Location','NorthEast');
subplot(2,2,4);
plot(t, abs(lambdaH)); hold on; plot(t, abs(lambdaS)); plot(t, 100*abs(lambdaHS)); ylim([0 0.4]);
legend('{\fontsize{15}\lambda_H}','{\fontsize{15}\lambda_S}','{\fontsize{15} 100\lambda_{HS}}','Location','NorthEast');
h = suptitle('m_N = v_\sigma = 10^{14} GeV, \lambda_S = 10^{-6}, \lambda_{HS} = 10^{-11}, Y_F = 0.05, Y_Q = 0.1');
set(h,'FontSize',20,'FontWeight','bold');