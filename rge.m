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
mt = 172.44;         % Top quark mass
mb = 4.18;           % Bottom quark mass
mh = 125.09;         % Higgs mass
mtau = 1.777;        % Tau mass
v = 246;             % SM Higgs VEV
vS = 1e11;           % Scalar singlet VEV
mn = 0.2*vS;         % Heavy neutrino mass
yt0 = mt*sqrt(2)/v;  % Top quark Yukawa coupling
yb0 = mb*sqrt(2)/v;  % Bottom quark Yukawa coupling
ytau0 = mtau*sqrt(2)/v;  % Tau Yukawa coupling
yn0 = mn*sqrt(2)/vS; % Majorana neutrino Yukawa coupling
yf0 = 0.05;           % Dirac neutrino Yukawa coupling
yq0 = 0.1;           % New quark Yukawa coupling
% Scalar potential parameters
lambdaH0 = mh^2/(2*v^2); % SM Higgs self-coupling at MZ
lambdaS0 = 1e-12;            % Scalar singlet self-coupling
lambdaHS0 = 2e-5;           % Scalar singlet-doublet coupling
muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter

% Coupling constants et cetera
g10 = 0.357;         % U(1) gauge coupling at MZ
g20 = 0.652;         % SU(2) gauge coupling at MZ
g30 = 1.221;         % SU(3) gauge coupling at MZ
q = -1/3;           % Extra quark-like field charge in e-units
nf = 6;             % Number of flavors
ng = nf/2;          % Number of generations

% Constants from Phys. Rev. D 46.3945
b1 = -4/3*ng-1/10;
b2 = 22/3-4/3*ng-1/6;
b3 = 11-4/3*ng;
b = diag([0 136/3 102]) - ng*[[19/5 1/5 11/30];[3/5 49/3 3/2];[44/15 4 76/3]] - [[9/50 3/10 0];[9/10 13/6 0];[0 0 0]];
C = [[17/10 1/2 3/2];[3/2 3/2 1/2];[2 2 0]];

% THIS IS WHERE THE ACTION BEGINS
% Initial values (all)
x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
% Initial value (no scalar masses)
%x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 yq0 lambdaS0 lambdaHS0 yn0];

% Energy scale: 10^Escale GeV, log10(MZ/GeV) = 1.96
Escale = [1.96 19]; 
% Relative and absolute error tolerances
%opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
[t, x] = ode45('rgeq',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
% All parameters
g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
yf = x(:,7);  lambdaH = x(:,8);           muH = sqrt(x(:,9));        yq = x(:,10);
lambdaS = x(:,11);  lambdaHS = x(:,12);   muS = sqrt(x(:,13));       yn = x(:,14);

% All parameters except scalar masses
%  g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
%  yf = x(:,7);  lambdaH = x(:,8);           yq = x(:,9);
%  lambdaS = x(:,10);  lambdaHS = x(:,11);   yn = x(:,12);

figure;
subplot(2,3,1);
plot(t,g1); hold on; plot(t,g2); plot(t,g3); 
ylim([0 2.5]); xlabel('log_{10} \mu/GeV');
legend('{\fontsize{15}g_1}','{\fontsize{15}g_2}','{\fontsize{15}g_3}','Location','NorthEast');
subplot(2,3,2);
plot(t,yt); hold on; plot(t, yn); xlabel('log_{10} \mu/GeV');
legend('{\fontsize{15}y_t}','{\fontsize{15}y_N}','Location','NorthEast');
subplot(2,3,3);
plot(t, yb); hold on; plot(t, ytau); plot(t, yf); plot(t, yq); xlabel('log_{10} \mu/GeV');
legend('{\fontsize{15}y_b}','{\fontsize{15}y_\tau}','{\fontsize{15}y_F}','{\fontsize{15}y_Q}','Location','NorthEast');
subplot(2,3,4);
plot(t, lambdaH); hold on; plot(t, 10^2*lambdaS); plot(t, 10^4*lambdaHS);
xlabel('log_{10} \mu/GeV'); % ylim([0 1.0]); 
legend('{\fontsize{15}\lambda_H}','{\fontsize{15}100\lambda_S}','{\fontsize{15} 10^4\lambda_{HS}}','Location','NorthEast');
subplot(2,3,5);
plot(t, muH); hold on; plot(t, muS/10^2);
ylabel('GeV'); xlabel('log_{10} \mu/GeV');
legend('{\fontsize{15}m_H}','{\fontsize{15}0.01m_S}','Location','NorthEast');
h = suptitle(['5m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0), ', Y_F = ', num2str(yf0) , ', Y_Q = ', num2str(yq0)]);
set(h,'FontSize',20,'FontWeight','bold');