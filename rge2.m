% SMASH-RGE project
% (C) Timo K?rkk?inen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
timeElapsed = cputime;
fs = 20;
% Masses, VEVs and Yukawa couplings at energy scale MZ
mt = 172.44;             % Top quark mass
mb = 4.18;               % Bottom quark mass
mh = 125.09;             % Higgs mass
mtau = 1.777;            % Tau mass
v = 246;                 % SM Higgs VEV
vS = 1e9;                % Scalar singlet VEV 1e10
nscale = 25;             % 28 or 100?
mn = vS/nscale;          % Heavy neutrino mass
yt0 = mt*sqrt(2)/v;      % Top quark Yukawa coupling
yb0 = mb*sqrt(2)/v;      % Bottom quark Yukawa coupling
ytau0 = mtau*sqrt(2)/v;  % Tau Yukawa coupling
yn0 = mn*sqrt(2)/vS;     % Majorana neutrino Yukawa coupling
yf0 = 0.001;              % Dirac neutrino Yukawa coupling 1e-6
yq0 = 0.1;              % New quark Yukawa coupling 1e-6
% Scalar potential parameters
lambdaH0 = mh^2/(v^2);   % SM Higgs self-coupling at MZ
lambdaHS0 = -1e-5;       % Scalar singlet-doublet coupling 1e-10, -0.34, -sqrt(0.05*lambdaS0)
delta = 3;
lambdaS0 = 2e-9;         % Scalar singlet self-coupling 5e-9, delta*lambdaHS0^2

muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter

% Coupling constants et cetera
g10 = 0.357;         % U(1) gauge coupling at MZ
g20 = 0.652;         % SU(2) gauge coupling at MZ
g30 = 1.221;         % SU(3) gauge coupling at MZ

% THIS IS WHERE THE ACTION BEGINS
% Initial values (all)
% Energy scale: 10^Escale GeV, log10(MZ/GeV) = 1.96
Escale = [1.96 log10(sqrt(2)*muS0)]; % Planck energy = 1.22e19 GeV 
x0 = [g10 g20 g30 yt0 yb0 ytau0 lambdaH0 muH0^2];
[t1, x] = ode45('rgeq_SM',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
% All parameters
g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);
yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
lambdaH = x(:,7);           muH = sqrt(x(:,8));
g10 = x(end,1);  g20 = x(end,2);  g30 = x(end,3);
yt0 = x(end,4);  yb0 = x(end,5);  ytau0 = x(end,6);
lambdaH0 = x(end,7); muH0 = sqrt(x(end,8));

lambdaH0 = lambdaH0 + lambdaHS0^2/lambdaS0;
%muH0 = sqrt(muH0^2 + lambdaHS0*vS^2); NO THRESHOLD CORRECTION TO THIS
Escale = [log10(sqrt(2)*muS0) log10(1.22e19)];
x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
[t2, x] = ode45('rgeq',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
% All parameters
g1 = [g1' x(:,1)'];  g2 = [g2' x(:,2)'];  g3 = [g3' x(:,3)'];
yt = [yt' x(:,4)'];  yb = [yb' x(:,5)'];  ytau = [ytau' x(:,6)'];
yf = x(:,7);  lambdaH = [lambdaH' x(:,8)']; muH = [muH' sqrt(x(:,9))'];        yq = x(:,10);
lambdaS = x(:,11);  lambdaHS = x(:,12);   muS = sqrt(x(:,13));       yn = x(:,14);

% figure;
% plot([t1' t2'],yt,'LineWidth',2); hold on; plot(t2, yn*100,'LineWidth',2); xlabel('log_{10} \mu/GeV');
% plot([t1' t2'], 10*yb,'LineWidth',2); plot([t1' t2'], 10*ytau,'LineWidth',2); plot(t2, yf,'LineWidth',2); plot(t2, 10*yq,'LineWidth',2); xlabel('log_{10} \mu/GeV');
% legend('{\fontsize{15}y_t}','{\fontsize{15}100y_N}','{\fontsize{15}10y_b}','{\fontsize{15}10y_\tau}','{\fontsize{15}y_F}','{\fontsize{15}10y_Q}','Location','NorthWest');
% h = suptitle([num2str(nscale,3),'m_N', '= v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0), ', Y_F = ', num2str(yf0) , ', Y_Q = ', num2str(yq0)]);
% set(h,'FontSize',20,'FontWeight','bold');
% set(gca,'XMinorTick','on','YMinorTick','on');
% set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
% set(gca,'FontSize',15);
% grid on;
% ylim([0 2]);

smcase = load('SMrunning2');
SMmuH = smcase.muH;
SMlambda = smcase.lambdaH;
SMt = smcase.t;

figure;
subplot(1,2,1);
plot([t1' t2'], muH,'LineWidth',2); hold on; plot(t2,muS/10,'LineWidth',2);
plot(SMt, SMmuH,'--b','LineWidth',2);
xlabel('log_{10} \mu/GeV'); ylabel('GeV');
legend('{\fontsize{15}m_H}','{\fontsize{15}m_\sigma/10}','{\fontsize{15}m_H (SM)}','Location','NorthWest');
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'FontSize',15);
ylim([0 200]);
%title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F} = ', num2str(yf0,3), ', Y_{Q} = ', num2str(yq0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
grid on;
fprintf('mH(Planck) = %.2f GeV\n',muH(end));

%figure;
subplot(1,2,2);
plot([t1' t2'], lambdaH,'LineWidth',2); hold on; plot(t2, 10^3*lambdaHS,'LineWidth',2);
plot(SMt, SMlambda, '--b','LineWidth',2);
xlabel('log_{10} \mu/GeV'); % ylim([0 1.0]);
legend('{\fontsize{15}\lambda_H}','{\fontsize{15} 10^3\lambda_{HS}}','{\fontsize{15}\lambda_H (SM)}','Location','NorthEast');
%title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F} = ', num2str(yf0,3), ', Y_{Q} = ', num2str(yq0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'FontSize',15);
grid on;
h = suptitle([num2str(nscale,3),'m_N', '= v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0,3), ', \lambda_{HS} = ', num2str(lambdaHS0,3), ', Y_F = ', num2str(yf0,3) , ', Y_Q = ', num2str(yq0,3)]);
set(h,'FontSize',fs,'FontWeight','bold');
set(gcf, 'Position', [0 0 800 500])
fprintf('lambdaH(Planck) = %.2g\n',lambdaH(end));

% figure;
% plot(t, g1,'LineWidth',2); hold on; plot(t, g2,'LineWidth',2); plot(t, g3, 'LineWidth', 2); 
% %ylim([0 0.025]);
% xlim([2 19]);
% xlabel('log_{10} \mu/GeV');
% legend('{\fontsize{15}g_1}','{\fontsize{15}g_2}','{\fontsize{15}g_3}','Location','NorthEast');
% set(gca,'XMinorTick','on','YMinorTick','on');
% set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
% set(gca,'FontSize',15);
% title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
% grid on;

% figure;
% plot(t, yb,'LineWidth',2); hold on; plot(t, ytau,'LineWidth',2); 
% ylim([0 0.025]); xlim([2 19]);
% xlabel('log_{10} \mu/GeV');
% legend('{\fontsize{15}y_b}','{\fontsize{15}y_\tau}','Location','NorthEast');
% set(gca,'XMinorTick','on','YMinorTick','on');
% set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
% set(gca,'FontSize',15);
% title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15); grid on;
