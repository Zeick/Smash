% SMASH-RGE project
% (C) Timo Karkkainen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
timeElapsed = cputime;
debug = true;        % False = save data, no output in console
lambdaOnly = true;
prefix = 'April15_test';
% Masses, VEVs and Yukawa couplings at energy scale MZ
mt = 172.44;         % Top quark mass
mb = 4.18;           % Bottom quark mass
mh = 125.09;         % Higgs mass
mtau = 1.777;        % Tau mass
v = 246;             % SM Higgs VEV
yt0 = mt*sqrt(2)/v;  % Top quark Yukawa coupling
yb0 = mb*sqrt(2)/v;  % Bottom quark Yukawa coupling
ytau0 = mtau*sqrt(2)/v;  % Tau Yukawa coupling

% Coupling constants et cetera
g10 = 0.357;         % U(1) gauge coupling at MZ
g20 = 0.652;         % SU(2) gauge coupling at MZ
g30 = 1.221;         % SU(3) gauge coupling at MZ

% Energy scale: 10^Escale GeV, log10(MZ/GeV) = 1.96,
% log10(1.22e19/GeV) = 19.09
Escale = [1.96 19.09];
nscale = 100; % 100
yf0 = 1e-3;             % Dirac neutrino Yukawa coupling
yq0 = yf0;              % New quark Yukawa coupling
vS = 1e10;               % Scalar singlet VEV 2e9
mn = 1/nscale*vS;       % Heavy neutrino mass 1/nscale*vS
yn0 = mn*sqrt(2)/vS;    % Majorana neutrino Yukawa coupling mn*sqrt(2)/vS

% Scalar potential parameters
lambdaH0 = mh^2/(v^2);% SM Higgs self-coupling at MZ
lambdaS0 = 5e-9;        % Scalar singlet self-coupling 5e-9
lambdaHS0 = 0.39;       % Scalar singlet-doublet coupling 7e-6
muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter
mtRange = 164:0.1:182;
mhRange = 110:0.1:140;
%mtRange = mt; mhRange = mh;

mtIndex = 0;
mhLimit = zeros(1,length(mtRange));
stable = true;
for mt = mtRange
    fprintf('mt = %.2f\n',mt);
    mtIndex = mtIndex + 1;
    mhLimit(mtIndex) = 200;
    mhIndex = 0;
    stable = true;
    for mh = mhRange
        if(~stable)
            break
        end
        mhIndex = mhIndex + 1;
        yt0 = mt*sqrt(2)/v;
        muH0 = mh;
        lambdaH0 = mh^2/(v^2);
        % THIS IS WHERE THE ACTION BEGINS
        % Initial values (all)
        x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
        [t, x] = ode45('rgeq',Escale,x0); % ODE-function, solution span, init-values, (opt.) error tolerance
        % All parameters
        g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5); ytau = x(:,6);
        yf = x(:,7);        lambdaH = x(:,8);     muH = sqrt(x(:,9));        yq = x(:,10);
        lambdaS = x(:,11);  lambdaHS = x(:,12);   muS = sqrt(x(:,13));       yn = x(:,14);
        % Check if the potential is stable, ie. if quartic couplings are
        % positive
       for k = 1:length(lambdaH)
           if(lambdaHS(k) > 1)
               stable = false;
               mhLimit(mtIndex) = mh;
               break;
           end
       end
   end
end
plot(mtRange,mhLimit,'LineWidth',2); hold on;

% SM area
mt = 172.44;         % Top quark mass
mh = 125.09;
dmt = 0.60;     dmh = 0.32;
mtLower = mt-dmt;  mtUpper = mt+dmt;
mhLower = mh-dmh;  mhUpper = mh+dmh;
h = rectangle;
h.Position = [mtLower mhLower 2*dmt 2*dmh];
h.FaceColor = 'magenta';

xlim([min(mtRange) max(mtRange)]); xlabel('m_t','FontSize', 20);
ylim([min(mhRange) max(mhRange)]); ylabel('m_H','FontSize', 20);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'FontSize',15);
text(172,115,'Non-perturbative region','FontSize',20);
text(166,125,'\lambda_{H\sigma} stable','FontSize',20);
grid on;
title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);



%legend('{\fontsize{15}\lambda_H}','Location','NorthEast');
% fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);