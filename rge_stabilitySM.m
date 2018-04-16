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

% Energy scale: 10^Escale GeV, log10(MZ/GeV) = 1.96
Escale = [1.96 19];
nscale = 100; % 100
yf0 = 0;             % Dirac neutrino Yukawa coupling 1e-3
yq0 = yf0;              % New quark Yukawa coupling 1e-3
vS = 0;               % Scalar singlet VEV 2e9
mn = 1/nscale*vS;       % Heavy neutrino mass
yn0 = 0;
%yn0 = mn*sqrt(2)/vS;    % Majorana neutrino Yukawa coupling

% Scalar potential parameters
lambdaH0 = mh^2/(2*v^2);% SM Higgs self-coupling at MZ
lambdaS0 = 0;        % Scalar singlet self-coupling 5e-9
lambdaHS0 = 0;       % Scalar singlet-doublet coupling 7e-6
muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter
mtRange = 164:182;
mhRange = 110:140;
%exprange = -7:0.5:-4; lhsRange = 10.^exprange;
%lhsRange = [-lhsRange lhsRange]; % Take also negative value into account
stabilityLimit = zeros(length(mtRange),length(mhRange));
%stabilityLimit = zeros(length(mtRange),length(mhRange), length(lhsRange), length(lhsRange));
%                               mtop           mhiggs         lambda_HS         lambda_HS
mtIndex = 0;
for mt = mtRange
    mtIndex = mtIndex + 1;
    mhIndex = 0;
    for mh = mhRange
        mhIndex = mhIndex + 1;
        stable = true;
        yt0 = mt*sqrt(2)/v;
        muH0 = mh;
        % THIS IS WHERE THE ACTION BEGINS
        % Initial values (all)
        x0 = [g10 g20 g30 yt0 yb0 ytau0 lambdaH0 muH0^2];
        [t, x] = ode45('rgeq_SM',Escale,x0); % ODE-function, solution span, init-values, (opt.) error tolerance
        % All parameters
        g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5); ytau = x(:,6);
        lambdaH = x(:,7);   muH = sqrt(x(:,8));
        % Check if the potential is stable, ie. if quartic couplings are
        % positive
        for k = 1:length(lambdaH)
           if(lambdaH(k) < 0 || lambdaH(k) > 1)
               stable = false;
               stabilityLimit(mtIndex, mhIndex) = t(k);
               break;
           end
        end

%         f = figure('visible', 'off');
%         xlabel('log_{10} \mu/GeV');
%         plot(t, lambdaH); hold on; plot(t, 10^7*lambdaS); plot(t, 10^4*lambdaHS);
%         ylim([0 0.5]);
%         h = suptitle([num2str(nscale), 'm_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0), ', Y_F = ', num2str(yf0) , ', Y_Q = ', num2str(yq0)]);
%         set(h,'FontSize',15,'FontWeight','bold');
%     %    title([num2str(100/(2*k),3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0)],'FontSize',20);
%         legend('{\fontsize{15}\lambda_H}','{\fontsize{15}10^7\lambda_S}','{\fontsize{15} 10^4\lambda_{HS}}','Location','NorthEast');
%         if(~debug)
%             if k < 10
%                 filename = sprintf('%s_0%d.png',prefix,k);
%             else
%                 filename = sprintf('%s_%d.png',prefix,k);
%             end
%             saveas(gcf,filename);
%         end
%         close(f);
    end
end
figure;
contour(mhRange, mtRange, stabilityLimit,'ShowText','On');
if(~debug)
    fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
end
