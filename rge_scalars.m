% SMASH-RGE project
% (C) Timo K?rkk?inen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
timeElapsed = cputime;
debug = false;        % False = save data, no output in console
prefix = 'April10_yF';
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
yf0 = 1e-3;             % Dirac neutrino Yukawa coupling
yq0 = 1e-3;             % New quark Yukawa coupling
vS = 2e9;               % Scalar singlet VEV 2e9
mn = 1/nscale*vS;       % Heavy neutrino mass
yn0 = mn*sqrt(2)/vS;    % Majorana neutrino Yukawa coupling

% Scalar potential parameters
lambdaH0 = mh^2/(2*v^2);% SM Higgs self-coupling at MZ
lambdaS0 = 5e-9;        % Scalar singlet self-coupling 5e-9
lambdaHS0 = 7e-6;       % Scalar singlet-doublet coupling 7e-6
muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter
for k=1:50
    yf0 = k*1e-3;
    yq0 = yf0;
    % THIS IS WHERE THE ACTION BEGINS
    % Initial values (all)
    x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
    % Initial value (no scalar masses)
    %x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 yq0 lambdaS0 lambdaHS0 yn0];

    % Relative and absolute error tolerances
    %opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
    
    [t, x] = ode45('rgeq',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
    % All parameters
    g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);  yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
    yf = x(:,7);  lambdaH = x(:,8);           muH = sqrt(x(:,9));        yq = x(:,10);
    lambdaS = x(:,11);  lambdaHS = x(:,12);   muS = sqrt(x(:,13));       yn = x(:,14);
    if(~debug)
        f = figure('visible', 'off');        
    else
        fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
        fprintf('%.2g\n',lambdaS(end));
    end
    subplot(1,2,1);
    xlabel('log_{10} \mu/GeV');
    plot(t, log(muH)); hold on; plot(t, log(muS));
    ylim([0 19]); 
    ylabel('log_{10} m_i /GeV'); xlabel('log_{10} \mu/GeV');
    legend('{\fontsize{15}m_H}','{\fontsize{15}m_S}','Location','NorthEast');
%    title([num2str(100/(2*k),3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0)],'FontSize',20);
%    title(['5m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0)],'FontSize',20);
%     if k < 10
%         filename = sprintf('%s3_0%d.png','scalars',k);
%     else
%         filename = sprintf('%s3_%d.png','scalars',k);
%     end
%     saveas(gcf,filename)
%     close(f);
    
%     g = figure('visible', 'off');
    subplot(1,2,2);
    xlabel('log_{10} \mu/GeV');
    plot(t, lambdaH); hold on; plot(t, 10^7*lambdaS); plot(t, 10^4*lambdaHS);
    ylim([0 0.5]);
    h = suptitle([num2str(nscale), 'm_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0), ', Y_F = ', num2str(yf0) , ', Y_Q = ', num2str(yq0)]);
    set(h,'FontSize',15,'FontWeight','bold');
%    title([num2str(100/(2*k),3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0)],'FontSize',20);
    legend('{\fontsize{15}\lambda_H}','{\fontsize{15}10^7\lambda_S}','{\fontsize{15} 10^4\lambda_{HS}}','Location','NorthEast');
    if(~debug)
        if k < 10
            filename = sprintf('%s_0%d.png',prefix,k);
        else
            filename = sprintf('%s_%d.png',prefix,k);
        end
        saveas(gcf,filename);
    end
    close(f);
end
if(~debug)
    fprintf('Time elapsed: %.2f seconds.\n', cputime - timeElapsed);
end