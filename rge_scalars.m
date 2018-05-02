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
prefix = 'May2_4-';
limit = -1;
lambdaOnly = false;
range = false;
% Masses, VEVs and Yukawa couplings at energy scale MZ
mt = 172.44;         % Top quark mass 172.44
mb = 4.18;           % Bottom quark mass
mh = 125.09;         % Higgs mass 125.09
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
yf0 = 1e-6;             % Dirac neutrino Yukawa coupling 1e-3
yq0 = yf0;             % New quark Yukawa coupling 1e-3
vS = 1e8;               % Scalar singlet VEV 2e9
mn = 1/nscale*vS;       % Heavy neutrino mass 1/nscale*vs
yn0 = mn*sqrt(2)/vS;    % Majorana neutrino Yukawa coupling mn*sqrt(2)/vS

% Scalar potential parameters
lambdaH0 = mh^2/(v^2);  % SM Higgs self-coupling at MZ
lambdaS0 = 5e-9;        % Scalar singlet self-coupling 5e-9
lambdaHS0 = 1e-6;       % Scalar singlet-doublet coupling 7e-6
muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter
if(range)
%    lhsRange = 0:0.005:0.5;
    lhsRange = -10:0.2:-0.4;
    vsRange = 8;
else
   lhsRange = log10(lambdaHS0);
   vsRange = 8:0.2:18;
end
limits = zeros(1,length(lhsRange));
n=0;
k=0;
for lambdaHS0 = 10.^lhsRange
    if range
       k = k+1; 
    end
    for vS = 10.^vsRange
        if ~range
            k = k+1;
        end
        mn = 1/nscale*vS;       % Heavy neutrino mass 1/nscale*vs
        yn0 = mn*sqrt(2)/vS;    % Majorana neutrino Yukawa coupling mn*sqrt(2)/vS
        muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter

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
        end
        if(~lambdaOnly)
    %        subplot(1,2,1);
            plot(t, log10(muH),'LineWidth',2); %hold on; plot(t, log(muS));
            xlim([2 20]); ylim([0 20]); 
            set(gca,'XMinorTick','on','YMinorTick','on');
            set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
            set(gca,'FontSize',15);
            xlabel('log_{10} \mu/GeV');
            ylabel('log_{10} m/GeV');
            grid on;
            title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
            hold on;
            plot([2 20],log10([mh 151]),'LineWidth',2);
            plot(t, log10(muS),'LineWidth',2);
            legend('{\fontsize{15}SMASH m_H}','{\fontsize{15}SM m_H}','{\fontsize{15} m_\sigma}','Location','NorthWest');
    %         subplot(1,2,2);
    %         xlabel('log_{10} \mu/GeV');
    %         plot(t, lambdaH); hold on; plot(t, 10^7*lambdaS); plot(t, 10^4*lambdaHS);
    %         ylim([-0.5 0.5]);
    %         h = suptitle([num2str(nscale), 'm_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0), ', Y_F = ', num2str(yf0) , ', Y_Q = ', num2str(yq0)]);
    %         set(h,'FontSize',15,'FontWeight','bold');
    %     %    title([num2str(100/(2*k),3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, \lambda_S = ', num2str(lambdaS0), ', \lambda_{HS} = ', num2str(lambdaHS0)],'FontSize',20);
    %         legend('{\fontsize{15}\lambda_H}','{\fontsize{15}10^7\lambda_S}','{\fontsize{15} 10^4\lambda_{HS}}','Location','NorthEast');
        else
            for k = 1:length(lambdaH)
                if(lambdaH(k) < 0 || lambdaH(k) > 1)
                    limits(n) = t(k);
                    limit = t(k);
                    break;
                end
            end
            plot (t, lambdaH); hold on;
            if(limit > 0)
                vline(limit,'r','{\fontsize{20}Stability bound}');
            end
            ylim([-0.1, max(lambdaH)]);
            set(gca,'XMinorTick','on','YMinorTick','on');
            set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
            set(gca,'FontSize',15);
            grid on;
            xlabel('log_{10} \mu/GeV');
            title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
            legend('{\fontsize{15}\lambda_H}','Location','NorthWest');
    %        fprintf('Limit = %f\n', limit);
        end
        if(~debug)
            if k < 10
                filename = sprintf('%s0%d.png',prefix,k);
            else
                filename = sprintf('%s%d.png',prefix,k);
            end
            saveas(gcf,filename);
        end
    end
end

if(lambdaOnly)
    nPoints = length(limits);
    for k = 1:length(limits)
       if(limits(k) == 0)
           nPoints = k;
           limits(k) = 20; 
           break;
       end
    end
    lhsRange = lhsRange(1:nPoints);
    limits = limits(1:nPoints);
    %plot(lhsRange, limits,'LineWidth',1); hold on;

    global xdata;
    global ydata;
    xdata = lhsRange;
    ydata= limits;
    lam0=[-1.8 1 10];         % Initial guess for lambda
    y0=fobj(lam0);        % Initial value of object function
    lam=fminsearch('fobj',lam0);
                          % lam is the fitted value for
                          % the parameter vector
    x=0:0.005:0.25;
    %x=(0:0.005:0.5)*1e-5;
    yfit=fmodel(lam, x);
    yfinal=fobj(lam);     % Final value of the object function

    % Fit y = Ae^Bx+C
    %P = polyfit(lhsRange, limits, 5);
    %fittedLimits = P(1)*lhsRange.^5 + P(2)*lhsRange.^4 + P(3)*lhsRange.^3 + P(4)*lhsRange.^2 + P(5)*lhsRange + P(6);
    plot(x, yfit,'LineWidth',2);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
    set(gca,'FontSize',15);
    grid on;
    xlabel('\lambda_{H\sigma}');
    ylabel('log_{10} \mu/GeV');
    ylim([12 19]);
    title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3)],'FontSize',15);
end