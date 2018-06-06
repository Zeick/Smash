% SMASH-RGE project
% (C) Timo K?rkk?inen 2017-2018
%%%%%%%%%%%%%%%%%%
% CONSTANTS      %
%%%%%%%%%%%%%%%%%%
% M = Y*v/sqrt(2), using GeV units
% v = mu/sqrt(lambda)
clear;
timeElapsed = cputime;
fs = 15;
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
yq0 = 0.001;              % New quark Yukawa coupling 1e-6
% Scalar potential parameters
lambdaH0 = mh^2/(v^2);   % SM Higgs self-coupling at MZ
lambdaHS0 = -1e8*(v/vS)^2;       % Scalar singlet-doublet coupling 1e-10, -0.34, -sqrt(0.05*lambdaS0)
lambdaS0 = 1e5*lambdaHS0^2;         % Scalar singlet self-coupling 5e-9

muH0 = mh;
muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter

% Coupling constants et cetera
g10 = 0.357;         % U(1) gauge coupling at MZ
g20 = 0.652;         % SU(2) gauge coupling at MZ
g30 = 1.221;         % SU(3) gauge coupling at MZ

% THIS IS WHERE THE ACTION BEGINS
% Initial values (all)
% Energy scale: 10^Escale GeV, log10(MZ/GeV) = 1.96
lsRange = -11:0.1:0;
lhsRange = -11:0.1:0;
warning off;
for lambdaHS0 = -10.^lhsRange
    fprintf('lhs = %g\n',lambdaHS0);
    for lambdaS0 = 10.^lsRange
        muS0 = vS*sqrt(lambdaS0); % Scalar singlet mu parameter
        Escale = [1.96 log10(sqrt(2)*muS0)]; % Planck energy = 1.22e19 GeV 
        x0 = [g10 g20 g30 yt0 yb0 ytau0 lambdaH0 muH0^2];
        [t1, x] = ode45('rgeq_SM',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
        % All parameters
        g1 = x(:,1);  g2 = x(:,2);  g3 = x(:,3);
        yt = x(:,4);  yb = x(:,5);  ytau = x(:,6);
        lambdaH = x(:,7);           muH = sqrt(x(:,8));
        g10 = x(end,1);  g20 = x(end,2);  g30 = x(end,3);
        yt0 = x(end,4);  yb0 = x(end,5);  ytau0 = x(end,6);
        lambdaH0 = x(end,7); muH0 = sqrt(x(end,8));

        lambdaH0 = lambdaH0 + lambdaHS0^2/lambdaS0;
%        muH0 = sqrt(muH0^2 + lambdaHS0*vS^2); NO THRESHOLD CORRECTION!
        Escale = [log10(sqrt(2)*muS0) log10(1.22e19)];
        x0 = [g10 g20 g30 yt0 yb0 ytau0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];
        [t2, x] = ode45('rgeq',Escale,x0); % ODE-function, solution span, initial values, (opt.) error tolerance
        % All parameters
        g1 = [g1' x(:,1)'];  g2 = [g2' x(:,2)'];  g3 = [g3' x(:,3)'];
        yt = [yt' x(:,4)'];  yb = [yb' x(:,5)'];  ytau = [ytau' x(:,6)'];
        yf = x(:,7);  lambdaH = [lambdaH' x(:,8)']; muH = [muH' sqrt(x(:,9))'];        yq = x(:,10);
        lambdaS = x(:,11);  lambdaHS = x(:,12);   muS = sqrt(x(:,13));       yn = x(:,14);

        if(lambdaH(end) > 0 && lambdaH(end) < 1 && ~isnan(lambdaH(end)) && real(muH(end)) > 0)
           fprintf('Solution: lam_HS = %g, lam_S = %g\n',lambdaHS0, lambdaS0); 
        end
    end
end