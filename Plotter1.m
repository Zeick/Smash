% Helper to plot data
clear;
S = load('StabilityMatrix.mat');
S = S.stabilityLimit;
mhRange = 110:5:140;
mtRange = 164:2:182;
nscale = 100; % 100
yf0 = 1e-3;             % Dirac neutrino Yukawa coupling
vS = 1e10;               % Scalar singlet VEV 2e9
lambdaS0 = 5e-9;        % Scalar singlet self-coupling 5e-9
lambdaHS0 = 0.5;       % Scalar singlet-doublet coupling 7e-6
mt = 172.44;    mh = 125.09;
dmt = 0.60;     dmh = 0.32;
mtLower = mt-dmt;  mtUpper = mt+dmt;
mhLower = mh-dmh;  mhUpper = mh+dmh;
figure;
cLevels = 6:18;
[cc hh] = contour(mhRange, mtRange, S,[18 18],'ShowText','Off'); hold on;
xlabel('m_H/GeV'); ylabel('m_t/GeV');
title([num2str(nscale,3),'m_N = v_\sigma = ', num2str(vS,3), ' GeV, Y_{F/Q} = ', num2str(yf0,3), ', \lambda_{\sigma} = ', num2str(lambdaS0,3), ', \lambda_{H\sigma} = ', num2str(lambdaHS0,3)],'FontSize',15);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'FontSize',15);
grid on;
xlim([110 140]);
ylim([164 182]);

% Fill with color
cTop = cc(2,2:12);      cTop2 = cc(2,14:end);
cHiggs = cc(1,2:12);    cHiggs2 = cc(1,14:end);
area(cHiggs2, 182*ones(length(cHiggs2)), 'FaceColor', 'Yellow', 'LineStyle', 'none');
area(cHiggs2, cTop2,'FaceColor','Green','LineStyle','none');
area([cHiggs2(1) cHiggs(end)], mtRange(end)*[1 1],'FaceColor','Green','LineStyle','none');
area([cHiggs mhRange(end)], [cTop mtRange(end)],'FaceColor','Cyan','LineStyle','none');

% Experimental range in (mHiggs, mTop) plane
h = rectangle;
h.Position = [mhLower mtLower 2*dmh 2*dmt];
h.FaceColor = 'magenta';
plot(mh,mt,'k.','MarkerSize',20);
text(135,166,'Unstable region','FontSize',20,'Rotation',90);
text(120,174','Stability region','FontSize',20,'Rotation',45);
text(111,177.2','Metastability region','FontSize',18,'Rotation',28);
