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
green = [0 1 0];
cTop = cc(2,2:12);      cTop2 = cc(2,14:end);
cHiggs = cc(1,2:12);    cHiggs2 = cc(1,14:end);
area(cHiggs2, 182*ones(length(cHiggs2)), 'FaceColor', 'Yellow', 'LineStyle', 'none');
area(cHiggs2, cTop2,'FaceColor',green,'LineStyle','none');
area([cHiggs2(1) cHiggs(end)], mtRange(end)*[1 1],'FaceColor',green,'LineStyle','none');
area([cHiggs mhRange(end)], [cTop mtRange(end)],'FaceColor','Cyan','LineStyle','none');
alfa = area(mhRange, 182*ones(length(mhRange)),'FaceColor', [0.91 0.41 0.17], 'LineStyle', 'none'); 
set(alfa,'facealpha',.75)

% Experimental range in (mHiggs, mTop) plane
h = rectangle;
h.Position = [mhLower mtLower 2*dmh 2*dmt];
h.FaceColor = 'magenta';
plot(mh,mt,'k.','MarkerSize',20);
fs = 20;
text(112,168,'Excluded due to \lambda_{H\sigma}','FontSize',fs+5,'Color','Red');
text(112,167,'non-perturbativity','FontSize',fs+5,'Color','Red');
text(135,165,'Non-perturbativity','FontSize',fs,'Rotation',90);
text(136.5,165,'region','FontSize',fs,'Rotation',90);
text(117,171','Stability region','FontSize',fs,'Rotation',50);
text(111,177.45','Metastability region','FontSize',fs-2,'Rotation',27);
