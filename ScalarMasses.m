% RGE-SMASH-project
% A small program to calculate the masses of SMASH scalars.
% (C) Timo K?rkk?inen
mh = 125.08;          % SM Higgs mass
v = 246;              % SM Higgs VEV
vs = 1e11;            % Scalar singlet VEV 2e9
lh = mh^2/(v^2);% SM Higgs self-coupling at MZ
ls = 5e-9;      % Scalar singlet self-coupling 5e-9
lhs = 0.39;     % Portal coupling
lhsRange = 10.^(-10:0.05:0);
lsRange = lhsRange;
vsRange = 10.^(8:0.1:18);
m1 = zeros(length(lhsRange));
m2 = m1;
k=0;
for vs = 1e9
for ls = 2e-11
    for lhs = lhsRange
        k = k+1;
        M = [[2*lh*v^2 2*lhs*v*vs];[2*lhs*v*vs 2*ls*vs^2]];
        eigenvalues = eig(M);
        m1(k) = eigenvalues(1); m2(k) = eigenvalues(2);
%         if(m2 < m1)
%            fprintf('Crossing with \lambda_HS = %.2g, \lambda_S = %.2g\n',lhs,ls);
%            break; 
%         end
    end
end
end
loglog(lhsRange, sqrt(m1),'-r','LineWidth',1); hold on;
loglog(lhsRange, sqrt(m2),'-b','LineWidth',1);
xlabel('log_{10}\lambda_{H\sigma}','FontSize',20);
ylabel('log_{10}m/GeV','FontSize',20);
title(['\lambda_H = ', num2str(lh,5), ', \lambda_S = ', num2str(ls,5),', v_\sigma = ', num2str(vs,5), 'GeV'],'FontSize',20);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'TickLength',[0.025 0.025]);
set(gca,'FontSize',15);
legend('{\fontsize{15}m_H}','{\fontsize{15}m_\sigma}','Location','NorthWest');