% SMASH-RGE project
% (C) Timo Karkkainen 2017-2018
% An incredibly large amount of typos removed and missing terms included
% 14.3.2018 >_<
% Repaired heavy Majorana neutrinos: now input is Y(1,1) and it creates
% a hierarchical Yukawa matrix Y = [Y(1,1) 3*Y(1,1) 3*Y(1,1)]
function z = rgeq( t, x )
% CONSTANTS
q = -1/3;           % Extra quark-like field charge in e-units -1/3 (or 2/3)
nf = 6;          % Number of flavors
ng = nf/2;       % Number of generations
BSM = true;      % Turn on/off BSM effects
% Phys. Rev. D 46.3945
b1 = -4/3*ng-1/10;
b2 = 22/3-4/3*ng-1/6;
b3 = 11-4/3*ng;
b = diag([0 136/3 102]) - ng*[[19/5 1/5 11/30];[3/5 49/3 3/2];[44/15 4 76/3]] - [[9/50 3/10 0];[9/10 13/6 0];[0 0 0]];
C = [[17/10 1/2 3/2];[3/2 3/2 1/2];[2 2 0]];
%         1   2   3   4   5     6   7        8      9  10       11        12     13  14
% x0 = [g10 g20 g30 yu0 yd0   ye0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];

% FULL SET HERE, note: mh0 = \mu_H^2 and ms0 = \mu_S^2, lhs0 = \lambda_{HS}
g10 = x(1); g20 = x(2); g30 = x(3);  yu0 = x(4);  yd0 = x(5);   ye0 = x(6);  yf0 = x(7);
lh0 = x(8); mh0 = x(9); yq0 = x(10); ls0 = x(11); lhs0 = x(12); ms0 = x(13); yn0 = x(14);

% NO SCALAR MASSES HERE
%g10 = x(1); g20 = x(2); g30 = x(3);  yu0 = x(4);   yd0 = x(5);   ye0 = x(6);  yf0 = x(7);
%lh0 = x(8); yq0 = x(9); ls0 = x(10); lhs0 = x(11); yn0 = x(12);

% The notations from Castano paper for SM RGE's
y2 = 3*yu0^2 + 3*yd0^2 + ye0^2;
y4 = (17/20*g10^2 + 9/4*g20^2 + 8*g30^2)*yu0^2 ...
    + (1/4*g10^2 + 9/4*g20^2 + 8*g30^2)*yd0^2 ...
    + 3/4*(g10^2 + g20^2)*ye0^2;
chi4 = 9/4*(3*yu0^4 + 3*yd0^4 + ye0^4 - 2/3*yu0^2*yd0^2);
betau1 = 3/2*(yu0^2 - yd0^2) + y2 - (17/20*g10^2 + 9/4*g20^2 + 8*g30^2);
betad1 = 3/2*(yd0^2 - yu0^2) + y2 - (1/4*g10^2 + 9/4*g20^2 + 8*g30^2);
betae1 = 3/2*ye0^2 + y2 - 9/4*(g10^2 + g20^2);
betau2 = 3/2*yu0^4 - yu0^2*yd0^2 - 1/4*yd0^2*yu0^2 + 11/4*yd0^4 ...
    + y2*(5/4*yd0^2 - 9/4*yu0^2) - chi4 + 3/2*lh0^2 - 2*lh0*(3*yu0^2 + yd0^2) ...
    + (223/80*g10^2 + 135/16*g20^2 + 16*g30^2)*yu0^2 ...
    - (43/80*g10^2 - 9/16*g20^2 + 16*g30^2)*yd0^2 + 5/2*y4 + (9/200 + 29/45*ng)*g10^4 ...
    - 9/20*g10^2*g20^2 + 19/15*g10^2*g30^2 - (35/4 - ng)*g20^4 + 9*g20^2*g30^2 ...
    - (404/3 - 80/9*ng)*g30^4;
betad2 = 3/2*yd0^4 - yd0^2*yu0^2 - 1/4*yu0^2*yd0^2 + 11/4*yu0^4 ...
    + y2*(5/4*yu0^2 - 9/4*yd0^2) - chi4 + 3/2*lh0^2 - 2*lh0*(3*yd0^2 + yu0^2) ...
    + (187/80*g10^2 + 135/16*g20^2 + 16*g30^2)*yd0^2 ...
    - (79/80*g10^2 - 9/16*g20^2 + 16*g30^2)*yu0^2 + 5/2*y4 - (29/200 + 1/45*ng)*g10^4 ...
    - 27/20*g10^2*g20^2 + 31/15*g10^2*g30^2 - (35/4 - ng)*g20^4 + 9*g20^2*g30^2 ...
    - (404/3 - 80/9*ng)*g30^4;
betae2 = 3/2*ye0^4 - 9/4*y2*ye0^2 - chi4 + 3/2*lh0^2 - 6*lh0*ye0^2 ...
    + (387/80*g10^2 + 135/15*g20^2)*ye0^2 + 5/2*y4 + (51/200 + 11/5*ng)*g10^4 ...
    + 27/20*g10^2*g20^2 - (35/4 - ng)*g20^4;
h = 3*yu0^4 + 3*yd0^4 + ye0^4;
betal1 = 12*lh0^2 - (9/5*g10^2 + 9*g20^2)*lh0 ...
    + 9/4*(3/25*g10^4 + 2/5*g10^2*g20^2 + g20^4) + 4*y2*lh0 - 4*h;
betal2 = -78*lh0^3 + 18*(3/5*g10^2 + 3*g20^2)*lh0^2 ...
    - ((313/8 - 10*ng)*g20^4 - 117/20*g10^2*g20^2 + 9/25*(229/4 + 50/9*ng)*g10^4)*lh0 ...
    + (497/8 - 8*ng)*g20^6 - 3/5*(97/24 + 8/3*ng)*g10^2*g20^4 ...
    - 9/25*(239/24 + 40/9*ng)*g10^4*g20^2 - 27/125*(59/24 + 40/9*ng)*g10^6 ...
    - 64*g30^2*(yu0^4 + yd0^4) - 8/5*g10^2*(2*yu0^4 - yd0^4 + 3*ye0^4) ...
    - 3/2*g20^4*y4 + 10*lh0*((17/20*g10^2 + 9/4*g20^2 + 8*g30^2)*yu0^2  ...
    + (1/4*g10^2 + 9/4*g20^2 + 8*g30^2)*yd0^2 + 3/4*(g10^2 + g20^2)*ye0^2) ...
    + 3/5*g10^2*((-57/10*g10^2 + 21*g20^2)*yu0^2 + (3/2*g10^2 + 9*g20^2)*yd0^2 ...
    + (-15/2*g10^2 + 11*g20^2)*ye0^2) - 24*lh0^2*y2 - lh0*h + 6*lh0*yu0^2*yd0^2 ...
    + 20*(3*yu0^6 + 3*yd0^6 + ye0^6) - 12*(yu0^2*(yu0^2 + yd0^2)*yd0^2);

% From 1303.4364, scaled! \hat{\lambda} = \lambda/2, factor of 2 difference
betah1 = -9/20*(5/3)*g10^2 - 9/4*g20^2 + 3*yu0^2 + 3*yd0^2 + ye0^2 + 6*lh0/2;
betah2 = (5/3)^2*g10^4*(ng/2 + 471/800) + g20^4*(5/2*ng - 385/32) + 36/5*(5/3)*g10^2*lh0/2 + 36*g20^2*lh0/2 ...
    - 30*lh0^2/4 - 36*yu0^2*lh0/2 - 36*yd0^2*lh0/2 - 12*ye0^2*lh0/2 + 9/16*(5/3)*g10^2*g20^2 + 17/8*(5/3)*g10^2*yu0^2 ...
    + 5/8*(5/3)*g10^2*yd0^2 + 15/8*(5/3)*g10^2*ye0^2 + 45/8*g20^2*yu0^2 - 9/4*ye0^4 + 45/8*g20^2*yd0^2 ...
    + 15/8*g20^2*ye0^2 + 20*g30^2*yu0^2 + 20*g30^2*yd0^2 - 27/4*yd0^4 - 21/2*yd0^2*yu0^2 ...
    - 27/4*yu0^4;

% SM contributions, from Castano paper
g1SM = -b1*g10^2/(16*pi^2) ...
    - b(1,1)*g10^2*g10^3/(16*pi^2)^2 ...
    - b(2,1)*g20^2*g10^3/(16*pi^2)^2 ...
    - b(3,1)*g30^2*g10^3/(16*pi^2)^2 ...
    - g10^3/(16*pi^2)^2*(C(1,1)*yu0^2 + C(1,2)*yd0^2 + C(1,3)*ye0^2);
g2SM = -b2*g20^2/(16*pi^2) ...
    - b(1,2)*g10^2*g20^3/(16*pi^2)^2 ...
    - b(2,2)*g20^2*g20^3/(16*pi^2)^2 ...
    - b(3,2)*g30^2*g20^3/(16*pi^2)^2 ...
    - g20^3/(16*pi^2)^2*(C(2,1)*yu0^2 + C(2,2)*yd0^2 + C(2,3)*ye0^2);
g3SM = -b3*g30^2/(16*pi^2) ...
    - b(1,3)*g10^2*g30^3/(16*pi^2)^2 ...
    - b(2,3)*g20^2*g30^3/(16*pi^2)^2 ...
    - b(3,3)*g30^2*g30^3/(16*pi^2)^2 ...
    - g30^3/(16*pi^2)^2*(C(3,1)*yu0^2 + C(3,2)*yd0^2 + C(3,3)*ye0^2);
yuSM = yu0*(betau1/(16*pi^2) + betau2/(16*pi^2)^2);
ydSM = yd0*(betad1/(16*pi^2) + betad2/(16*pi^2)^2);
yeSM = ye0*(betae1/(16*pi^2) + betae2/(16*pi^2)^2);
lhSM = betal1/(16*pi^2) + betal2/(16*pi^2)^2;
mhSM = 2*mh0*(betah1/(16*pi^2) + betah2/(16*pi^2)^2);

% BSM contributions added to them here, from SMASH paper
ynMatrix = diag([yn0 3*yn0 3*yn0]); % Heavy neutrino matrix hierarchy created here
g1 = g1SM + 12/(80*pi^2)*g10^3*q^2 + 1/(16*pi^2)^2*(108/25*g10^5*q^4 ...
    + q^2*(48/5*g10^3*g30^2 - 18/5*g10^3*yq0^2) - 3*g10^3/10*yf0^2);
g2 = g2SM - g20^3/(2*(16*pi^2)^2)*yf0^2;
g3 = g3SM + g30^3/(24*pi^2) ...
    + 1/(16*pi^2)^2*(38*g30^5/3 - g30^3*yq0^2 + 6/5*g10^2*g30^3*q^2);
%g3 = g3SM + 1/(16*pi^2)^2*( - g30^3*yq0^2 + 6/5*g10^2*g30^3*q^2);
% Top Yukawa
yu = yuSM + yu0/(16*pi^2)*yf0^2 + 1/(16*pi^2)^2*(40/9*g30^4*yu0 ...
    + 2*lhs0^2*yu0 + 29/25*g10^4*q^2*yu0 + 1/8*yf0^2*(10*yd0^2*yu0 ...
    + 3*g10^2*yu0 + 15*g20^2*yu0 - 18*yu0^3) + ye0^2*yf0^2*yu0/2 ...
    - 9*yf0^4*yu0/4 - 3*trace(ynMatrix^2)*yf0^2*yu0/4);
% Bottom Yukawa
yd = ydSM + yd0/(16*pi^2)*yf0^2 + 1/(16*pi^2)^2*(40/9*g30^4*yd0 ...
    + 2*lhs0^2*yd0 - 1/25*g10^4*q^2*yd0 + 1/8*yf0^2*(3*g10^2*yd0 ...
    + 15*g20^2*yd0 + 10*yd0*yu0^2 - 18*yd0^3) + ye0^2*yf0^2*yd0/2 ...
    - 9*yf0^4*yd0/4 - 3*trace(ynMatrix^2)*yf0^2*yd0/4);
% Tau Yukawa
ye = yeSM + 1/(16*pi^2)*(ye0*yf0^2 - 3*yf0^2*ye0/2) ...
    + 1/(16*pi^2)^2*(15*yf0^2*ye0*yd0^2/4 + 2*ye0*lhs0^2 ...
    + g10^2*(3*ye0*yf0^2/8 - 27*yf0^2*ye0/16) ...
    + g20^2*(15*ye0*yf0^2/8 + 9*yf0^2*ye0/16) ...
    - 9*ye0*yf0^4/4 + ye0*yf0^2*ye0^2/2 + 15*yf0^2*ye0*yu0^2/4 ...
    + 5*yf0^2*ye0*yf0^2/4 + 5*yf0^2*ye0*ye0^2/4 ...
    + 11*yf0^4*ye0/4 - yf0^2*ye0^3 - 9*ye0^3*yf0^2/4 - ye0^2*yf0^2*ye0/4 ...
    + 7*yf0*trace(ynMatrix^2)*yf0*ye0/8 - 3*ye0*trace(ynMatrix^2)*yf0^2/4 + 99/25*ye0*g10^4*q^2);
% SM Higgs self-coupling
lh = lhSM + 1/(16*pi^2)*(4*lhs0^2 + 4*yf0^2*lh0 - 2*yf0^4) ...
    + 1/(16*pi^2)^2*(18/125*g10^4*q^2*(25*lh0 - 6*g10^2 - 10*g20^2) ...
    - 40*lh0*lhs0^2 - 32*lhs0^3 - 24*yq0^2*lhs0^2 + g10^2*(3*yf0^2*lh0/2 ...
    - 3*g20^2*yf0^2/10) + 15/2*g20^2*yf0^2*lh0 - 9*g10^4/100*yf0^2 ...
    - 3*g20^4*yf0^2/4 - 14*ye0^2*yf0^2*lh0 - 48*yf0^2*lh0^2 - yf0^4*lh0 ...
    - 2*yf0^2*ye0^4 - 2*yf0^4*ye0^2 + 10*yf0^6 - 3*trace(ynMatrix^2)*yf0^2*lh0 ...
    - 4*trace(ynMatrix^2)*lhs0^2 + 2*yn0*yf0^2*yn0*yf0^2 + 2*trace(ynMatrix^2)*yf0^4);
% SM Higgs quadratic coupling
mhBSM = 1/(16*pi^2)*(4*lhs0*ms0 + 2*yf0^2*mh0) ...
   + 1/(16*pi^2)^2*(9/5*g10^4*q^2*mh0 - ms0*(16*lhs0^2 + 24*yq0^2*lhs0) ...
   - 4*mh0*lhs0^2 + mh0*(3*g10^2*yf0^2/4 + 15*g20^2*yf0^2/4 - 24*yf0^2*lh0 ...
   - 7*ye0^2*yf0^2 - 9*yf0^4/2) - 3/2*trace(ynMatrix^2)*yf0^2*mh0 - 4*trace(ynMatrix^2)*lhs0*ms0);
mh = mhSM + mhBSM;
% BSM parameters
% Dirac neutrino mass Yukawa coupling
yf = 1/(16*pi^2)*(3*yf0*yd0^2 - 9/20*yf0*g10^2 - 9*yf0*g20^2/4 + 3*yf0*yu0^2 ...
    + yf0*yf0^2 + yf0*ye0^2 + 3*yf0^3/2 - 3*ye0^2*yf0/2 + yf0*trace(ynMatrix^2)/2) ...
    + 1/(16*pi^2)^2*(g10^2*(5*yf0*yd0^2/8 - 27/20*yf0*g20^2 + 17*yf0*yu0^2/8 ...
    + 3*yf0*yf0^2/8 + 15*yf0*ye0^2/8 + 279*yf0^3/80 - 243*ye0^2*yf0/80) ...
    + g20^2*(45*yf0*yd0^2/8 + 45*yf0*yu0^2/8 + 15*yf0*yf0^2/8 + 15*yf0*ye0^2/8 ...
    + 135*yf0^3/16 + 9*ye0^2*yf0/16) + g30^2*(20*yf0*yd0^2 + 20*yf0*yu0^2) ...
    + yu0^2*(3*yf0*yd0^2/2 - 27*yf0^3/4 + 15*ye0^2*yf0/4) - 27*yf0*yd0^4/4 ...
    + yd0^2*(15*ye0^2*yf0/4 - 27*yf0^3/4) + 21*yf0*g10^4/40 - 23*yf0*g20^4/4 ...
    + 6*yf0*lh0^2 + 2*yf0*lhs0^2 - 27*yf0*yu0^4/4 - 9*yf0*yf0^4/4 + yf0*yf0^2*ye0^2/2 ...
    - 9*yf0*ye0^4/4 - 12*yf0^3*lh0 - 9*yf0^3*yf0^2/4 - 9*yf0^3*ye0^2/4 + 3*yf0^5/2 ...
    - yf0^2*ye0^2*yf0/4 + 5*ye0^2*yf0*yf0^2/4 + 5*ye0^2*yf0*ye0^2/4 - ye0^2*yf0^3 ...
    + 11*ye0^4*yf0/4 - 3*yf0*yf0^2*trace(ynMatrix^2)/4 - 4*yf0*trace(ynMatrix^2)*lhs0 + 7*yf0*yn0*yf0*yn0*yf0/4 ...
    - yf0*trace(ynMatrix^2)*yf0^2/8 - yf0*trace(ynMatrix^4)/8 - 9/4*yf0*trace(ynMatrix^2)*yq0^2 - 3*yf0*trace(ynMatrix^2)*trace(ynMatrix^2)/8 ...
    + 9/25*yf0*g10^4*q^2);
% Singlet scalar self-coupling
ls = 1/(16*pi^2)*(8*lhs0^2 + 20*ls0^2 - 6*yq0^4 + 12*yq0^2*ls0 + 2*trace(ynMatrix^2)*ls0 - trace(ynMatrix^4)) ...
    + 1/(16*pi^2)^2*(lhs0^2*(48*g20^2 + 48*g10^2/5 - 48*yd0^2 - 80*ls0 - 48*yu0^2) ...
    + yq0^4*(6*ls0 - 32*g30^2) + yq0^2*(80*g30^2*ls0 - 120*ls0^2) - 64*lhs0^3 ...
    - 240*ls0^3 + 24*yq0^6 + 36/5*g10^2*q^2*yq0^2*(5*ls0 - 2*yq0^2) - 16*ye0^2*lhs0^2 ...
    - 16*yf0^2*lhs0^2 + 4*yf0^2*trace(ynMatrix^4) + ls0*(trace(ynMatrix^4) - 6*trace(ynMatrix^2)*yf0^2) + 4*trace(ynMatrix^6) - 20*trace(ynMatrix^2)*ls0^2);
% Singlet-doublet coupling
lhs = 1/(16*pi^2)*(lhs0*(6*yd0^2 - 9*g20^2/2 - 9*g10^2/10 + 12*lh0 + 8*ls0 + 6*yu0^2 ...
    + 6*yq0^2) + 8*lhs0^2 + lhs0*(2*ye0^2 + 2*yf0^2) + trace(ynMatrix^2)*lhs0 - 2*trace(ynMatrix^2)*yf0^2) ...
    + 1/(16*pi^2)^2*(lh0*(lhs0*(72*g10^2/5 - 72*yd0^2 + 72*g20^2 - 72*yu0^2) ...
    - 144*lhs0^2) + lhs0^2*(6*g10^2/5 - 24*yd0^2 + 6*g20^2 - 96*ls0 - 24*yu0^2 ...
    - 24*yq0^2) + lhs0*(g30^2*(40*yd0^2 + 40*yu0^2 + 40*yq0^2) + g10^2*(5*yd0^2/4 ...
    + 9*g20^2/8 + 17*yu0^2/4) + g20^2*(45*yd0^2/4 + 45*yu0^2/4) - 21*yd0^2*yu0^2 ...
    - 27*yd0^4/2 + 1671*g10^4/400 - 145*g20^4/16 - 40*ls0^2 - 27*yu0^4/2 - 9*yq0^4 ...
    - 48*yq0^2*ls0) - 60*lh0^2*lhs0 - 44*lhs0^3 - 9/25*g10^2*q^2*(g10^2*(18*yq0^2 ...
    - 5*lhs0) - 50*yq0^2*lhs0) + lhs0*(g10^2*(15*ye0^2/4 + 3*yf0^2/4) ...
    + g20^2*(15*ye0^2/4 + 15*yf0^2/4) - 9*ye0^4/2 - 7*ye0^2*yf0^2 - 9*yf0^4/2 ...
    - lh0*(24*ye0^2 + 24*yf0^2) + 7*trace(ynMatrix^2)*yf0^2/2 - 3*trace(ynMatrix^4)/2 - 8*trace(ynMatrix^2)*ls0) ...
    - lhs0^2*(4*trace(ynMatrix^2) + 8*ye0^2 + 8*yf0^2) + 5*yf0^2*trace(ynMatrix^4) + 4*yn0*yf0^2*yn0*yf0^2 ...
    + 7*trace(ynMatrix^2)*yf0^4 - trace(ynMatrix^2)*yf0*ye0^2*yf0);
% Neutrino Majorana mass Yukawa coupling
yn = 1/(16*pi^2)*(yf0^2*yn0 + yn0*yf0^2 + trace(ynMatrix^3) + 3*yn0*yq0^2 + trace(ynMatrix^3)/2) ...
    + 1/(16*pi^2)^2*(g10^2*(51*yf0^2*yn0/40 + 51*yn0*yf0^2/40) - yd0^2*(9*yf0^2*yn0/2 ...
    + 9*yn0*yf0^2/2) + g20^2*(51*yf0^2*yn0/8 + 51*yn0*yf0^2/8) + 20*g30^2*yn0*yq0^2 ...
    - lhs0*(8*yf0^2*yn0 + 8*yn0*yf0^2) + 4*yn0*lhs0^2 - yf0^4*yn0/4 + 4*yf0^2*yn0*yf0^2 ...
    - (9*yf0^2*yn0/2 + 9*yn0*yf0^2/2)*yu0^2 - 3*yf0^2*yn0*yf0^2/2 - 3*yf0^2*yn0*ye0^2/2 ...
    - yf0*ye0^2*yf0*yn0/4 - yn0*yf0^2*trace(ynMatrix^2)/4 - 3*yn0*yf0^2*yf0^2/2 - 3*yn0*yf0^2*ye0^2/2 ...
    - yn0*yf0^4/4 - yn0*yf0*ye0^2*yf0/4 - trace(ynMatrix^2)*yf0^2*yn0/4 + 7*trace(ynMatrix^5)/4 ...
    - 8*trace(ynMatrix^3)*ls0 - 9*trace(ynMatrix^3)*yq0^2/2 - 3*trace(ynMatrix^5)/4 + 4*yn0*ls0^2 - 9/2*yn0*yq0^4 ...
    - 3*yn0*yf0^2*trace(ynMatrix^2)/2 - 3*yn0*trace(ynMatrix^4)/4 + 9*g10^2*q^2*yn0*yq0^2);
% New quark Yukawa coupling
yq = 1/(16*pi^2)*(4*yq0^3 - 8*g30^2*yq0 - 18/5*g10^2*q^2*yq0 + 1/2*trace(ynMatrix^2)*yq0) ...
    + 1/(16*pi^2)^2*(yq0*(4*lhs0^2 + 4*ls0^2 - 932/9*g30^4) + yq0^3*(92*g30^2/3 - 8*ls0) ...
    - 29*yq0^5/4 + 3/50*g10^2*q^2*yq0*(g10^2*(102*q^2 + 211) - 80*g30^2 + 230*yq0^2) ...
    - yq0*(3*trace(ynMatrix^4)/4 + 3*trace(ynMatrix^2)*yf0^2/2) - 3*trace(ynMatrix^2)*yq0^3/4);
% Quadratic singlet scalar coupling
ms = 1/(16*pi^2)*(8*mh0*lhs0 + ms0*(8*ls0 + 6*yq0^2 + trace(ynMatrix^2))) ...
   + 1/(16*pi^2)^2*(mh0*(lhs0*(48*g10^2/5 - 48*yd0^2 + 48*g20^2 - 48*yu0^2) ...
   - 32*lhs0^2) + ms0*(40*g30^2*yq0^2 - 8*lhs0^2 - 40*ls0^2 - 9*yq0^4 ...
   - 48*yq0^2*ls0) + 18*g10^2*q^2*yq0^2*ms0 - mh0*lhs0*(16*ye0^2 ...
   + 16*yf0^2) - ms0*(3/2*trace(ynMatrix^4) + 3*trace(ynMatrix^2)*yf0^2 + 8*trace(ynMatrix^2)*ls0));
if (~BSM)
   g1 = g1SM;
   g2 = g2SM;
   g3 = g3SM;
   yu = yuSM;
   yd = ydSM;
   ye = yeSM;
   lh = lhSM;
   mh = mhSM;
   yf = 0;
   yq = 0;
   ls = 0;
   lhs = 0;
   ms = 0;
   yn = 0;
end


% All except scalar masses
%z = [g1 g2 g3 yu yd ye yf lh yq ls lhs yn]';
% All parameters
z = [g1 g2 g3 yu yd ye yf lh mh yq ls lhs ms yn]';

end

