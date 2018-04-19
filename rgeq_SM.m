% SMASH-RGE project
% (C) Timo Karkkainen 2017-2018
% An incredibly large amount of typos removed and missing terms included
% 14.3.2018 >_<
function z = rgeq_SM( t, x )
% CONSTANTS
q = -1/3;           % Extra quark-like field charge in e-units
nf = 6;             % Number of flavors
ng = nf/2;          % Number of generations
% Phys. Rev. D 46.3945
b1 = -4/3*ng-1/10;
b2 = 22/3-4/3*ng-1/6;
b3 = 11-4/3*ng;
b = diag([0 136/3 102]) - ng*[[19/5 1/5 11/30];[3/5 49/3 3/2];[44/15 4 76/3]] - [[9/50 3/10 0];[9/10 13/6 0];[0 0 0]];
C = [[17/10 1/2 3/2];[3/2 3/2 1/2];[2 2 0]];
%         1   2   3   4   5     6   7        8      9  10       11        12     13  14
% x0 = [g10 g20 g30 yu0 yd0   ye0 yf0 lambdaH0 muH0^2 yq0 lambdaS0 lambdaHS0 muS0^2 yn0];

% FULL SET HERE, note: mh0 = \mu_H^2 and ms0 = \mu_S^2, lhs0 = \lambda_{HS}
g10 = x(1); g20 = x(2); g30 = x(3);  yu0 = x(4);  yd0 = x(5);   ye0 = x(6);
lh0 = x(7); mh0 = x(8);

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
    - (43/80*g10^2 -9/16*g20^2 + 16*g30^2)*yd0^2 + 5/2*y4 + (9/200 + 29/45*ng)*g10^4 ...
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

z = [g1SM g2SM g3SM yuSM ydSM yeSM lhSM mhSM]';

end

