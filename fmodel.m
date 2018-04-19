function y=fmodel(lam,x)
%y= exp(-lam(1)*x)+lam(2)*exp(-lam(3)*x);
y = lam(1)+lam(2)*exp(lam(3)*x);
% FILE  fmodel.m  ends.
