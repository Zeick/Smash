function y=fobj(lambda)
global xdata;
global ydata;
y=norm(fmodel(lambda,xdata)-ydata);
% FILE  fobj.m  ends.
