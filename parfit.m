% FILE  parfit.m  begins.
% USES: fmodel, fobj
% MME99
global xdata;
global ydata;
xdata= 0:0.05:1.1;
lam1=[0.2 1.5 0.7];
y=fmodel(lam1,xdata); % Generate syntetic data 
                      % with parameters lam1
ydata= y.*(0.97+0.15*rand(1,length(xdata)));
                      % Add some "errors"
%fprintf('%8.4f',xdata); fprintf('\n');
%fprintf('%8.4f',ydata);
lam0=[-1.8 1 10];         % Initial guess for lambda
y0=fobj(lam0);        % Initial value of object function
% OLD MATLAB:lam=fmins('fobj',lam0);
lam=fminsearch('fobj',lam0);
                      % lam is the fitted value for
                      % the parameter vector
x=0:0.05:1.5;
yfit=fmodel(lam, x);
yfinal=fobj(lam);     % Final value of the object function
clf;
axes('FontSize',[15],'FontWeight','bold'); hold on;
title(['Object function values: start =' num2str(y0) ', final = ' ...
       num2str(yfinal)])
plt=plot(x,yfit,xdata,ydata,'k.','MarkerSize',30);  grid;
txt1=' {\bf Fitted curve (solid)}';
text(x(5),yfit(5),txt1,'FontWeight','bold','FontSize',[16]);
txt2=' {\bf Data point (dots)}';
text(xdata(8),ydata(8),txt2,'FontWeight','bold','FontSize',[16]);
ylabel('ydata '); xlabel('xdata (parfit/MME99)')
set(plt,'LineWidth',2.5);
% FILE  parfit.m  ends.
