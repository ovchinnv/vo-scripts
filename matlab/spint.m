% test case for evaluating an integral using splines
% we are interested in functions that look like power laws (c.f. confinement method)
close all;


a=1;
b=1;
n=1;
syms x;

fun='b + a / (x ^ n)';

%fun='b*exp(- a * x )'; % has a maximum at 0

f=inline(eval(fun),'x');
fs=sym(eval(fun));

x0=0.00002;
x1=10 ;

F=int(fs,x0,x1); F=eval(F)
ezplot(fs,x0,x1); 
%return
%set(gca,'xscale','log','yscale','log');

% compute the integral using various methods
% in each case the curve is sampled with an exponentially distributed series

np=15 ;
x=exp(linspace(log(x0),log(x1),np));
fx=f(x);
plot(x,fx,'k.-'); hold on;

set(gca,'xscale','log','yscale','log');


%%% method 1 : integrate using trapezoid rule in log space (Tyka's method)
clear F1 ; F1(1)=0;
for i=1:np-1
 df(i)=log( x(i+1)/x(i) ) * ( x(i+1)*fx(i+1)-x(i)*fx(i) ) / ( log(x(i+1)*fx(i+1))-log(x(i)*fx(i)) );
 F1(i+1) = F1(i) + df(i);
end

Ferr1=F-F1(end)

%%% method 0 (least informed) : integrate using  trapezoid rule in real space

clear F0 ; F0(1)=0;
for i=1:np-1
 df(i) = (x(i+1)-x(i)) * 0.5 * ( fx(i+1) + fx(i)) ;
 F0(i+1) = F0(i) + df(i);
end

Ferr0=F-F0(end)

%%% method 3 (most informed) : increased order of integration in log space
% I find that with this method the integration is more accurate near zero if the function is bounded
% this is actually the case for the confinement simulations because at zero force constant the standard potential 
% prevents extremely large displacements (i.e. the RMSD is bounded)
% This appears to be in agreement with Tyka's statements

xl=log(x);
fl=log(fx);
pol=spline(xl,fl);

npfine=2000;
xfine=linspace(x0,x1,npfine);
xlfine=log(xfine);
flfine=ppval(pol,xlfine);
ffine=exp(flfine);

plot(xfine,ffine,'r.')

% now integrate the fine curve using trapezoid

clear F2 ; F2(1)=0;
for i=1:npfine-1
 df(i) = (xfine(i+1)-xfine(i)) * 0.5 * ( ffine(i+1) + ffine(i)) ;
 F2(i+1) = F2(i) + df(i);
end

Ferr2=F-F2(end)

% method 4 %% simmilar to above but integration in log space 
%
npfine=2000;
xlfine=linspace(log(x0),log(x1),npfine);
flfine=ppval(pol,xlfine);
xfine=exp(xlfine);
ffine=exp(flfine);

plot(xfine,ffine,'g.')

% now integrate the fine curve in log space

clear F3 ; F3(1)=0;
for i=1:npfine-1
 df(i)=log( xfine(i+1)/xfine(i) ) * ( xfine(i+1)*ffine(i+1)-xfine(i)*ffine(i) ) / ( log(xfine(i+1)*ffine(i+1))-log(xfine(i)*ffine(i)) );
 F3(i+1) = F3(i) + df(i);
end

Ferr3=F-F3(end)



% Will perhaps try this method for the confinement integration

