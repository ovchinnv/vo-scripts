
d=load('rmsd-name CA.dat');

figure(98); clf; hold on; box on;
plot(d,'k'); 


%%%%%%%%% fit RMSD to a saturating curve
% assume the following model :
%
% 1  -  y   = exp ( -a x )
%      ---
%      ymax
%
% find ymax, a using Least squares, NR
%
% make sure the residual is defined correctly :
%
% y=ymax*(1-exp(-ax))
%
%
% variables :
y=d ;
x=(1:length(d))' ;
% initial guess :
ymax = max( y ) ;
z=y/ymax;
a =  -sum ( log ( 1 - z(find(z<1)))) / sum (x(find(z<1)))

% plot initial guess solution :
yini=ymax * (1-exp(-a*x)) ;
plot(x,yini,'r.');

% iteratively refine a ymax in the least squares sense:
nmax=100 ;
damp=0.1; % damping factor to slow down solution evolution
amin= 1e-6 ; % smallest allowed exponent
ymaxmin = mean(y) % smallest allowed maximum

for n=1:nmax
 e=exp(-a*x);
 f=ymax*(1-e);
 R=y-f;
%
 dRdy=e-1;
 dRda=-x.*e*ymax;
 d2Rdy=0;
 d2Rdyda=-x.*e;
 d2Rda=x.^2.*e*ymax;
 
 F=[sum(R.*dRdy) ; sum(R.*dRda)];
%
 GradF=[ sum(dRdy.^2), sum( dRdy.*dRda + R.*d2Rdyda) ; ...
         sum( dRdy.*dRda + R.*d2Rdyda) , sum(dRda.^2 + R.*d2Rda) ];

% solve for new parameters:
 dpar=GradF\F ;
 ymax=max(ymaxmin,ymax-dpar(1)*damp);
 a   =max(amin,a-dpar(2)*damp);

 fprintf([' Iteration : ',num2str(n), '; Residual : ',num2str(norm(R)),'\n']);

end

[ymax,a]

% plot final solution :
yfin=ymax * (1-exp(-a*x)) ;
plot(x,yfin,'g.');

legend({'rmsd','initial fit','final fit'},4);

% convergence is assumed when the model function is within 99 of ymax:
xconv = - log ( 1 - 0.99 ) / a ;
fprintf([' 99%% convergence location :',num2str(round(xconv)),'\n']);
% draw a line there :
plot([xconv,xconv],[0 max(y)],'k--');

set(gcf, 'paperpositionmode','auto');
print(gcf,'-depsc2', 'rmsdfit.eps');

