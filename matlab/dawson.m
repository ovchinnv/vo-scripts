% Dawson's integral (and optional derivative) using Rybicki's method in Numerical Recipes
% accuracy is about single precision
%
%             2     x     2
%           -x      /   -u
%  D(x) =  e     .  |  e   du
%                   /
%                   0
%
% Victor Ovchinnikov, Harvard, 2016
%
function [d,dp]=dawson(x) ;
%
if (nargin<1)
 fprintf('??? Error using ==> dawson\n');
 fprintf('Not enough input arguments\n\n');
 return
end
% constants
NMAX=13;
h=0.235;
A1=2/3;
A2=0.4;
A3=2/7;
%oosqpi=1./sqrt(pi)
oosqpi=0.564189583547756;
c = exp( -( (2.*[1:NMAX]-1.)*h).^2); % see below :
% note that c's depend on NMAX and on h
%return
%c=[0.946272212751303   0.608337775176690   0.251421365590295   0.066801816880469   0.011410457722565   0.001252986972712   0.000088454257852   0.000004014390612 ...
%   0.000000117124894   0.000000002196886   0.000000000026491   0.000000000000205   0.000000000000001];
%
% series expansion solution (default; accurate only for small x -- see below) :
x2=x.*x;
d = x .* (1. - A1*x2 .* (1.-A2*x2 .* (1.-A3 .* x2)));

% sampling theorem representation
xx=abs(x);
n0=2*round( 0.5 * (xx./h) ) ;
xp=xx-n0*h;
e1=exp(2.*xp.*h);
e2=e1.^2.;
d1=n0+1.;
d2=d1-2.;
s=zeros(size(x));
for i=1:NMAX
 s = s + c(i) * (e1 ./ d1 + 1./(d2.*e1));
 d1=d1+2.;
 d2=d2-2.;
 e1=e1.*e2;
end
dsampling =oosqpi * (sign(x).*exp(-xp.^2).*s);

% choose the more accurate solution for larger x :
inds=find(abs(x)>0.2) ;
d(inds)=dsampling(inds) ;

% compute derivative if requested
if (nargout>1)
 dp=1.-2*x.*d; % derivative w.r.t. x
end
