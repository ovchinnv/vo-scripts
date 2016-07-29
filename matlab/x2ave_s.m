% compute average value from _steady-state_ PDF given by Fokker-Planck equation for Brownian motion
% in 1D with reflective boundary conditions, assuming a quadratic potential (below)
% The FPE is given by
% .
% P = D [(γ + 2xε)P + P']' ; γ = α / (kBT) ; ε = δ / (kBT)  assuming the potential has the form U = αx + δx^2 + C
% the domain is I=(0..1) and the boundary conditions are P'(0,t) = P'(1,t) = 0
% for the initial condition, take δ(x-x0) where x0 is in I; (here, δ is the delta distribution)
%
function [x,dx]=x2ave_s(g,e)
%
% g    : gamma above
% e    : epsilon above
%
tole=1e-5 ; % tolerance at which to switch to taylor expansion of steady-state contribution to average
tolg2oe=50 ; % maximum tolerance of g^2/e (otherwise exponential explodes and erf is essentially constant)
tolg=1e-7 ;

sqpi=1.772453850905516;

if ( (abs(e) > tole) && ( g^2 < tolg2oe * abs(e) ) )
 oe=1./e;
 if (e>0) % error function form :
fprintf(' ********************* Erf \n')
  osqe=1./sqrt(e);
  gh=0.5*g*osqe ;
% partition function
  Z = 0.5 * sqpi * osqe * exp ( (gh)^2 ) * ( erf (gh + osqe*e) - erf(gh) ) ;
% partition function derivatives :
  ooZ=1./Z ;
  a=exp(-(g+e)) ;
  b=0.5*g*oe ;
% -dlogZ/dg = <x>:
  ddg = - b - ooZ * 0.5 * oe * ( a - 1 ) ;
% -dlogZ/de = <x^2>:
  dde = (0.5*g*oe)^2 + 0.5*oe - ooZ * 0.5 * oe * ( a * ( 1 - b) + b ) ;
  if (nargout>1) % output derivatives w.r.t. g and e respectively
  
  end
 else % Dawson integral form
fprintf(' ********************* Dawson \n')
  osqe=1./sqrt(-e);
  gh=0.5*g*osqe;
%
  [D1,D1p] = dawson( - (osqe*e+gh) ) ;
  [D2,D2p] = dawson( gh ) ;
  a=exp(-(g+e)) ;
  b=0.5*g*oe ;
%
  Z=osqe * ( a * D1 + D2 )
  ooZ=1./Z ;
% partition function derivatives :
% -dlogZ/dg = <x>:
  ddg = ooZ * osqe * ( a * D1 + 0.5 * osqe * (a*D1p-D2p)) ;
% -dlogZ/de = <x^2>:
  dde = 0.5*oe + ooZ * ( osqe * a * D1 + 0.5 * oe * ( D1p*a*(b-1)-D2p*b ) ) ;
  if (nargout>1) % output derivatives w.r.t. g and e respectively
  
  end
 end
else % use expansions

% expansion for e --> zero
 if ( abs(g) > tolg )
fprintf(' ********************* Intermediate exp. \n')
  eg=exp(-g);
  g2=g*g;
  g3=g*g2;
  g4=g2*g2;
  g5=g3*g2;
  og=1/g;
% note that we do not actually need the partition function
%  Z=og*(1-eg) + e*og^3*(2-(g2+2*g+2)*eg) + (e/g2)^2*og*(12-(g2^2+4*g2*g+12*g2+24*g+24)*0.5*eg)
% -dlogZ/dg = <x>:
  ddg = (og - eg/(1-eg)) ... % 0th order
  - ( (4+2*g)*eg^2+4-eg*(g3+2*g2+2*g+8) )/(g3*(1-eg)^2)*e ... %1st order
  - ( eg*(eg^2*(40+4*g*(6 + g)) + ((120 + 16*g2 + 6*g3 + 2*g4 + 0.5*g5) - eg*(- 0.5*g5 - 2*g4 + 2*g3 + 20*g2 + 48*g + 120) + 24*g)) - 40)/(g5*(1 - eg)^3) * 0.5*e^2 %2nd
% -dlogZ/de = <x^2>:
  dde = (2*og^2 - eg*(g+2)*og/(1-eg)) ... %0th order
  - ( (20+16*g+4*g2)*eg^2+20-eg*(g2^2+4*g*g2+8*g2+16*g+40))/(g2*(1-eg))^2 * e ; %1st order
  %2nd order
  if (nargout>1) % output derivatives w.r.t. g and e respectively
  
  end
 else % use expansion for g --> 0
fprintf(' ********************* Lowest exp. \n')
% -dlogZ/dg = <x>:
  ddg=0.5 - (g+e)/12 + e * g / 180 ; % 1/12 and 1/180
% -dlogZ/de = <x^2>:
  dde=1/3 - g/12 - e*((8-g)/90);
 end

end

x=[ddg dde];
