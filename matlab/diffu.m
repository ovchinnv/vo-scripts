% 6/16/2015 : starting from windowed time series, compute diffusion coefficient assuming
% a linear PMF in each window.
% formally correct for short times only, and therefore does not make good use of simulation statistics
%clear;

addpath('../');

%qnoplot=1 ; 

potential;

flag='k100_a0';

rhist ; % load coordinate data
if ~exist('qplanar') 
 if (qplanar==1)
  error('qplanar = 1, but qwindow = 1 required');
 end
end
%
return
% string coords
pstring ;
if ~exist('xstr')
 xstr=xi ; % initial coordinates defined in potential
 ystr=yi ;
end
%
% subsample series if desired
nstep=1 ;
nsamp=floor((niter-1)/nstep)+1;
%
% compute displacement coordinate (delta)
fac=[ 0.5 ones(1,np-2) 0.5 ] ;
% string tangent vector
dx=[xstr(2)-xstr(1) xstr(3:np)-xstr(1:np-2) xstr(np)-xstr(np-1) ];
dy=[ystr(2)-ystr(1) ystr(3:np)-ystr(1:np-2) ystr(np)-ystr(np-1) ];
ds2=1./(dx.^2 + dy.^2); ds=sqrt(ds2);

% vectorize calculations for speed
% target string nodes
xstrv=ones(nsamp,1)*[xstr(1) xstr(1:np-1)] ;
ystrv=ones(nsamp,1)*[ystr(1) ystr(1:np-1)] ;
%
dxv=ones(nsamp,1)*dx;
dyv=ones(nsamp,1)*dy;
ds2v=ones(nsamp,1)*ds2 ;
dsv=ones(nsamp,1)*ds ;
xstrv=ones(nsamp,1)*[xstr(1) xstr(1:np-1)] ;
ystrv=ones(nsamp,1)*[ystr(1) ystr(1:np-1)] ;
delta0v=ones(nsamp,1)*delta0 ;
eps0v=ones(nsamp,1)*eps0 ;

dxx = xhist(1 : nstep : niter,:) - xstrv ; %x-displacement to string node
dyy = yhist(1 : nstep : niter,:) - ystrv ; %y-displacement to string node
delta = (dxx.*dxv + dyy.*dyv).*ds2v ;
%
% loop over each window and compute drift and diffusion coefficients
%
clear Diff Dgamma;

n=length(delta) ;
maxrun=n; % maximum run length
maxrun=500;
%
dall=zeros(n,maxrun) - 9999 ; % initialize to a negative value
d2all=zeros(n,maxrun) - 9999 ; % initialize to a negative value

%
for ind=1:np
 disp(ind);
 dd=delta(:,ind) - delta0(ind); % displacement from restraint center
 dbox = dwidth(ind) ; % diffusion box size (half-window to which RW is restrained)
%
 dini = dbox/4 ; % size of window from which to sample trajectory starting points
 inds = find ( abs(dd(1:end-1))  < dini ) ; % will only consider these points for calculations
 % loop over all indices and record displacements
 % create displacement matrix
 clear dinst;
%
 for i=1:length(inds)
 % inst. trajectory
  ibeg = inds(i) ;
  nmax=n-ibeg; % remaining samples in time series
  iend=min(nmax,maxrun) ;

  d  = ( dd(ibeg+1:ibeg+iend)-dd(ibeg) )     ;
  dall(i,1:iend)=d;

  d2all(i,1:iend)=d.^2;

 end % over all runs i
%
 nsample = sum(d2all > -1) ; % number of entries ne 0 ;
%
 clear rho err dmean ;
 for j=1:length(nsample) %  compute statistics for fitting (over all times)
  k=find(d2all(:,j)>-1) ;
% compare probability distribution with a Gaussian using KSL test
% if data not reasonably Gaussian, stop fitting at this j
  ddist=dall(k,j)/ds(ind) ; % convert back to absolute distance units
  dmean(j)=mean(ddist) ;
  sdist=std(ddist) ;
%
%% for error computation :
% compute approximate number of samples using correlation function :
%%%%%%
  alpha = 3.;
  dn=(ddist-dmean(j))/sdist;
  dhat=fft(dn);
  dhat2=dhat.*conj(dhat);
  dc=ifft(dhat2);
  dc=dc(1:floor(length(dc)/2)) ;
  dc=dc/dc(1) ; % should be same as below
% find the point at which dc crosses zero; set block length to alpha times that interval
  k=find(dc<0,1); k=floor(k*alpha) ; iblock=min(k, floor(length(dn)/2));
%
  rho(j)=sdist.^2; % var computed about mean

  nsamp=length(ddist)/iblock ; %effective number of samples
  rhoerr(j)= ( var(ddist.^2) + (rho(j)^2)/nsamp )/nsamp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% err(j)=var(d2s(k,j)) ; % unclear how this should be computed
% err(j)=var(d2s(k,j)) / (nsample(j)) ;
% apply Kolmogorov-Smirnoff-Lilliefors to determine whether sample is normally-distributed
% [h, p(j)]=kstest( (ddist-mean(ddist))/std(ddist) )
%  if ( kstest( (ddist-dmean(j))/sdist) == 1)
%  if ( lillietest( (ddist-mean(ddist))/std(ddist), 0.5) == 1)
  if (j>0)
   jmax=j;
   break;
  end % KSL test
 end % over all times

 range=1:jmax;

 tstep = nstep * dtau * qsave ;
 t = tstep * [1:maxrun] ;

% NOTE : we are working in absolute distance units -- see multiplication of series by dx(i) above
% Dgamma(ind) = - (dmean(range)*t(range)')/(t(range)*t(range)')
% Diff(ind)=(rho(range)*t(range)')/(t(range)*t(range)')/2 
% using a single index
 mind=1;

 Dgamma(ind) = - dmean(mind)/(mind*tstep) 

 Diff(ind)   = 0.5 * rho(mind)/(mind*tstep) 
 Diffe(ind)  = 0.25 / (mind*tstep)^2 * rhoerr(mind) 
% NOTE: may need to correct Dgamma by adding D'

 dfdx(ind)= Dgamma(ind) / Diff(ind) * kBT ; % this gives alpha -- dF/dx but in internal units

end % all point indices

Diff
Dgamma

% plot FE
dfdxc = [ 0.5 * ( dfdx(2:end) + dfdx(1:end-1 ) ) ] ;

dstr= sqrt ( diff(xstr).^2 + diff(ystr).^2 ) ;
fe=[0 cumsum(dfdxc.*dstr) ];

% correct for nonuniform diffusion
fecorr= [0 cumsum( 2*kBT*diff(Diff)./(Diff(1:end-1) + Diff(2:end)) ) ] ; % differentiate wrt. s, but then integrate wrt s, so omit s altogether

figure('position',[300 300 600 200]) ;  hold on ; box on;
plot(fe+fecorr, 'k.') ;
plot(fe, 'k:') ; % uncorrected
% show potential energy curve
%return
%v=inline(vstr);
plot(v(xstr,ystr)-v(xstr(1),ystr(1)),'r');
xlabel('\it x');
ylabel('\it \beta U(\phi)');
set(gca, 'fontsize',14);
set(gcf,'paperpositionmode','auto');
print(gcf,'-depsc2','diffu_fe.eps');

figure('position',[300 300 600 200]) ;  hold on ; box on;
%plot([0:np-1]/(np-1), Diff, 'k.-') ;
errorbar(Diff, sqrt(Diffe),'k.-') ;
% show exact curve
if (0)
 D=inline(Dstr);
 plot(xx, D(xx),'k');
 axis([-pi pi 0 .35]);
end
xlabel('\it x');
ylabel('\it D(x)');
set(gca, 'fontsize',14);
set(gcf,'paperpositionmode','auto');
print(gcf,'-depsc2','diffu_diff.eps');

% save fe profile and diffusion coeff. profile
ftag='diffu' ;
fname=[ftag,'_fe',num2str(np),'.mat'];
save(fname, 'fe','fecorr','Diff','Diffe', 'np');
