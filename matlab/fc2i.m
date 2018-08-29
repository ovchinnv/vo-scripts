% examine forces from double half-harmonic window simulations
%
close all;

if ~exist('styles')
 styles={'r-*','g-x','b-o','m-','c-','k-','r--','g--','b--','m--','c--','k--', 'r-*','g-s','b-v'};
end
lw=1;
leg={};

addpath '~/scripts/matlab'; %for xave_

kboltz=1.987191d-3;        % boltzmann constant
Temp=300;

if ~exist('nofig')
 nofig=0;
end

if ~exist('read')
 read=1;
end

if (read)
%%%%%%%%%% process windows
 fbw=0.5; % only applies to the first position component
 iwin=1;
 nwin=15;
 rcind=1 ; %cv index corresponding to the reaction coordinate
 kforce=100; % need force constant(s)
 dt=100 * 4 / 1000 / 1000 ; % time between frames in ns
%
 rc=zeros(1,nwin+1-iwin); % reaction coordinate
%
 for j=1:nwin+1-iwin
%   fprintf([' => Processing window ', num2str(j),'\n']);
   ncv=3;
   fname=['data/hist',num2str(j-1+iwin),'.dat'];
%   fprintf([' => Loading file ', fname,'\n']);
   d=load(fname) ;
   data=reshape(d,ncv,[])' ;
%
% extract the only relevant CV
   cvs=data(:,rcind) ;
%
   [niter(j),ncv]=size(cvs); % redefine ncv for rest of script
%
% sample limits and number of boxes
   ie   = niter
   ib   = 1;
   ib   = 10000 ;
%   ib=max(1,round(niter * 0.5));
%
   if (~exist('fc')); % free energy derivative
    fc=zeros(ie-ib+1,nwin);
    drc=zeros(ie-ib+1,nwin);
   end

% loop over all cvs, and compute PMF derivative
% first, need the reference position -- open cv.dat file :
   cvs0 = load(['cv',num2str(j-1+iwin),'.dat']);
   cvs0 = cvs0(rcind);
   rc(j) = cvs0;
%
% compute force
   dd=cvs(ib:ie) - cvs0; % displacement from restraint center
   drc(1:ie-ib+1)=dd ; %store displacement from RC
   fc(1:ie-ib+1,j) = kforce * ( min ( dd + 0.5*fbw , 0 ) + max ( dd - 0.5*fbw , 0 ) );
 end; % j -- windows
 read=0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot some statistics
% histograms, smoothed forces, correlations, series


figure(3) ; hold on; box on;
leg={};
for iwin=1:nwin ;

dr=drc(:,iwin) ;
f=fc(:,iwin);
n=length(fc);
t=[1:n]*dt ;

% running force average : 
%
skip=5000;
ms=4;
frunave=cumsum(f)'./[1:n]; % this is informative -- it shows how quickly/slowly the force converges
plot(t(1:skip:end),frunave(1:skip:end),char(styles(iwin)),'markersize',ms);
leg=[leg {['#',num2str(iwin),', x=',num2str(rc(iwin))]}];
xlabel('t(ns)')
ylabel('Force cumulative average (kcal/mol/A)') ;

end

legend(leg,-1)

set(gcf,'paperpositionmode','auto');
print(gcf,'-depsc2','3ngb-r71a-fc')


