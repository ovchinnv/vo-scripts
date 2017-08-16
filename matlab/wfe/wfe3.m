% compute FE from double half-harmonic window simulations
% In this version (2) we compute the fe derivative from the average force 
% obtained form histogram data
%
close all;

if ~exist('styles')
 styles={'r.-','g-x','b-o','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
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
 nwin=15;
 nbox =2; % number of statistical samples
 rcind=1 ; %cv index corresponding to the reaction coordinate
 kforce=100; % need force constant(s)

 clear df ;
 rc=zeros(1,nwin); % reaction coordinate
%
 for j=1:nwin
%   printf([' => Processing window ', num2str(j),'\n']);
   ncv=3;
   fname=['data/hist',num2str(j),'.dat'];
%   printf([' => Loading file ', fname,'\n']);
   d=load(fname) ;
   data=reshape(d,ncv,[])' ;
%
% extract the only relevant CV
   cvs=data(:,rcind) ;
%
   [niter,ncv]=size(cvs); % redefine ncv for rest of script
%
   if (~exist('df')); % free energy derivative
    df=zeros(nwin,ncv,nbox);
   end
% load sample weights from tempering
   wgtname=['pmfsep1_',num2str(j),'.temp_series.txt'];
   w=load(wgtname);
   logw=w(:,end); % log of weights

% sample limits and number of boxes
   ie   = niter
   ie   = 11000 ;
   ib   = 1;
   ib=round(niter * 0.25);
   bsize=ceil( (ie-ib+1)/nbox);
   nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));

% loop over all cvs, and compute PMF derivative
% first, need the reference position -- open cv.dat file :
   cvs0 = load(['cv',num2str(j),'.dat']);
   cvs0 = cvs0(rcind);
   rc(j) = cvs0;
%
   fprintf([' => Computing average force...\n']);
   for k=1:ncv
    ibox=1;
    for i=ib:bsize:ie
     inds=[i:i+min(bsize-1,ie-i)] ; % indices for this box
% compute force
     dd=cvs(inds,k) - cvs0(k); % average displacement from restraint center
     force = kforce * ( min ( dd + 0.5*fbw , 0 ) + max ( dd - 0.5*fbw , 0 ) );
     logwgt=logw(inds);
     logwgt=logwgt-max(logwgt);
     wgt=exp(logwgt);
     df(j, k, ibox) = (force')*wgt / sum(wgt) ;

     ibox=ibox+1;
    end % over sample boxes
   end % over cvs
%
 end; % j -- windows
 read=0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we know that only one cv is present
dfc=reshape(df(:,1,:), nwin, nbox)' ;
rc0=rc ; % copy
if (0)
% append 0th point, at which we assume df to be zero (stable eq. simulation)
% however, this may not be true is the initial condition is not perfect, or the CVs are not perfect
% in some cases it is best to omit 0th point
 dfc=[zeros(nbox,1) dfc];
 rc0 =[2*rc0(1) - rc0(2), rc0]; % assume uniform interval
end

% trapezoidal rule
% compute center derivative
dfc = 0.5 * ( dfc(:,1:end-1) + dfc(:,2:end) );
%integral
% note that I am assuming that the force is zero in the first [0th] cv value
% this may not be true
fe = [ zeros(nbox,1)  cumsum(-dfc.*(ones(nbox,1)*diff(rc0)), 2) ];

%============= PLOT FE =============
if ~nofig
 close all;
 figure('position',[200,200,450,350]); hold on; box on;
end
%
for i=1:nbox
 plot(rc0,fe(i,1:end),[char(styles(mod(i-1,length(styles))+1)),''], 'linewidth', lw)
end
%
fave=mean(fe,1);
fstd=std(fe,1);
%mean
%plot(rc0,fave,'k--','linewidth',lw);
%std
%plot(alpha,fave+fstd,'k--','linewidth',1);
%plot(alpha,fave-fstd,'k--','linewidth',1);
%leg=[leg {['Average']}];

legend(leg,4);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it x', 'fontsize',14);
%
%xlim([0 1]);
%ylim([0 9]);
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'wfe2.eps');
minfe=min(fe,[],2);
fe-repmat(minfe,1,size(fe,2))
mean(ans)
%pause(10)
