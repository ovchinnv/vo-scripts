% compute FE from double half-harmonic window simulations
% in this version, plot only the forces, not the FE integral
% to see which windows are converging
%

close all;

if ~exist('styles')
 styles={'r-x','g-v','b-s','m-*','c-','k-','r--','g--','b--','m--','c--','k--'};
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
 iwin=0; % can be 0 or 1 depending on whether the equilibrium point is included
 nwin=7;
% nsamples=4;
% [status, result]=system('grep "will quit" pmf3.out | tail -n1 | awk ''{print $3}'''); nsamples=str2num(result)-1 ; nsamples=nsamples-31 ; % screwed up counts due to crash
% [status, result]=system('grep "will quit" pmf2.log | tail -n1 | awk ''{print $3}'''); nsamples=str2num(result)-251 ;
 nbox=2; % number of statistical samples
 rcind=1 ; %cv index corresponding to the reaction coordinate
 kforce=100;

 clear df ;
 rc=zeros(1,nwin+1-iwin); % reaction coordinate
 for j=1:nwin+1-iwin;

   ncv=4;
   fname=['data/fbwin',num2str(j-1+iwin),'.dat'];
   dnew=load(fname) ;
   if exist('nsamples')
    d=dnew(1:2*ncv*nsamples,:); % i.e. 2*ncv lines per samples
   else
    d=dnew;
   end
%
   data=reshape(d,ncv,[])' ;
%
% extract the only relevant CV
   cvs=data(1:2:end,rcind) ;
   nsamp=data(2:2:end,rcind) ;
%
   [niter,ncv]=size(cvs); % redefine ncv for rest of script
%
   if (~exist('df'))
    df=zeros(nwin,ncv,nbox);
   end

% sample limits and number of boxes
   ie   = niter
   ib   = 1;
   ib=max(1,round(niter * 0.5));
   bsize=ceil( (ie-ib+1)/nbox);

% loop over all cvs, and compute PMF derivative
% first, need the reference position -- open cv.dat file :
   cvs0 = load(['cv',num2str(j-1+iwin),'.dat']);
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
     df(j, k, ibox) = mean (force) ;
     ibox=ibox+1;
    end % over sample boxes
   end % over cvs
%
 end; % j -- windows
 read=0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we know that only one cv is present
dfc=reshape(df(:,1,:), nwin+1-iwin, nbox)' ;
rc0=rc;
if (0)
% append 0th point, at which we assume df to be zero (stable eq. simulation)
% however, this may not be true is the initial condition is not perfect, or the CVs are not perfect
% in some cases it is best to omit 0th point
 dfc0=[zeros(nbox,1) dfc0];
 rc0 =[2*rc0(1) - rc0(2), rc0]; % assume uniform interval
end

%============= PLOT =============
if ~nofig
 close all;
 figure('position',[200,200,450,350]); hold on; box on;
end
%
for i=1:nbox
% plot(rc,dfc(i,1:end),[char(styles(mod(i-1,length(styles))+1)),''], 'linewidth', lw)
end
%
dfave=mean(dfc,1);
dfstd=std(dfc,1);
plot(rc,dfstd,'ko-','linewidth',lw);
% approximate error in the FE : 
['RMSD error : ',num2str(norm(dfstd)),' kcal/mol']

%
legend(leg,4);
box on;
ylabel('\it dF/dx(\alpha) (kcal/mol/A)', 'fontsize',14);
xlabel('\it x', 'fontsize',14);
%
%xlim([0 1]);
%ylim([0 9]);
%set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'dfe.eps');
%
% save force std to a text file :
%dfstd=dfstd'; % transpose
save -ascii dfe.dat dfstd
gwin ;% calculate windows to sample
