% compute FE from double half-harmonic window simulations
%

fefile='wfe.mat' ;

kboltz=1.987191d-3;        % boltzmann constant
Temp=300;

if ~exist('nofig')
 nofig=0;
end

if ~exist('styles')
 styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
end
lw=1;
leg={};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%% load multiple files %%%%%%%%%%%%
 basename='./window'
 basext='.dat';

 fnames={};
 
 ibeg=72;
 iend=82;
 for i=ibeg:iend
  fnames=[fnames {[basename,num2str(i),basext]}];
 end

 fnames
%fnames={'window.dat'} ;
% specify window width (must match simulation file)
dwin=0.5 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('read'))
 read=1;
end
read=1 ;
if (read==1)
%
clear data;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  data=load(fname);
 else
  data=[data; load(fname)];
 end 
end
%
% separate averages from number of samples
nsamp=data(2:2:end,:) ;
delta=data(1:2:end,:) ;
%
read=0 ;
end
%
%
[niter,m]=size(delta);
% restraint centers for each replica
delta0=[0 0.5*ones(1,m-2) 1] ;

% local distance metric in terms of average ds
fac=[ 1 2*ones(1,m-2) 1 ] ;
dlocal=dwin./fac ; % dwin expressed in units of local metric

% sample limits and number of boxes
ie   =niter;
ib   =1;
%ib=round(niter * 0.2);
nbox =1;
%
dfe=zeros(nbox,m);
% loop over all replicas, and compute PMF derivative

for i=1:m
% consider only entries with nonzero samples
 inds=find(nsamp(ib:ie,i)>0)+ib-1; % sample for consideration
 if (isempty(inds)) % simulation too short to produce samples within window
  gam(:,i)=NaN; % will deal with this later
  warning(['Window ',num2str(i),' has no valid samples. Setting F''(',num2str(i),') to NaN.']);
  continue
 end
 nn=nsamp(inds,i);
 dd=delta(inds,i) - delta0(i); % average displacement from restraint center
 dbox = 0.5*dlocal(i) ; % half-width in the local metric
 dd = (dd+dbox)/(2*dbox) ; % normalize data to the interval [0 1]
%
% check to make sure there are no outliers
 if ( find ( abs(dd - 0.5) > 0.5 , 1) )
  error ['ERROR : one or more averages for replica ', num2str(i),' is outside the [0,1] interval. Aborting.'];
 end
% cumulative number of samples
 nnc=cumsum(nn);
 bsize=ceil( nnc(end)/nbox ); %block size
 nnc=nnc/bsize ; % normalize sample count by block size and look for block indices, 1, 2, 3 ... etc
%
 ibeg=1;
 for ibox=1:nbox
  if (i==1) ; leg=[leg {['iteration ',num2str(ibeg)]}]; end ;

  iend=find(nnc>ibox,1)-1 ; % last index
  if (isempty(iend))
   if ibox < nbox
    error ['End of sample reached in box ',num2str(ibox),' of ',num2str(nbox),' of replica ',num2str(m)];
    return
   else
    iend=length(nnc) ;%take the last sample
   end
  end % iend
% combined average
  ddcomb=sum(dd(ibeg:iend).*nn(ibeg:iend))/sum(nn(ibeg:iend)) ;
% define ibeg for next block (if any)
  ibeg=iend+1;
% solve for FE slope
  f=@(x) xave_s(x)-ddcomb ;
  g=fzero(f,0) ;
% stop if a problem with gamma
  if(isnan(g)) ; error('Cannot find gamma: got NaN'); return ; end
%
  gam(ibox,i) = g ;
%
 end % ibox

end % over all replicas
% compute dfe from gamma
dfe=gam*kboltz*Temp/dwin ; % scale to inter-replica distance and kBT
%
% deal with NaNs:
indnan=find(isnan(dfe)) ;
indok=setdiff([1:m],indnan) ;
for i=indnan
 dfe(i)=interp1(indok,dfe(indok),i);
end
%
% compute center derivative
dfec = 0.5 * ( dfe(:,1:end-1) + dfe(:,2:end) );
fe = [ zeros(nbox,1)  cumsum(dfec,2) ];

% save data in a file
alpha=[0:m-1]; alpha=alpha/alpha(end);
save(fefile, 'gam', 'dwin', 'kboltz', 'Temp', 'fe', 'alpha');
%============= PLOT FE =============
if ~nofig
 close all;
 figure('position',[200,200,450,350]); hold on; box on;
end
%
indi=2;
inde=29;

alpha2=(alpha(indi:inde) - alpha(indi)) / (alpha(inde)-alpha(indi));
%
for i=1:nbox
 fe2=fe(i,indi:inde)-fe(i,indi)
 plot( alpha2 ,fe2, [char(styles(mod(i-1,length(styles))+1)),'o'], 'linewidth', lw)
end
%
fave=mean(fe,1);
fstd=std(fe,1);
%mean
%plot(alpha,fave,'k--','linewidth',lw);
%std
%plot(alpha,fave+fstd,'k--','linewidth',1);
%plot(alpha,fave-fstd,'k--','linewidth',1);
leg=[leg {['Average']}];

legend(leg,2);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
%
xlim([0 1]);
%ylim([0 9]);
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'wfe.eps');
