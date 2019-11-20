% examine forces from double half-harmonic window simulations
%
close all;

if ~exist('styles')
 styles={'r-*','g-x','b-o','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
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
   [niter,ncv]=size(cvs); % redefine ncv for rest of script
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

iwin=3 ; 
nhist = 25 ; % number of bins for force histogram

dr=drc(:,iwin) ;
f=fc(:,iwin);
n=length(fc);
t=[1:n]*dt ;

figure(1) ; plot(t,f) ; % this is noisy ; good to superpose smoothed force
dsmooth=200 ; 
fs=smooth2(t,f,dsmooth);
figure(1) ; hold on ; plot(t,fs,'r','linewidth',2) ;
xlabel('t(ns)');
ylabel('Force (kcal/mol/A)');
legend('Instantaneous',[num2str(dsmooth*dt),'-ns average']);

[h,x]=hist(fs,nhist) ; % can histogram f or fs (fs should tend to Gaussian at convergence)
dx=x(2)-x(1);
figure(2) ; hold on ; box on ; bar(x,h,'facecolor','white'); xlabel('kcal/mol/A'); % hist(f) not too useful because most samples fall within flat botton window
gfit=exp(-0.5*((x-mean(fs))/std(fs)).^2)*n*dx/sqrt(2*pi*std(fs)^2);
plot(x,gfit,'r.-','linewidth',2)
sk=skewness(fs); % crude measure of deviation from gaussianity
% compute scalar product ; comes out too high
hg = h*gfit' / norm(h) / norm(gfit) ;
% compute difference between pdfs ; comes out not sensitive enough !
hdiff = sqrt ( sum ( ((h-gfit)/n/dx).^2 ) *dx ) ;
legend(['skewness=',num2str(sk)],'Gaussian',-1);
%legend(['similarity=',num2str(hg)],'Gaussian',-1);
%legend(['RMSD=',num2str(hdiff)],'Gaussian',-1);
%
% running force average : 
%
frunave=cumsum(f)'./[1:n]; % this is informative -- it shows how quickly/slowly the force converges
figure(3); plot(t,frunave); 
xlabel('t(ns)')
ylabel('Force cumulative average (kcal/mol/A)') ;
%
% correlation plots -- smoothed and unsmoothed plots
% unsmoothed :
ds=f; 
dsm=mean(ds); dss=std(ds);
dn=(ds-dsm)/dss;
dhat=fft(dn);
dhat2=dhat.*conj(dhat);
dc=ifft(dhat2);
dcl=floor(length(ds)/2);
dc=dc(1:dcl); dc=dc/dc(1) ; % normalization
% smoothed :
ds=fs;
dsm=mean(ds); dss=std(ds);
dn=(ds-dsm)/dss;
dhat=fft(dn);
dhat2=dhat.*conj(dhat);
dc2=ifft(dhat2);
dcl2=floor(length(ds)/2);
dc2=dc2(1:dcl2); dc2=dc2/dc2(1) ; % normalization
% find the point at which dc crosses zero; set block length to alpha times that interval
alpha = 1.;
k=find(dc<0,1); k=floor(k*alpha) ; iblock=min(k, floor(length(dn)/2)); t(iblock)
k=find(dc2<0,1); k=floor(k*alpha) ; iblock2=min(k, floor(length(dn)/2)); t(iblock2)
% plot
figure(4) ; hold on; box on
plot(t(1:dcl), dc,'k',t(1:dcl2), dc2,'r');
set(gca,'xscale','log')
xlabel('t(ns)');


