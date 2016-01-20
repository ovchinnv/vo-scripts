% compute FE using BAR
%
close all;

btol=1e-6;
format long ;
temp=300;
kboltz = 1.987191e-3;
bet=1/temp/kboltz ;


ffw='rab-m5a0-1.mat';
fbk='rab-m5a1-0.mat';

load(ffw);
defw=dall;
nfw=nsamp;

load(fbk);
debk=dall;
nbk=nsamp;

clear dall nsamp;


nw=length(nfw); % number of windows

sif=zeros(nw,1); % statistical inefficiency
sib=zeros(nw,1); % statistical inefficiency (backward)
dA=zeros(nw+1,1); % free energy as a function of window
dAe=zeros(nw+1,1); % error in the above free energy


%%%%%%%%% first, need to produce uncorrlated data %%%%%%%%%
alpha = 0.5;
%alpha = 1;
for i=1:nw
% forward
% compute correlation length to get block size
 de=defw(1:nfw(i),i);
 dm=mean(de);
 dv=var(de);
 ds=sqrt(dv);
 dn=(de-dm)/ds;
 dhat=fft(dn);
 dhat2=dhat.*conj(dhat);
 dc=ifft(dhat2);
 dc=dc(1:floor(length(dc)/2)) ;
 dc=dc/dc(1) ; % should be same as below

% find the point at which dc crosses zero; set block length to alpha times that interval
 ind=find(dc<0); ind=floor(ind(1)*alpha) ; iblock=min(ind, floor(length(de)/2));
 fprintf([num2str(i),': Computed forward block length : ',num2str(iblock),'\n']);
 sif(i)=iblock;

% plot(dc,'r-'); hold on
%return
% ibeg=1; % can be any number
 ibeg=floor(1+1.*iblock); % can be any number
 def=bet*de(ibeg:sif(i):end) ; nf=length(def);
%
%
% backward
% compute correlation length to get block size
% note that order of backward windows is reversed with respect to the order of forward ones
 
 j=nw-i+1;

 de=debk(1:nbk(j),j);
 dm=mean(de);
 dv=var(de);
 ds=sqrt(dv);
 dn=(de-dm)/ds;
 dhat=fft(dn);
 dhat2=dhat.*conj(dhat);
 dc=ifft(dhat2);
 dc=dc(1:floor(length(dc)/2)) ;
 dc=dc/dc(1) ; % should be same as below

% find the point at which dc crosses zero; set block length to alpha times that interval
 ind=find(dc<0); ind=floor(ind(1)*alpha) ; iblock=min(ind, floor(length(de)/2));
 fprintf([num2str(j),': Computed backward block length : ',num2str(iblock),'\n']);
 sib(j)=iblock;
%
% subsample uncorrelated time series
%
 ibeg=floor(1+1.*iblock); % can be any number

 deb=bet*de(ibeg:sib(j):end) ; nb=length(deb);

% BAR iterations
 C=0 ;
 dC=Inf;
 while dC > btol
% bennett equations :
  logsize = log(nb/nf); % note : "forward" corresponds to simulations sampling state _0

  dAest = log ( sum ( 1./(1+exp(deb+C)) ) ) - log ( sum ( 1./(1+exp(def-C)) ) ) + C - logsize ;

  Cold=C ;
  C=dAest+logsize;
  dC=abs(C-Cold);
 end
 dA(i+1) = dA(i) + dAest ;
% compute std error :
 ddAe(i+1) =      ( ( mean ( 1./(1+exp(def-C)).^2 ) / mean ( 1./(1+exp(def-C) ) )^2 - 1 ) / nf + ...
                    ( mean ( 1./(1+exp(deb+C)).^2 ) / mean ( 1./(1+exp(deb+C) ) )^2 - 1 ) / nb ) ;
 dAe(i+1) = dAe(i) + ddAe(i+1);

end


dA  = dA / bet  ;% scale by thermal energy
dAe = sqrt(dAe) / bet ;% standard error

for i=2:nw+1
 fprintf([' ',num2str(i),' : ',num2str(dA(i)),' +/- ', num2str(dAe(i)),'\n']);
end

plot(dA,'r.-') ;
errorbar([1:nw+1], dA, dAe);


