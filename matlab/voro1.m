
his=load('his1.dat');
%prevhis=load('voro6.dat');
%his=his-prevhis;

[m,n]=size(his);
nstep=his(end,:);
%his=his(n+1:2*n,1:n);
his=his(1:n,1:n);

% fix diagonals if zero
du=diag(his,1);  du=min(1,du); du=1-du; %ones in place of zeros; elsewhere -- zeros
dl=diag(his,-1); dl=min(1,dl); dl=1-dl;
his2=diag(du,1)+diag(dl,-1) ;
% add matrix with fixes
his=his+1*his2;

nstep(:)=max(max(his));%100; % ad-hoc
nstep(:)=1 ;
%%%%%%%%%%%%%  reconstruct his matrix
r=zeros(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=1;% allow collisions between replicas separated by at most d-1 cells
%
for i=1:n
% compute indices
  ibeg = max(1,i-d) ; iend=min(i+d,n) ; inds=[ibeg:iend]; 
  r(i,inds)= his(inds,i)./nstep(inds)';
  r(i,i)= - sum(his(i,inds))/nstep(i);
end
%%%%%%%%%%%%%% set c(1)=1 %%%%%%%%%%%
c=zeros(n,1);
c(1)=1;
f=zeros(n,1);
f(:)=f(:)-r(:,1)*c(1);
c(2:end)=r(1:end-1,2:end)\f(1:end-1);
c=abs(c); % badly scaled matrices can erroneously produce small negative probabilities
%%%%%%%%%%%%%%%% solve for FE
beta=0.59582;
f=-beta*log(c);
close;

alpha=[0:length(f)-1]; alpha=alpha/alpha(end);
d=1;
fs=smooth2(alpha,f,d);

figure('position',[200,200,600,250]);
%plot(alpha,fs,'r*-');hold on;
plot(alpha,f,'r*-');hold on;

