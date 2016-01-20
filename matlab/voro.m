
his=load('his3.dat');
%his0=load('voro130.dat');
%his=his-his0;

[m,n]=size(his);
nstep=his(end,:);
%his=his(n+1:2*n,1:n);
his=his(1:n,1:n);

nstep(:)=max(nstep); % ad hoc

% fix diagonals if zero
du=diag(his,1);  du=min(1,du); du=1-du; %ones in place of zeros; elsewhere -- zeros
dl=diag(his,-1); dl=min(1,dl); dl=1-dl;
his2=diag(du,1)+diag(dl,-1) ;
% add matrix with fixes
his=his+his2;

%%%%%%%%%%%%%  reconstruct his matrix

r=zeros(n);
for i=1:n
  r(i,:)= his(:,i)./nstep(:);
  r(i,i)= - sum(his(i,:))/nstep(i);
end
%%%%%%%%%%%%%% set c(1)=1
c=zeros(n,1);
c(1)=1;
f=zeros(n,1);
f(:)=f(:)-r(:,1)*c(1);
c(2:end)=r(1:end-1,2:end)\f(1:end-1);
%%%%%%%%%%%%%%%% solve for FE
beta=0.59582;

f=-beta*log(c);
close;

alpha=[0:length(f)-1]; alpha=alpha/alpha(end);
d=2;
fs=smooth2(alpha,f,d);
lw=1.2;

figure('position',[200,200,600,250]);
plot(alpha,fs,'r.-', 'linewidth',lw,'markersize',14);hold on;

%legend(leg,2);
xlim([0 1]);
%ylim([-15 22]);
box on;
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
xlabel('\it \alpha', 'fontsize',14);
%
%
%
%set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'vfe.eps');
%