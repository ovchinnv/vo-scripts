% 4.09: caculate work from the cv.dat and force.dat files
% The results show that the two methods of computing work agree (i.e. this, external, and the one in CHARMM, internal)

close all;

styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%% load multiple force files %%%%%%%%%%%%
 basename='data/force'
 basext='.dat';

 fnames={};
 
 ibeg=1;
 iend=72;
 for i=ibeg:iend
  fnames=[fnames {[basename,num2str(i),basext]}];
 end

 fnames
%
%fnames={'data/force1.dat'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('read'))
 read=1;
end
read=1
if (read==1)
%
clear data;
for i=1:length(fnames)
%for i=1:7
 fname=char(fnames(i));
 fprintf('%s %s\n','Reading file: ',fname);
 if (i==1)
  fc=load(fname);
 else
  fc=[fc; load(fname)];
 end 
end
%
[niter,m]=size(fc);
ncv=3;
niter=niter/ncv;

cv=[]; for i=1:ncv ; cv=[cv ; 1:m]; end ;

[ncv,m]=size(cv);

r=zeros(ncv,m,niter);
f=zeros(ncv,m,niter);

row=1;
for i=1:niter
 for j=1:ncv
  r(j,:,i)=cv(j,:);
  f(j,:,i)=fc(row,:); row=row+1;
 end
end
read=0;
%f(2,:,:)=-f(2,:,:); % flip sign of curvature force (debugging)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=[1:2]; %including curvature at index 2
ind=[1:1]; %planar force only
%ind=[2:2]; %curvature only
%ind=[3:3]; %normal force only (not a true profile)

dr=r(ind,2:end,:)-r(ind,1:end-1,:); 
fc=0.5*(f(ind,1:end-1,:)+f(ind,2:end,:)); 

dw=-squeeze(sum(dr.*fc,1));
%dw=dw(:);

work=zeros(m,niter);
for i=2:m
 work(i,:)=work(i-1,:)+dw(i-1,:);
end

% perpendicular force correction
workp=zeros(m,niter);

for i=1:m
 workp(i,:)=-squeeze(f(3,i,:))';
end
workp(:,:)=workp(:,:)-ones(m,1)*workp(1,:); % reset to zero

%work=work+workp;
%work=workp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bin and plot average work curves %%%%%%%%%
% code taken from cfe.m
% NOTE : string fixed as of simulation #46 (iteration 104)

ie   =niter ;
ib   =1 ;
ib=round(niter * 0.5);
nbox =1;
bsize=ceil( (ie-ib+1)/nbox);

lw=1.2;

%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
fe=zeros(nsample,m);
j=1;
for i=ib:bsize:ie
 fe(j,:)=mean(work(:, i : i+min(bsize-1,ie-i) ),2);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

figure('position',[200,200,600,250]); hold on; box on;
alpha=[0:m-1]; alpha=alpha/alpha(end);

for i=1:nsample
 plot(alpha,fe(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', lw)
end 

fave=mean(fe,1);
fstd=std(fe,1);

%mean
%plot(alpha,fave,'rx-','linewidth',lw);
%std
%plot(alpha,fave+fstd,'k--','linewidth',1);
%plot(alpha,fave-fstd,'k--','linewidth',1);
leg=[leg {['Average']}];

%legend(leg,2);
box on;
grid on ;
set(gca,'ytick',[-2:1:7]);
ylabel('\it F(\alpha) (kcal/mol)', 'fontsize',14);
%xlabel('\it \alpha', 'fontsize',14);
xlabel('\it Reaction coordinate (\alpha)', 'fontsize',14);
%
% label states
text(0.04,0,'\it \bf M5B');
text(0.88,1.5,'\it \bf M5A');
text(0.48,4,'\it \bf FE barrier');

%
xlim([0 1]);
%ylim([-15 20]);
set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'fe_contrib.eps');
print(gcf, '-dpsc', 'fe-rm5b-rm5a-32.eps');
%
%
return
% plot evolution of FE drop
iw=10;
dfe=zeros(1,nsample);
for i=1:nsample
 dfe(i)=min(fe(i,end-iw+1:end) - fe(i,1:iw));
end
figure(3);
plot(dfe,'k-x');
ylabel('\it Endpoint Free Energy Difference(kcal/mol)', 'fontsize',14);
xlabel('\it Iteration \# ', 'fontsize',14);

