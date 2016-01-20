% 4.09: caculate work from the cv.dat and force.dat files
% The results show that the two methods of computing work agree (i.e. this, external, and the one in CHARMM, internal)
% use Simpson's Rule

close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if (~exist('read'))
 read=1;
end
if (read==1)


%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'force_22.dat_48','force_23.dat_48','force_24.dat_48','force_25.dat_48','force_26.dat_48',...
        'force_27.dat_48','force_28.dat_48','force_29.dat_48','force_30.dat_48','force_31.dat_48'};
%
clear data;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  fc=load(fname);
 else
  fc=[fc; load(fname)];
 end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cv=load('./cv.dat');
[ncv,m]=size(cv);

[n,m]=size(fc);
niter=n/ncv;

r=zeros(ncv,m,niter); 
f=zeros(ncv,m,niter);

row=1;
for i=1:niter
 for j=1:ncv
  r(j,:,i)=cv(j,:); % replicate cv data (cv's do not change in this calc.)
  f(j,:,i)=fc(row,:); row=row+1;
 end
end
read=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind=[1:2];

dr=r(ind,2:end,:)-r(ind,1:end-1,:); 

% compute work (a particularly inelegant coding)
a=zeros(ncv,m,niter);
b=zeros(ncv,m,niter);
c=zeros(ncv,m,niter);
d=zeros(ncv,m,niter);
alpha=zeros(ncv,m,niter);
beta=zeros(ncv,m,niter);
gamma=zeros(ncv,m,niter);
work=zeros(m+1-mod(m,2),niter);

for i=2:2:m-1
 a(:,i,:)=dr(:,i,:).*dr(:,i,:);
 b(:,i,:)=dr(:,i,:).*dr(:,i-1,:);
 c(:,i,:)=dr(:,i-1,:).*dr(:,i-1,:);
 d(:,i,:)=dr(:,i,:)+dr(:,i-1,:);
 alpha(:,i,:)=(2*a(:,i,:)+b(:,i,:)-c(:,i,:))./dr(:,i,:);
 beta(:,i,:)=d(:,i,:).^3./b(:,i,:);
 gamma(:,i,:)=(2*c(:,i,:)+b(:,i,:)-a(:,i,:))./dr(:,i-1,:);
 work(i+1,:)=work(i-1,:) + squeeze ( sum ( alpha(:,i,:).*f(:,i+1,:) + beta(:,i,:).*f(:,i,:) + gamma(:,i,:).*f(:,i-1,:) , 1) )';
end
 rep=[1:2:m];

if (mod(m,2)==0) % if the number of points is even, there is a hanging half-interval; appx. with cubic and integrate analytically on this half-interval
 i=m;
 alpha(:,i,:)=dr(:,i-1,:).*(3-dr(:,i-1,:)./(dr(:,i-1,:)+dr(:,i-2,:)));
 beta (:,i,:)=dr(:,i-1,:).*(3+dr(:,i-1,:)./dr(:,i-2,:));
 gamma(:,i,:)=-dr(:,i-1,:).^3./( dr(:,i-2,:) .* (dr(:,i-2,:)+dr(:,i-1,:)));
 work(i,:)=work(i-1,:) + squeeze ( sum ( alpha(:,i,:).*f(:,i,:) + beta(:,i,:).*f(:,i-1,:) + gamma(:,i,:).*f(:,i-2,:) , 1) )';

 rep=[rep m];
end
work=work(rep,:)/6;
%return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bin and plot average work curves %%%%%%%%%
% code taken from cfe.m
bsize=100;
ib=1;
ie=niter;

%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
fe=zeros(nsample,m/2+1-mod(m,2));
j=1;
for i=ib:bsize:ie
 fe(j,:)=mean(work(:, i : i+min(bsize-1,ie-i) ),2);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

figure; hold on;box on;

for i=1:nsample
 plot(rep,fe(i,1:end),[char(styles(mod(i-7,length(styles))+1)),'x'], 'linewidth', 0.1)
end 

%

fave=mean(fe,1);
fstd=std(fe,1);

%mean
%plot([1:m],fave,'k-*','linewidth',3);
%std
%plot([1:nrep-1],fave+fstd,'k:','linewidth',3);
%plot([1:nrep-1],fave-fstd,'k:','linewidth',3);
%leg=[leg {['Average']}];


legend(leg,2);
box on;
ylabel('\it Free Energy (kcal/mol)', 'fontsize',14);
xlabel('\it Replica ', 'fontsize',14);

%set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'fe.eps');
%print(gcf, '-djpeg100', 'fe.jpg');

work=[rep' fe'];

