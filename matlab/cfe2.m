% 4.09: caculate work from the cv.dat and force.dat files
% The results show that the two methods of computing work agree (i.e. this, external, and the one in CHARMM, internal)

%close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%% load multiple force files %%%%%%%%%%%%
fnames={'force_22.dat_48','force_23.dat_48','force_24.dat_48','force_25.dat_48','force_26.dat_48',...
        'force_27.dat_48','force_28.dat_48','force_29.dat_48','force_30.dat_48','force_31.dat_48'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('read'))
 read=1;
end
%
if (read==1)
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
%
cv=load('cv.dat');
ncv=2; % number of CV; this is specified by the user

[ncv,m]=size(cv);

[n,m]=size(fc);
niter=n/ncv;


r=zeros(ncv,m,niter); % (numpos, numrep, numiter)
f=zeros(ncv,m,niter);

row=1;
for i=1:niter
 for j=1:ncv
  r(j,:,i)=cv(j,:);
  f(j,:,i)=fc(row,:); row=row+1;
 end
end
read=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dr=r(:,2:end,:)-r(:,1:end-1,:); 
fc=0.5*(f(:,1:end-1,:)+f(:,2:end,:)); 

dw=squeeze(sum(dr.*fc,1));

work=zeros(m,niter);
for i=2:m
 work(i,:)=work(i-1,:)+dw(i-1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bin and plot average work curves %%%%%%%%%
% code taken from cfe.m
bsize=20;
ib=1;
ie=niter;

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

figure; hold on;box on;

for i=1:nsample
 plot(fe(i,1:end),[char(styles(mod(i-7,length(styles))+1)),'x'], 'linewidth', 0.1)
end 

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






