% calculate average forces using Kastner & Thiel 05
% for Kwangho: nonuniform force constants

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
 fname=['./',char(fnames(i))];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note: there are hardwirings below
% force constants:
%NOTE: these vary as a function of replica (NOT CV!)
%9-10: 125; 
%11-12: 175; 
%13-14: 250; 
%15-18: 300; 
%19-20: 250;  
%21-22: 174; 
%23-24: 125
kf=100*[ones(1,48)]; kf(9:10)=125; kf(11:12)=175; kf(13:14)=250; kf(15:18)=300; kf(19:20)=250; kf(21:22)=174; kf(23:24)=125;
%kf=[ones(1,45), 10*ones(1,5),1]; % originally per CV!
%
%divide by kf to get z-theta(x):
dz=zeros(size(f));
theta=zeros(size(f));
%for j=1:ncv
% dz(j,:,:)=f(j,:,:)/kf(j);
%end
%
for j=1:m
 dz(:,j,:)=f(:,j,:)/kf(j);
end
%
% subtract bias center;
theta=r-dz;
% calculate variance
varz=var(theta,1,3);
avez=mean(theta,3);
% beta
beta=0.59582;

%%%%%%%%%%%%%%%% create an interpolated grid %%%%%%%%%%%%%%%
%
np=48;
r0=squeeze(r(:,:,1));
rnew=zeros(ncv,np);
% interpolate CV:
for i=1:ncv
 rnew(i,:)=interp1([0:m-1]/(m-1),r0(i,:),[0:np-1]/(np-1));
end
% go over all bins and all z values:
fnew=zeros(ncv,np,m);
w=zeros(np,m);
for i=1:np
 fprintf([num2str(i),'\n']);
 for j=1:m
  fnew(:,i,j) = (rnew(:,i)-avez(:,j))./varz(:,j)*beta - kf(j)*(rnew(:,i)-r0(:,j));
% using means:
  w(i,j)=prod(1./sqrt(varz(:,j)*2*pi).*exp(-0.5*(rnew(:,i)-avez(:,j)).^2./varz(:,j)));
% using z:
%  w(i,j)=prod(1./sqrt(varz(:,j)*2*pi).*exp(-0.5*(rnew(:,i)-r0(:,j)).^2./varz(:,j)));
 end
% normalize weights
  w(i,:)=w(i,:)/sum(w(i,:));
end
% average force over all bins
fave=zeros(ncv,np);
for i=1:ncv
 fave(i,:)=sum(squeeze(fnew(i,:,:)).*w(:,:),2);
end
%
%
% trapezoid rule:
dr=rnew(:,2:end)-rnew(:,1:end-1); 
fc=0.5*(fave(:,1:end-1)+fave(:,2:end)); 

dw=squeeze(sum(dr.*fc,1));

work=zeros(np,1);
for i=2:np
 work(i)=work(i-1)+dw(i-1);
end

%plot(cv(2,:)-cv(1,:), work)
plot(work)

set(gcf, 'paperpositionmode', 'auto');
%print(gcf, '-dpsc', 'fe.eps');
%print(gcf, '-djpeg100', 'fe.jpg');

%save corrected_fe.dat -ASCII work
