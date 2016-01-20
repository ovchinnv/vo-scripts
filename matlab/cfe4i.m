% calculate average forces using Kastner & Thiel 05

if (~exist('read'))
 read=1;
end
if (read==1)


%%%%%%%%%% load multiple files %%%%%%%%%%%%
%fnames={'force_ss5.dat','force_ss6.dat','force_ss7.dat','force_ss8.dat','force_ss9.dat','force_ss10.dat','force_ss11.dat','force_ss12.dat'};
%fnames={'force_ss25.dat','force_ss26.dat','force_ss27.dat'};
fnames={'force_ss27.dat'};
%
clear data;
for i=1:length(fnames)
 fname=['../',char(fnames(i))];
 if (i==1)
  fc=load(fname);
 else
  fc=[fc; load(fname)];
 end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cv=load('../cv_4.dat');
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
kf=[ones(1,45), 10*ones(1,5),1];
%divide by kf to get z-theta(x):
dz=zeros(size(f));
theta=zeros(size(f));
for j=1:ncv
 dz(j,:,:)=f(j,:,:)/kf(j);
end
% subtract bias center;
theta=r-dz;
% calculate variance
varz=var(theta,1,3);
avez=mean(theta,3);
% beta
beta=0.59582;

%%%%%%%%%%%%%%%% create an interpolated grid %%%%%%%%%%%%%%%
% do the interpolation consistently with the original spacing distribution
%
ifact=4; % factor by which the number of points will be increased (roughly)
np=(m-1)*ifact+1;

r0=squeeze(r(:,:,1));
rnew=zeros(ncv,np);

% interpolate CV:
% intermediate points:
for j=1:ifact
  rnew(:,j:ifact:j+ifact*(m-1)-1) = ((j-1) * r0(:,2:end) + (ifact-j+1) * r0(:,1:end-1))/(ifact);
end
% endpoint
rnew(:,end) = r0(:,end);
% check:
%plot([0:m-1]/(m-1),r0(1,:),'kx'); hold on
%plot([0:np-1]/(np-1),rnew(1,:),'r.');
%return

% go over all bins and all z values:
fnew=zeros(ncv,np,m);
w=zeros(np,m);
for i=1:np
 fprintf([num2str(i),'\n']);
 for j=1:m
  fnew(:,i,j) = (rnew(:,i)-avez(:,j))./varz(:,j)*beta - kf'.*(rnew(:,i)-r0(:,j));
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

