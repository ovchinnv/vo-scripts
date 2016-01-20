% preprocessor to calculate average forces using WHAM (Dave Mihn implementationi using Benoit's scheme)

if (~exist('read'))
 read=1;
end
if (read==1)


%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'force1.dat'};
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

cv=load('./cv1.dat');
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
kf=250*[ones(1,32)];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%% now, theta contains the raw time series %%%%%%%%%%%%%%%%%%%%

samples=squeeze(theta);
bins=cv';
npt=32;
x=linspace(0,1,npt);
[fe,x]=wham2(samples,bins,kf,300,x); % need to scale because wham script uses different units

plot(x,fe)
