% 2/2013 :  look at mean force profiles
% The results show that the two methods of computing work agree (i.e. this, external, and the one in CHARMM, internal)

close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%% load multiple force files %%%%%%%%%%%%
fnames={'force5.dat' 'force6.dat' 'force7.dat' 'force8.dat' 'force9.dat' 'force10.dat' 'force11.dat' 'force12.dat' 'force13.dat' ...
 'force14.dat' 'force15.dat' 'force16.dat' 'force17.dat' 'force18.dat'  'force19.dat' 'force20.dat' 'force21.dat' 'force22.dat' ...
 'force23.dat' 'force24.dat' 'force25.dat' 'force26.dat' 'force27.dat'  'force28.dat' 'force29.dat' 'force30.dat' 'force31.dat' 'force32.dat' ...
 'force33.dat' 'force34.dat' 'force35.dat' 'force36.dat' 'force37.dat' 'force38.dat' ...
 'force39.dat' 'force40.dat' 'force41.dat' 'force42.dat' 'force43.dat' 'force44.dat' ...
 'force45.dat' 'force46.dat' 'force47.dat' 'force48.dat' 'force49.dat' 'force50.dat' 'force51.dat' 'force52.dat' ...
 'force53.dat' 'force54.dat' 'force55.dat' 'force56.dat' 'force57.dat' 'force58.dat' 'force59.dat' 'force60.dat' ...
 'force61.dat' 'force62.dat' 'force63.dat' 'force64.dat' 'force65.dat' ...
 'force66.dat' 'force67.dat' 'force68.dat' 'force69.dat' 'force70.dat' ...
 'force71.dat' 'force72.dat' 'force73.dat' 'force74.dat' 'force75.dat' 'force76.dat' 'force77.dat' 'force78.dat' 'force79.dat' ...
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist('read'))
 read=1;
end
read=1
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
[niter,nstring]=size(fc);
ndim=3; % number of force components
niter=niter/ndim;

f=zeros(ndim,nstring,niter);

row=1;
for i=1:niter
 for j=1:ndim
  f(j,:,i)=fc(row,:); row=row+1;
 end
end
read=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=[1:2]; %including curvature at index 2
ind=[1:1]; %planar force only
%ind=[2:2]; %curvature only
%ind=[3:3]; %normal force only (not a true profile)

ftotal=squeeze(sum(f(ind,:,:),1)); % total force

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bin and plot curves %%%%%%%%%
% code taken from cfe.m
%
nbox =2;
ib   =1;
ib   =300;
%ib   =2500;
ie   =niter;
bsize=ceil( (ie-ib+1)/nbox);

%bsize=90

%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
fav=zeros(nsample,nstring);
j=1;
for i=ib:bsize:ie
 fav(j,:)=mean(ftotal(:, i : i+min(bsize-1,ie-i) ),2);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

figure; hold on;box on;


favs=zeros(nsample,nstring);
for i=1:nsample
% smooth force 
 dfilter=1;
 favs(i,:)=smooth2([1:nstring],fav(i,:),dfilter);
%
% plot(fav(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', 2)
 plot(favs(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'.'], 'linewidth', 2)
end 

fave=mean(fav,1);
fstd=std(fav,1);

% integrate force to obtain FE profile
fe=zeros(nsample,nstring);


fe(:,1)=0 ;
for i=2:nstring
 fe(:,i) = fe(:,i-1) - 0.5 * ( favs(:,i-1) + favs(:,i) ) ;
end

figure; hold on;box on;
for i=1:nsample
 plot(fe(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', 2)
end 


