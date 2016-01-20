%
close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};

if (~exist('read'))
 read=1;
end

if (read==1)
%%%%%%%%%% load multiple files %%%%%%%%%%%%
fnames={'curv5.dat' 'curv6.dat' 'curv7.dat' 'curv8.dat' 'curv9.dat' 'curv10.dat' 'curv11.dat' 'curv12.dat' 'curv13.dat' ...
 'curv14.dat' 'curv15.dat' 'curv16.dat' 'curv17.dat' 'curv18.dat'  'curv19.dat' 'curv20.dat' 'curv21.dat' 'curv22.dat' ...
 'curv23.dat' 'curv24.dat' 'curv25.dat' 'curv26.dat' 'curv27.dat'  'curv28.dat' 'curv29.dat' 'curv30.dat' 'curv31.dat' 'curv32.dat' ...
 'curv33.dat' 'curv34.dat' 'curv35.dat' 'curv36.dat' 'curv37.dat' 'curv38.dat' ...
 'curv39.dat' 'curv40.dat' 'curv41.dat' 'curv42.dat' 'curv43.dat' 'curv44.dat' ...
 'curv45.dat' 'curv46.dat' 'curv47.dat' 'curv48.dat' 'curv49.dat' 'curv50.dat' 'curv51.dat' 'curv52.dat' ...
 'curv53.dat' 'curv54.dat' 'curv55.dat' 'curv56.dat' 'curv57.dat' 'curv58.dat' 'curv59.dat' 'curv60.dat' ...
 'curv61.dat' 'curv62.dat' 'curv63.dat' 'curv64.dat' 'curv65.dat' ...
 'curv66.dat' 'curv67.dat' 'curv68.dat' 'curv69.dat' 'curv70.dat' ...
 'curv71.dat' 'curv72.dat' 'curv73.dat' 'curv74.dat' 'curv75.dat' 'curv76.dat' 'curv77.dat' 'curv78.dat' 'curv79.dat' ...
 'curv81.dat' 'curv82.dat' 'curv83.dat' 'curv84.dat' 'curv85.dat' 'curv86.dat' 'curv87.dat' 'curv88.dat' 'curv89.dat' ...
 'curv91.dat' 'curv92.dat' 'curv93.dat' 'curv94.dat' 'curv95.dat' 'curv96.dat' 'curv97.dat' 'curv98.dat' 'curv99.dat' ...
 };

fnames={'curv94.dat'};

%
clear crv;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  crv=load(fname);
 else
  crv=[crv; load(fname)];
 end 
end
crv=crv(:,2:end)'; % leave off inxed column
read=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nstring,niter]=size(crv);
%
%
nbox =1;
ib   =1;
%ib   =300;
%ib   =3000;
%ib   =niter;
ie   =niter;
bsize=ceil( (ie-ib+1)/nbox);

%bsize=90

%
leg={};
nsample=floor((ie-ib+1)/bsize)+sign(mod(ie-ib+1,bsize));
cav=zeros(nsample,nstring);
j=1;
for i=ib:bsize:ie
 cav(j,:)=mean(crv(:, i : i+min(bsize-1,ie-i) ),2);
 j=j+1;
 leg=[leg {['iteration ',num2str(i-1)]}];
end 

%figure; hold on;box on;
figure('position',[200,200,600,250]); hold on; box on;


alpha=([1:nstring]-1)/(nstring-1);

cavs=zeros(nsample,nstring);
for i=1:nsample
% smooth curvature
 dfilter=3;
 cavs(i,:)=smooth2(alpha,cav(i,:),dfilter);
%
% plot(cav(i,1:end),[char(styles(mod(i-1,length(styles))+1)),'x'], 'linewidth', 2)
 plot(alpha,cavs(i,:),['*',char(styles(mod(i-1,length(styles))+1))], 'linewidth', 2)
end 

box on;
ylabel('\kappa (Ang^{-1})', 'fontsize',14);
xlabel('\it \alpha ', 'fontsize',14);

axis([0 1 0 4.5]);
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'curvp.eps');
%print(gcf, '-djpeg100', 'fe64_p1.jpg');
