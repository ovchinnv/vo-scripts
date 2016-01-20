close all;
%
% read m tensor

if (~exist('read'))
 read=1;
end

if (read==1)

uuref=load('mmat_ref.dat');

uu=load('mmat.dat');
[m,n]=size(uu);

nstring=m/n;

uu=reshape(uu',[n,n,nstring]);
%uave=mean(uu,3);
uave=uu(:,:,end);

%d=max(abs(diag(uave)));
%uave=uave/d;

%for i=1:n
% uave2(i,:)=uave(i,:) / uave(i,i);
% uave2(i,i)=0;
%end
%sum(abs(uave2))/(n-1)

uuref=reshape(uuref',[n,n,nstring]);
uref=uuref(:,:,1);
%d=max(max(abs(uref)));

du=uave-uref;
% normalize this somehow !
for i=1:n
 for j=1:n
  du(i,j)=du(i,j) / sqrt ( uref(i,i) * uref(j,j) );
 end
end

read=0;
end ;% read

aind=[1:n]/3;

pcolor(aind, aind, abs(du)); b=colorbar; set(b,'box','on'); %caxis([0 0.01]);
colormap gray; shading flat ;% faceted
c=colormap;
off=0.6; y=1/(1-off); x=y-1;
c= (c+x)/y ; colormap(c)
box on; set(gca,'TickDir','out','fontsize',10, 'fontangle','italic')
xlabel('Atom index');
ylabel('Atom index');
set(gcf, 'paperpositionmode','auto');
print(gcf, '-dpsc', 'mmat.eps');

mean(mean(du))

%pcolor(abs(uave-uref)); colorbar; caxis([0 0.01]);


%pcolor(uave); colorbar;
%pcolor(abs(uave)); colorbar;

%[v,e]=eig(uave);
%for i=4:n
% e(i,i)=1/e(i,i);
%end
%ui=v*e*v';

%max(max(abs(ui-uave)))

%pcolor(ui);
