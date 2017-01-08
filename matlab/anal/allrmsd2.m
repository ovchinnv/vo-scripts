% plot RMSD
addpath('~/scripts/matlab');
close all;

load('rmsds.mat') ;
hc=figure('position',[200 200 500 300]); hold on; box on;
lc=figure('position',[200 200 500 300]); hold on; box on;
bb=figure('position',[200 200 500 300]); hold on; box on;

handles=[ hc lc bb ];
styles={'r' 'g' 'b', 'm', 'c', 'k', };
d=50;

ibeg=1 ; %3bnc
%ibeg=10 ; %ch103
%ibeg=19 ;%pgt121

iend=ibeg+3*3-1 ; % look at one class at a time

for offset=0:2

 j=0;
 for i=ibeg+offset:3:iend+offset
  j=j+1;
%  i
  figure(handles(offset+1));
  rmsd=allrmsd(i,:);
  time=[1:length(rmsd)]*40/1000;
  rmsds=smooth2(time,rmsd,d);
  plot(time,rmsds,char(styles(j)))
 end
 legend(labels(ibeg+offset:3:iend+offset),2,'interpreter','none')
 xlabel('\it t(ns)');
 ylabel('\it RMSD(\AA)','interpreter','latex')
 ylim([0 3]);
 set(gcf,'paperpositionmode','auto');
 print(gcf,'-depsc2',[char(labels(ibeg+offset)),'rmsd.eps']);
end
