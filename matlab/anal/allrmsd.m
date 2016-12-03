% plot RMSD
addpath('~/scripts/matlab');
close all;

load('rmsds.mat') ;
hc=figure('position',[200 200 500 300]); hold on; box on;
lc=figure('position',[200 200 500 300]); hold on; box on;
bb=figure('position',[200 200 500 300]); hold on; box on;

handles=[ hc lc bb ];
styles={'r' 'g' 'b'};
d=20;

for ibeg=1:3

 j=0;
 for i=ibeg:3:length(labels)
  j=j+1;
  figure(handles(ibeg));
  rmsd=allrmsd(i,:);
  time=[1:length(rmsd)]*40/1000;
  rmsds=smooth2(time,rmsd,d);
  plot(time,rmsds,char(styles(j)))
 end
 legend(labels(ibeg:3:end),2,'interpreter','none')
 xlabel('\it t(ns)');
 ylabel('\it RMSD(\AA)','interpreter','latex')
 ylim([0 2.5]);
 set(gcf,'paperpositionmode','auto');
 print(gcf,'-depsc2',[char(labels(ibeg+6)),'.eps']);
end
