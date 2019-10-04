% plot RMSD
addpath('~/scripts/matlab');
close all;

load('rmsds.mat') ;

figure('position',[200 200 500 300]); hold on; box on;

handles=[ bb ];

styles={'r' 'g' 'b', 'm', 'c', 'k', };
d=50;

 for i=1:length(labels)
  rmsd=allrmsd(i,:);
  time=[1:length(rmsd)]*40/1000;
  rmsds=smooth2(time,rmsd,d);
  plot(time,rmsds,char(styles(i)))
 end
 legend(labels(i),2,'interpreter','none')
 xlabel('\it t(ns)');
 ylabel('\it RMSD(\AA)','interpreter','latex')
 ylim([0 3.5]);
 set(gcf,'paperpositionmode','auto');
 print(gcf,'-depsc2','rmsd.eps']);

