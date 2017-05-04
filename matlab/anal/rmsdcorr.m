% autocorrelation function of RMSD
addpath('~/scripts/matlab');
close all;

load('rmsds.mat') ;
h=figure('position',[100 100 750 525]); hold on; box on; grid on;

handles=[ h ];
styles={'r' 'g' 'b', 'm', 'c', 'k', 'k--'};

d=20;

j=0;

inds=1:2:length(labels)
leg={};

for i=inds
  j=j+1;
  figure(handles(1));
  rmsd=allrmsd(i,:);
  time=[1:length(rmsd)]*20/1000 * dcdstep ;
  rmsds=smooth2(time,rmsd,d/dcdstep);

% compute approximate number of samples using correlation function :
  ds=rmsd;
  dsm=mean(rmsd);
  dss=std(rmsd);
  dn=ds-dsm;
  dn=(ds-dsm)/dss;
  dhat=fft(dn);
  dhat2=dhat.*conj(dhat);
  dc=ifft(dhat2);
%
%  dc=xcorr(dn,'unbiased'); alternative matlab calc
  dcl=floor(length(ds)/2);
%  dc=autocorr(dn,dcl)
  dc=dc(1:dcl);
  dc=dc/dc(1) ; % normalization
% find the point at which dc crosses zero; set block length to alpha times that interval
  alpha = 1.;
  k=find(dc<0,1); k=floor(k*alpha) ; iblock=min(k, floor(length(dn)/2));
  plot(time(1:length(dc))-time(1), dc, char(styles(j)))
  leg=[ leg {[char(labels(i)),', t_{corr}=',num2str(time(iblock)),'ns']} ] ;
end

legend(leg,1,'interpreter','tex')
xlabel('\it t(ns)', 'fontsize',14);
ylabel('\it RMSD(\AA)','interpreter','latex','fontsize',14)
%ylim([-1 1  ]);
%xlim([0 300/2]);
set(gcf,'paperpositionmode','auto');
print(gcf,'-depsc2','g120-gbis-rmsd-acf.eps');
