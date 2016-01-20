
skip=10;
p=100 ;
fname='rmsd-ehee';

d=load([fname,'.dat']);

dt=40 ; %interval between frames in picoseconds

t=(1:length(d)) * dt / 1000 ; %ns

plot(t,d);

xlabel('\it t(ns)', 'fontsize',14);
ylabel('\it RMSD(Ang)', 'fontsize',14);
set(gca,'fontsize',14)
title('RMSD from initial structure');
legend(fname);

print(gcf,'-depsc2',[fname,'.eps']);

