% crossing frequency as a function of time
% requires that vlog*.m be executed prior.
% basically, number of reflections per ns
fid=10;

sbin=10000; % bin size
nbin= ceil ( (max(time)-min(time))/sbin );
xf=hist(time,nbin);

figure(fid);
tfac=1000000/2; %conversion from simulation timestep (2fs) to ns
tt=[0:length(xf)-1]*sbin/tfac;
d=5;
xfs=smooth2(tt,xf,d);
plot(tt, xf*tfac/sbin, 'k', tt, xfs*tfac/sbin,'r');
xlabel('\it t(ns)', 'FontSize',15);
ylabel('\it # reflections', 'FontSize',15);
%
save xs_k10new.mat tt xf xfs
