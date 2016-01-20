%
close all;

ener=load('ener_all.dat');

enav=mean(ener,1);
enav=enav-enav(1);

nrep=length(enav);
alpha=([1:nrep]-1)/(nrep-1);

d=3;

enavs=smooth2(alpha,enav,d);

lw=1.2;

figure('position',[200,200,600,250]);
plot(alpha,enavs,'k.-', 'linewidth',lw, 'markersize',15);

%legend(leg,2);
xlim([0 1]);
ylim([-15 22]);
box on;
%
xlabel('\alpha','fontsize',14);
ylabel('\it <E(\alpha)> (kcal/mol) ','fontsize',14);

%
set(gcf, 'paperpositionmode', 'auto');
print(gcf, '-dpsc', 'pote.eps');
