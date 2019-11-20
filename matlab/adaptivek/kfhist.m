% plot temperature history
%
close all;
kboltz=1.98e-3 ;
temp=300;

files={ 'adaptive.series.txt'};

%files={};
%for i=10:19
% files=[files {['hhh_4.0_01.temp',num2str(i),'.series.txt']}];
%end

for i=1:length(files)
 file=char(files(i));
 if (i==1)
  d=load(file);
 else
  d=[d; load(file)];
 end 
end

step=d(:,1);
param=d(:,2);
dual=d(:,3);

dt = 1 * 4 / 1000000 ; % plugin_freq x timestep (fs) x (ns / fs) => to convert to time (ns)

plot(step*dt, param,'k-','linewidth',0.5) ; box on; hold on;
xlabel('\it t(ns)', 'Fontsize',14) ;
ylabel('\it k_{fc}(kcal/mol/A^2)', 'Fontsize',14);

%ylim([250 500]);

set(gcf,'position',[100 100 800 200]);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'kfc_hist');

% RMS
figure;
plot(step*dt, 2*dual*kboltz*temp,'k-','linewidth',0.5) ; box on; hold on;
xlabel('\it t(ns)', 'Fontsize',14) ;
ylabel('\it \rho^2(A^2)', 'Fontsize',14);

set(gcf,'position',[100 100 800 200]);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'dual_hist');

% compute and plot density (this can also be done from the restart file itself of course)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid :
npt=100;
[pdf,bins]=hist(param,npt); pdf=pdf/length(step);
figure ; hold on; grid on; box on
plot(bins, pdf,'k.-')
% superpose powerlaw
%plot(bins,4*bint.^0,'r')

set(gca, 'xscale','log')
set(gca, 'yscale','log')

xlabel('\it k_{fc}(kcal/mol/A^2)', 'fontsize',14)
ylabel('\it PDF', 'fontsize',14);
legend('pdf','power-law fit');
set(gcf, 'paperpositionmode','auto') ;
%print(gcf, '-depsc2', 'pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mapped pdf
npb=100;
[pdf2,bins2]=hist(log(param),npb); pdf2=pdf2/length(step);
figure ; hold on; grid on; box on
plot(bins2, pdf2,'k.-')
% superpose powerlaw
%plot(binb,0.013*binb.^-1,'r')

set(gca, 'xscale','log')
set(gca, 'yscale','log')

xlabel('\it \beta(kcal/mol)^{-1}', 'fontsize',14)
ylabel('\it PDF(\beta) a.u.', 'fontsize',14);
legend('pdf','power-law fit');
set(gcf, 'paperpositionmode','auto') ;
%print(gcf, '-depsc2', 'mapped_pdf');


