% plot energy vs temperature
%
close all;
kboltz=1.98e-3 ;


files={};
for i=10:19
 files=[files {['eeh_1.0_01.temp',num2str(i),'.series.txt']}];
end

for i=1:length(files)
 file=char(files(i));
 if (i==1)
  d=load(file);
 else
  d=[d; load(file)];
 end 
end

step=d(:,1);
temp=d(:,2);
ener=d(:,3);

dt = 10 * 4 / 1000000 ; % plugin_freq x timestep (fs) x (ns / fs) => to convert to time (ns)

% plot for prtions of trajectory
i1 = floor(300 / dt)+1 ; i1=find(step>i1,1)-1 
i2 = round(600 / dt)+1 ; i2=find(step>i2,1)-1 ; % from step # find the index #

plot(temp(i1:i2), ener(i1:i2), 'k.','linewidth',0.5) ; box on; hold on;

% plot the folded portion :
i3 = i2+1 ;
i4 = length(step) ; % remainder
plot(temp(i3:i4), ener(i3:i4), 'r.','linewidth',0.5)



xlabel('\it T(K)', 'Fontsize',14);
ylabel('\it E(kcal/mol)', 'Fontsize',14) ;

xlim([250 650]);

set(gcf, 'paperpositionmode','auto') ;
print(gcf, '-depsc2', 'evt');

