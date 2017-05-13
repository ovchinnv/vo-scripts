% plot energy vs temperature
%
close all;
kboltz=1.98e-3 ;

irun=1;
erun=35;
irep=1;
basename=['../pmfsep1_',num2str(irep),'.temp_'];

files={};
for i=irun:erun
 files=[files {[basename,num2str(i),'.series.txt']}];
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
t1=0;
t2=1000;
i1 = floor(t1 / dt) ; i1=find(step>i1,1,'first') ;
i2 = floor(t2 / dt) ; i2=find(step<=i2,1,'last') ; % from step # find the index #

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

