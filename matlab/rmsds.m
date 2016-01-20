%
styles={'r-','g-','b-','m-','c-','k-','y-'};

clear rms dom ndom
close all;

fname='rmsd_bpgm_nvt.dat';

[dom,d]=textread(fname,'%s%[^\n]', 'bufsize', 100000); % bufsize : max string length

ndom=length(dom);


if (ndom>0) 
 rms=[str2num(char(d(1)))];
end
%
for i=2:ndom
 rms=[rms ; str2num(char(d(i))) ];
end
%
rms=rms(:,1:end-1); % leave off last column
niter=length(rms);
%
% associate domain names with indices
for i=1:ndom
 eval([char(dom(i)),'=',num2str(i),';']);
end

%
ind=[core cap active loop1 helix2 hirms1 hirms2];

time=[1:niter]/100*2;
figure; hold on; box on; lw=1.5;

%d=60;
d=5;
for i=1:length(ind);
 plot(time,smooth2([1:niter],rms(ind(i),:),d), [char(styles(mod(i-1,length(styles))+1))], 'linewidth', lw);
% plot(rms(ind(i),:));
end

legend(dom(ind),-1);

xlabel('\it t(ns)');
ylabel('\it RMSD(x,x_T)');
set(gca,'linewidth',lw,'FontSize',9);

set(gcf,'paperpositionmode','auto');
%print(gcf,'-depsc2','rmsd_bpgm_nvt.eps');

