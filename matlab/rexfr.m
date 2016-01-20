% replica diffusion analysis
% following the replica (environment)
% calculate approximate acceptance frequency (swap of replicas on adjacent cv sites)
%
% NOTE : HARDWIRE OFFSET (CORRECT LOGGING MISTAKE); works : maps agree !

close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
leg={};

rfreq=100;
nrep=32;
fnames={'rex.dat'};

tfac=rfreq*2/1000000; %output interval in nanoeconds (2fs time step; 1M fs in ns)

if (~exist('read')); read=1; end
if (read==1)

%
clear data;
for i=1:length(fnames)
 fname=char(fnames(i));
 if (i==1)
  data=load(fname);
 else
  data=[data; load(fname)];
 end 
end
%
data(:,1:2)=data(:,1:2)+1; % replica IDs start from 1
read=0;
end
%
% exchange freq.
nrep=max(data(:,2)); % will fail if some are zero!
nrep=32;
time=data(:,3)/rfreq;
numex=zeros(1,nrep-1);
%
for i=1:length(time)
 from=data(i,1);
 numex(from)=numex(from)+1;
end

tmax=time(end);
ar=numex/tmax*2;

% save
plot([0:nrep-2]/(nrep-2),ar,'k.-','linewidth',1.5);
%
box on;
ylabel('\it Acceptance Rate', 'fontsize',14);
xlabel(' \alpha', 'fontsize',14);
set(gcf, 'paperpositionmode', 'auto');
%
%print(gcf, '-dpsc', 'ar1_32rex.eps');
%print(gcf, '-djpeg100', 'ar1_32rex.jpg');


save ar1_32rex.dat ar numex -ascii

