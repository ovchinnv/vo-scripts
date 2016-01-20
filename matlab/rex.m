% replica diffusion analysis
% following the replica (environment)
%
% NOTE : HARDWIRE OFFSET (CORRECT LOGGING MISTAKE); works : maps agree !

close all;
styles={'r-','g-','b-','m-','c-','k-','r--','g--','b--','m--','c--','k--'};
leg={};

rfreq=100;
nrep=32;

tfac=rfreq*2/1000000; %output interval in nanoeconds (2fs time step; 1M fs in ns)

fnames={'rex.dat'};
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

data=sortrows(data,3);
nrep=max(data(:,2))+1;
data(:,1:2)=data(:,1:2)+1; % replica IDs start from 1
time=data(:,3)/rfreq;
rmap =zeros(nrep,time(end)+1); % map following CV, and seeing which replicas pass through
rmapa=rmap;                    % map following replica
% generate map
tprev=0*[1:nrep]+1;
tpreva=0*[1:nrep]+1;
rmap(:,1) =[1:nrep];
rmapa(:,1)=[1:nrep];
for i=1:length(time)
% fprintf([num2str(i),' of ',num2str(length(time))]);
 from=data(i,1);
 to  =data(i,2);
 who1=rmap(from,tprev(from));
 who2=rmap(to,  tprev(to));
 when=time(i)+1;
%
 rmap(from,tprev(from)+1:when-1)=who1;
 rmap(to,tprev(to)+1:when-1)    =who2;
%
 rmap(from,when)                =who2;
 rmap(to,when)                  =who1;
%
 tprev(from)=when; tprev(to)=when;
% adjoint
 rmapa(who1,tpreva(who1)+1:when-1)=rmapa(who1,tpreva(who1));
 rmapa(who2,tpreva(who2)+1:when-1)=rmapa(who2,tpreva(who2));
 rmapa(who1,when)=to;
 rmapa(who2,when)=from;
%
 tpreva(who1)=when; tpreva(who2)=when;
%
end
%
for i=1:nrep
 rmap(i,tprev(i)+1:end)=rmap(i,tprev(i));
 rmapa(i,tpreva(i)+1:end)=rmapa(i,tpreva(i));
end
%
map=rmap(:,end)';
mapa=rmapa(:,end)';

t=[1:time(end)+1]*tfac;

mapa2=load('rex.map');
mapa2=mapa2(2,:)+1;
