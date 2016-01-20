% plot cell dimensions to decide when the system cell has been quilibrated
% 

inds=[6:20];

nind=length(inds) ;

celld=zeros(0,4); %cell vectors

for ind=inds
 f=['scratch/emr',num2str(ind),'npt.xstfile']; % file name
 d=load(f) ; %read file
% make sure we only take the latest run's data : 
 i=size(d,1) ;
 tstep = d(i,1);
 for j=i-1:-1:1
  tprev = d(j,1);
  if (tprev >= tstep)
   break;
  else
   tstep=tprev;
  end
 end
 d=d(j+1:end,[1,2,6,10]); % take data below this line
 celld = [celld ; d];
%
% ind
% d
%"===="
end

% now anneal time indices and compute time
i=size(celld,1);
dstep = celld(2,1)-celld(1,1) ; % assuming this is correct;
for j=2:i
 celld(j,1)=celld(j-1,1) + dstep ;
end

% compute time (ns)
dt = 2 / 1000000 ;
t=celld(:,1) * dt ;  %steps => 2fs /step =>ns

% plot
figure ; hold on;
plot(t, celld(:,2),'r+' ) ; %xsize (NOTE : xsize = ysize)
%plot(t, celld(:,3),'go' ) ; %ysize
%plot(t, celld(:,4),'bs' ) ; %zsize

ylabel('L(Ang)');
xlabel('t(ns)')
pause(0);

% area : 
celld(:,5) = celld(:,2).* celld(:,3);

nlipid = 75.5 ;% average between two leaflets
celld(:,5) = celld(:,5) / nlipid ; % note : this includes the area occupied by protein; so will be an large UPPER bound


