% compute msas for all tables
tabs=[90 95 97 98]; % omit 99 for now
tabs=[97 ]; %
%
for tab=tabs
 qmsa=1;
 pca;
 time(i)=t;
end
clf;
plot(maxseqs,time(time>0),'k*-')
hold on ;
xp=1.5 ; % exponent
plot(maxseqs, maxseqs.^xp / maxseqs(1)^xp * time(find(time,1)) ,'r' )

set(gca, 'xscale','log')
set(gca, 'yscale','log')

