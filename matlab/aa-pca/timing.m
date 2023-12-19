% plot MSA time cost vs number of sequences
maxseqs=[10:10:100];
for maxseq=maxseqs
 qmsa=1;
% pca;
% time(i)=t;
end
clf;
plot(maxseqs,time(time>0),'k*-')
hold on ;
xp=1.5 ; % exponent
plot(maxseqs, maxseqs.^xp / maxseqs(1)^xp * time(find(time,1)) ,'r' )

set(gca, 'xscale','log')
set(gca, 'yscale','log')

