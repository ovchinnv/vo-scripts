% compare quasiharmonic entropy
%

close all;

fname='entropy-ca.mat' ;
load(fname);

labels={'3BNC60/HC' 'CH103/HC' 'PGT121/HC'} ;
leg={'Mature','Intermediate','Germline'};

tshc=reshape(tsclass(:,1),3,3);

b=bar(tshc');

ylabel('-TS^{cg}_{quasi}')
set(gca,'xticklabel',labels)

legend(leg,-1);

set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc','quasi_hc.eps');

% same for light chain

labels={'3BNC60/LC' 'CH103/LC' 'PGT121/LC'} ;

tshc=reshape(tsclass(:,2),3,3);

figure;
b=bar(tshc');

ylabel('-TS^{cg}_{quasi}')
set(gca,'xticklabel',labels)


legend(leg,-1);

set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc','quasi_lc.eps');
