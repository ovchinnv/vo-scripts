%
names={'3bnc60glt', '3bnc60at', '3bnc60t'};

names={'ch103t', 'ch103-i3.2t', 'ch103ucat'};

names={'gl121t', '3h109lt', 'pgt121t'};

flags={'-ca-hc', '-ca-lc' };
flags={'-cg-hc', '-cg-lc' };
clrs={'b','g','r'};
aligns={'malign.mat' 'malign-lc.mat'};

%c='r'; % mature
%c='g'; % p60a
%c='b'; % germline
%flag='-cg-hc';
%flag='-cg-lc';

for ii=1:3
 name=char(names(ii));
 clr=char(clrs(ii));
 for jj=1:2
  flag=char(flags(jj));
  align=char(aligns(jj));

qpdb=0 ; qpsf=0 ; qdcd=0; clear xall;% to reread

pdbfile=[name,flag,'.pdb'];

struc ;
% define some selections
select ;

% dcd trajectory names :
dcdnames={ [name,flag,'.dcd'] } ;
%
% read dcds :
traj;


figure('position',[100 300 900 300]); hold on; box on;
leg={};
lw=2;
% load sequence alignment data
%load malign.mat ;% for heavy chains
%load malign-lc.mat ; % for light chains
load(align);
% indicate domains
showdom;

% compute rmsf for five subtrajectories
nrun=2500;
fluc=zeros(natom,5);
for i=1:5
 irun=(i-1)*nrun+1;
 erun=i*nrun;
 fluc(:,i)=calc_rmsf(xall(:,irun:erun), yall(:,irun:erun), zall(:,irun:erun));
end
flucall=calc_rmsf(xall,yall,zall);
plot(flucall,clr,'linewidth',lw) ;hold on ; box on ;
flucstd=std(fluc');
plot(flucall+flucstd,[clr,'--'],'linewidth',1) ;hold on ; box on ;
plot(flucall-flucstd,[clr,'--'],'linewidth',1) ;hold on ; box on ;


xlabel('\it residue');
ylabel('$\it RMSF(\AA)$', 'interpreter','latex');

legend(['Average RMSF : ',name],'RMSF +/- 1xSTDev');

% print
xlim([0 125])

set(gcf,'paperpositionmode','auto');
print(gcf, '-depsc2', [name,'-rmsf',flag,'.eps']);
% save data
save([name,'-rmsf',flag,'.mat'], 'fluc', 'flucall','flucstd','resid');

end %jj
end %ii
