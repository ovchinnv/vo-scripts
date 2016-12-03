% covariance matrices
names={'3bnc60glt', '3bnc60at', '3bnc60t'};
flags={'-ca-hc', '-ca-lc' };
flags={'-cg-hc', '-cg-lc' };
clrs={'b','g','r'};
aligns={'malign.mat' 'malign-lc.mat'};

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

% compute covariance matrix :
CRX=calc_covar(xall,yall,zall,mass); % optional mass-weighting

figure('position',[ 200 200 500 500 ]);

pcolor(CRX) ; shading flat ; hold on ; box on ;
xlabel('\it residue');
ylabel('\it residue');
caxis([-0.5 1]);
colorbar('h'); % not needed

% load sequence alignment data
load(align);
%load malign-lc.mat
showdom2;

% print
set(gcf,'paperpositionmode','auto');
print(gcf, '-depsc2', [name,'-covar',flag,'.eps']);
% save data
save([name,'-covar',flag,'.mat'], 'CRX');

end %jj
end %ii
