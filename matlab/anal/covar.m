% covariance matrices
%
name='3bnc60glt';
flag='-ca-hc';
%flag='-ca-lc';

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
colorbar('h'); % not needed

% load sequence alignment data
load malign.mat
load malign-lc.mat
showdom2;

% print
set(gcf,'paperpositionmode','auto');
print(gcf, '-depsc2', [name,'-covar',flag,'.eps']);
% save data
save([name,'-covar',flag,'.mat'], 'CRX');
