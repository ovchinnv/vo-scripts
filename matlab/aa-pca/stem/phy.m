addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')
addpath('../')
% phylogenetic tree

if ~exist('tab')
 tab=90 ; % table in the strains database (90, 95, 97-99)
end
% exclude residue positions with this fraction of missing res. (1 keeps all) :
pmissing=0.99;
qgrantham=0; % whether to use aa assignment from Grantham, 1974 (vs. Atchley 05, default)
q2d=0;

flag='';
%flag='-1'; % select different PCA file
%flag='-1-2-5-6-8-9-11-12-13-16-17-18';

fname=['pca',num2str(tab),flag,'stem.mat'];
load(fname);

dmethod='jukes-cantor';
dist=seqpdist(msa,'method','jukes-cantor','indels','pairwise-delete') ;
ptree=seqlinkage(dist,'single',msa);
view(ptree)