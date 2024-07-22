% compare distances in sequence space : alignment distance (MATLAB), numerical distance, distance in PCA coordinates
%
addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')
addpath('../')
% phylogenetic tree

if ~exist('tab')
 tab=97 ; % table in the strains database (90, 95, 97-99)
end
% exclude residue positions with this fraction of missing res. (1 keeps all) :
pmissing=0.99;
qgrantham=0; % whether to use aa assignment from Grantham, 1974 (vs. Atchley 05, default)

flag='';
%flag='-1'; % select different PCA file
%flag='-1-2-5-6-8-9-11-12-13-16-17-18';

fname=['pca',num2str(tab),flag,'stem.mat'];
load(fname);

msamat=char( {msa.Sequence}' ) ;
% exclude missing residues (see above)
pmiss=mean((msamat=='-'),1);
imiss=find(pmiss>pmissing);
msamat(:,imiss)='' ; % delete those residue positions from alignment matrix
%showalignment(msamat);
% convert to coordinates
[coor,ndist]=aln2coor(msamat, qgrantham); % also compute pairwise distances

% (2) compute aligment distance
%adist=seqpdist(msa,'method','jukes-cantor','indels','pairwise-delete') ; ~90
adist=seqpdist(msa,'method','p-distance','indels','pairwise-delete') ; % ~97 for tab90 ; ~98 for tab 95/97 w/o grantham ; w/ grantham get close to 90% in some cases
%adist=seqpdist(msa,'method','alignment-score','indels','pairwise-delete','scoringmatrix',@pam250) ; % ~78

% (3) compute distances in PCA space ; expect to have a correlation to ndist, though not to adist
covmat=cov(coor,1) ; % default is to normalize by n-1 ; including 1 normalizes by n
% diagonalize :
if ~exist('qeig') ; qeig = true ; end
if (qeig)
 [v,evd]=eig(covmat) ;
 ev=diag(evd);
 qeig=false;
end
% project
[nseq,nvar]=size(coor) ;
vind=[0 -1 -2] + nvar ; % take the last evecs
pcc=msa2proj(msamat, v, vind, qgrantham);
% compute distances based on principal component coordinates (pcc)
pdist=zeros(1,nseq*(nseq-1)/2);
npc=3;
ind=1;
for i=1:nseq-1
 di=nseq-i ;
 pdist(ind:ind+di-1)=sqrt(sum(bsxfun(@minus,pcc(1:npc,i+1:end),pcc(1:npc,i)).^2,1));
 ind=ind+di;
end

corr(ndist', adist') % ~90%
corr(pdist', adist') % ~80%
corr(ndist', pdist') % ~80%
