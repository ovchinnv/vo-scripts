addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')
addpath('../')
% note that eigenvectors are invariant w.r.t reflection ; so two sets could be essentially the same, but up to -/+ reclections
% the correspondence between the 90 and 95% databases is not terrible

if ~exist('tab')
 tab=95 ; % table in the strains database (90, 95, 97-99)
end
% exclude residue positions with this fraction of missing res. (1 keeps all) :
pmissing=0.99;
qgrantham=0; % whether to use aa assignmenef from Grantham, 1974 (vs. Atchley 05, default)
q2d=0;

if ~exist('flag','var')
 flag='ha2';
%flag='-1'; % select different PCA file
% flag='-1-2-5-6-8-9-11-12-13-16-17-18stem';
end

fname=['pca',num2str(tab),flag,'.mat'];
load(fname);

msamat=char( {msa.Sequence}' ) ;
% exclude missing residues (see above)
pmiss=mean((msamat=='-'),1);
imiss=find(pmiss>pmissing);
msamat(:,imiss)='' ; % delete those residue positions from alignment matrix
%showalignment(msamat);
% convert to coordinates
coor=aln2coor(msamat, qgrantham);

cave=mean(coor,1) ;
covmat=cov(coor,1) ; % default is to normalize by n-1 ; including 1 normalizes by n
% subtract average sequence
nseq=size(coor,1);
coor=coor-repmat(cave,nseq,1) ;
% compute covariance
cvmat = coor' * coor / nseq ;
vmat = diag(cvmat) ; % variance
% max(abs(cvmat-covmat)) % same thing
% compute fluctuations (ala RMSF)
nres=size(msamat,2);
fluc=zeros(nres, nres);
nvar=numel(coor)/numel(msamat);% number of cordinates per residue
for i=1:nvar % number of coordinates
 indi = i : nvar : nvar*nres ;
 fluc = fluc + cvmat(indi, indi) ; % combine fluctuations ; treat the five components as orthogonal, which makes them (technically) uncorrelated
end
rmsf=sqrt(diag(fluc)); % take square root of the diagonal
nrm=rmsf*rmsf'; % outer product
%pcolor(fluc./nrm) ; shading interp; % not sure whether (and how) the covariance here is meaningful

%figure;
%plot(rmsf); % fluctuations
%return
% diagonalize :
if ~exist('qeig') ; qeig = true ; end
if (qeig)
 [v,evd]=eig(cvmat) ;
 ev=diag(evd);
% check diagonalization :
 dcheck = v*evd*inv(v) - cvmat ;
 max(abs(dcheck(:))) % get 1e-13 ; 
 qeig=false;
end

% project
vind=[0 -1 -2] + nvar * nres ; % take the last evecs
pcc=msa2proj(msamat, v, vind, qgrantham);

q2d=0;
if ~exist('q2d') ; q2d = false ; end
ms=5 ;
mc='black';
%mc='red';
mt='.' ;
if (q2d)
 scatter(pcc(1,:), pcc(2,:), ms, mc, mt) ; view(2);
else
 scatter3(pcc(1,:), pcc(2,:), pcc(3,:), ms, mc, mt) ; view(3);
 zlabel('PC3');
end
xlabel('PC1');
ylabel('PC2');
leg={'all'};
hold on;

% add projection of subsets :
%pcahost
pcatype;
%pcayears; % not very illuminating

% fix axes
l=30;
%axis([-l l -l l -l l]);
