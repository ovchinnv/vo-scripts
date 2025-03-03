addpath('/home/taly/scripts/matlab/sqlite3/mex-sqlite3')
% note that eigenvectors are invariant w.r.t reflection ; so two sets could be essentially the same, but up to -/+ reclections
% the correspondence between the 90 and 95% databases is not terrible

db='/home/taly/flurepo/flu_strains.db'
tab=95 ; % table in the strains database (90, 95, 97-99)

% exclude residue positions with this fraction of missing res. (1 keeps all) :
pmissing=0.99;
qgrantham=0; %whether to use aa assignmenef from Grantham, 1974 (vs. Atchley 05, default)
q2d=0;

table=['flu_strains',num2str(tab)];
%strains=sqlite3(db, ['select *, len(sequence) as slen, from flu_strains',num2str(tab)],' where nunk=0') ;
cmd= ['select rowid, *, length(sequence) as slen from ',table,' where Nunk=0 and sequence not LIKE ''%J%'''] ; % ignore sequences with 'J' because Matlab does not score J's
strains=sqlite3(db, cmd) ;
%
msa=struct([]);
nseq=length(strains);
% sadly, there is no way to vector-copy fields
for i=1:nseq
  msa(i).Header = strains(i).strain_id ; 
  msa(i).Name = strains(i).name ; 
  msa(i).Sequence = strains(i).msa ; 
end
%
msamat=char( {msa.Sequence}' ) ;
% exclude missing residues (see above)
pmiss=mean((msamat=='-'),1);
imiss=find(pmiss>pmissing);
msamat(:,imiss)='' ; % delete those residue positions from alignment matrix
%showalignment(msamat);

% convert to coordinates
coor=aln2coor(msamat,qgrantham);

cave=mean(coor,1) ;
covmat=cov(coor,1) ; % default is to normalize by n-1 ; including 1 normalizes by n
% subtract average sequence
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

%plot(rmsf); % fluctuations

% diagonalize :
if ~exist('qeig') ; qeig = true ; end
if (qeig)
 [v,evd]=eig(cvmat) ;
 ev=diag(evd);
% check diagonalization :
 dcheck = v*evd*inv(v) - cvmat ;
 max(abs(dcheck(:))); % get 1e-13 ; 
 qeig=false;
end

close all ;
% project
vind=[0 -1 -2] + nvar * nres ; % take the last evecs
pcc=msa2proj(msamat, v, vind, qgrantham);

if ~exist('q2d') ; q2d = false ; end
ms=5 ;
mc='k';
%mc='red';
mt='.' ;
if (q2d)
 scatter(pcc(1,:), pcc(2,:), ms, mc, mt) ; view(2);
else
 scatter3(pcc(1,:), pcc(2,:), pcc(3,:), ms, mc, mt) ; view(3);
end
leg={'all'};
hold on;

% add projection of subsets :
%pcahost
pcatype;
%pcayears; % not very illuminating

% fix axes
if (~qgrantham) ; l=50; else ; l=20 ; end

%axis([-l l -l l -l l]);

xlabel('PCA1'); xlim([-30 50]);
ylabel('PCA2'); ylim([-50 40]);
zlabel('PCA3'); zlim([-30 30]);

campos=[-509 -55 342];

set(gca,'CameraPosition',campos);
set(gcf, 'paperpositionmode','auto');
print(gcf,'-depsc2', 'ha-pca.eps');

