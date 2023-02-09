function CR=calc_covar(x,y,z,w);
% compute covariance matrix :
 [natom,nframes]=size(x);
 if (~exist('w'))
  sw=ones(size(x));
 else
  sw=sqrt( repmat( w, 1, nframes ));
 end
%
 X=[ sw.*x ; sw.*y ; sw.*z; ]' ;% concatenate coordinates ; the variables must be in different columns; the observations in different rows
%
 disp(['Computing covariance matrix ...']);
 fprintf('Found %d frames and %d atoms ...\n',nframes, natom);
 CX=cov(X); % note that this is weighted covariance matrix (unless weights are 1)

% compute variances for normalization
 VX=var(X,0); % likewise, a scaled variance
% test:
%SX=sqrt(VX);
%W=CX./(SX'*SX); % to remove mass weighting
%W2=corrcoef(X);
%pcolor(W-W2) ;shading interp ; colorbar % this difference is 1e-15 or 1e-5 depending on how variance is normalized
%return

% combine correlation components of each atom
 CR=zeros(natom,natom);
% ====== standard combination procedure
 for i=1:3
  indi = (i-1)*natom + 1 : i*natom ;
% below, only considering correlations between same components
% this is the traditional way of computing correlations between atoms
   CR=CR+CX(indi,indi);
%  for j=1:3
%   indj = (j-1)*natom + 1 : j*natom ; 
%   CR=CR+CX(indi,indj);
%  end
 end
% ====
 return
% omit normalization, because it can be done from CR later ; on the other hand, it cannot be undone without SR (diagonals)
% same for variances (diagonal of CR) :
 VR=zeros(1,natom);
 for i=1:3 % over components
  VR=VR+VX( (i-1)*natom + 1 : i*natom );
 end

 SR=sqrt(VR); % these are also root-mean-square fluctuations
 NR=SR'*SR;   % normalization matrix
 CRX=CR./NR ; % correlation coefficient matri
