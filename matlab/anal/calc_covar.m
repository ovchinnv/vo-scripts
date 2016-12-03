function CRX=calc_covar(x,y,z,w);
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
 CX=cov(X);

% compute variances for normalization
 VX=var(X,0);
% test:
%SX=sqrt(VX);
%W=CX./(SX'*SX);
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
% same for variances : 
 VR=zeros(1,natom);
 for i=1:3 % over components
  VR=VR+VX( (i-1)*natom + 1 : i*natom );
 end

 SR=sqrt(VR); % these are also root-mean-square fluctuations
 NR=SR'*SR;   % normalization matrix
 CRX=CR./NR ;
