function [xall,yall,zall,rmsd]=bestfit(xall,yall,zall,xref,yref,zref,wgt,inds)
% perform alignment using procrustes
% does not implement mass weighting correctly but error is small for most purposes
% not quite the same as VMD/CHARMM because weights cannot be passed into the procrastes routine
% thus, until Procrustes gets fixed, only a uniformly-weighted fit is possible
%
 if (~exist('inds')) ; inds=[1:length(xref)]' ;
 else
  disp(['Superposition uses a subset of atoms : ',num2str(inds')]);
 end
%
 if (~exist('wgt') || isempty(wgt)) ; wgt=ones(size(xref)); end
 wgt=wgt/sum(wgt(inds));
 swgt=sqrt(wgt);
 oswgt=1./swgt(:);
 comref=[wgt(inds)'*xref(inds) wgt(inds)'*yref(inds) wgt(inds)'*zref(inds)];
% mass weighting after translation to COM
 Xref=[ swgt(:).*( xref(:) - comref(1) )  swgt(:).*( yref(:) - comref(2) )  swgt(:).*( zref(:) - comref(3) ) ];
%
 nall=(size(xall,2)); %; number of frames
 rmsd=zeros(1,nall);
 for i=1:nall
  if (~mod(i,100))
   disp(['==> Aligning coordinate set ',num2str(i)]);
  end
%
  com=[wgt(inds)'*xall(inds,i) wgt(inds)'*yall(inds,i) wgt(inds)'*zall(inds,i)];
%
  X=[ swgt(:).*(xall(:,i)-com(1)) swgt(:).*(yall(:,i)-com(2)) swgt(:).*(zall(:,i)-com(3)) ]; % also subtract COM
  [d,Xnew,trans]=procrustes(Xref(inds,:), X(inds,:), 'scaling',0,'reflection',0);
% conpute _all_ rotated coordinates (without translation component)
  Xnew=X*trans.T ; % pure rotation (only in this case recover correct magnitudes below)
% compute RMSD :
%  rmsd(i)=norm(Xref(inds,:)-Xnew(inds,:)) ; this is the p=2 matrix norm
  rmsd(i)=norm(Xref(inds,:)-Xnew(inds,:), 'fro') ; % Frobenius norm is consistent with VMD and other codes
% compute rotated coordinates
% put aligned structure into the same array to save memory
  xall(:,i)=Xnew(:,1).*oswgt + comref(1);
  yall(:,i)=Xnew(:,2).*oswgt + comref(2);
  zall(:,i)=Xnew(:,3).*oswgt + comref(3);
 end
%
end
