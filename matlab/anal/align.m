% perform alignment using procrustes

% cannot use mass weighting
% not quite the same as VMD/CHARMM because weights cannot be passed into the procrastes routine
% thus, until Procrustes gets fixed, only a uniformly-weighted fit is possible
%
% nomass=1;
 if (~exist('qmass')) ; qmass=0; end;
 if (qmass) ; wgt=mass; else ; wgt=ones(size(mass)); end
%
 wgt=wgt/sum(wgt);
 swgt=sqrt(wgt);
 oswgt=1./swgt(:);

 if (~exist('xref'))
  xref=xpdb; yref=ypdb; zref=zpdb;
 end
% prepare reference coordinates
 comref=[wgt*xref(:) wgt*yref(:) wgt*zref(:)];
% mass weighting after translation to COM
 Xref=[ swgt(:).*( xref(:) - comref(1) )  swgt(:).*( yref(:) - comref(2) )  swgt(:).*( zref(:) - comref(3) ) ];
%
 rmsd=zeros(1,nall);
 for i=1:nall
  if (~mod(i,100))
   disp(['==> Aligning coordinate set ',num2str(i)]);
  end
%
  com=[wgt*xall(:,i) wgt*yall(:,i) wgt*zall(:,i)];
%
  X=[ swgt(:).*(xall(:,i)-com(1)) swgt(:).*(yall(:,i)-com(2)) swgt(:).*(zall(:,i)-com(3)) ]; % also subtract COM
  [d,Xnew,trans]=procrustes(Xref, X, 'scaling',0,'reflection',0);
% recompute coordinates without translation component
  Xnew=X*trans.T ; % pure rotation (only in this case recover correct magnitudes below)
% compute RMSD :
  rmsd(i)=norm(Xref(:)-Xnew(:));
% compute rotated coordinates
% put aligned structure into the same array to save memory
  xall(:,i)=Xnew(:,1).*oswgt + comref(1);
  yall(:,i)=Xnew(:,2).*oswgt + comref(2);
  zall(:,i)=Xnew(:,3).*oswgt + comref(3);
 end
%
