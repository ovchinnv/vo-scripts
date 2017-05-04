function rmsd=calc_rmsd(xall,yall,zall,xref,yref,zref,wgt,inds)
% based on bestfit, except with no alignment
%
 if (~exist('inds')) ; inds=[1:length(xref)] ; 
 else
  disp(['Computing RMSD using a subset of atoms : ',num2str(inds)]);
 end
%
 if (~exist('wgt') || isempty(wgt)) ; wgt=ones(size(xref)); end
 wgt=wgt/sum(wgt(inds));
 swgt=sqrt(wgt);
 oswgt=1./swgt(:);
% mass weighting
 Xref=[ swgt(inds).*( xref(inds) )  swgt(inds).*( yref(inds) )  swgt(inds).*( zref(inds) ) ];
%
 nall=(size(xall,2)); %; number of frames
 rmsd=zeros(1,nall);
 for i=1:nall
  if (~mod(i,100))
   disp(['==> Computing RMSD for set ',num2str(i)]);
  end
%
  X=[ swgt(inds).*(xall(inds,i)') swgt(inds).*(yall(inds,i)') swgt(inds).*(zall(inds,i)') ];
% compute RMSD :
  rmsd(i)=norm(Xref(:)-X(:),'fro'); % Frobenius norm is consistent with VMD and other codes
 end
%
end
