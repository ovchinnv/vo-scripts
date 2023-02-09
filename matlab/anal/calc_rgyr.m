function rgyr=calc_rgyr(xall,yall,zall,wgt,inds)
%
 if (~exist('inds')) ; inds=[1:size(xall,1)] ; 
 else
  disp(['Computing radius of gyration (Rgyr) using a subset of atoms : ',num2str(inds')]);
 end
%
 if (~exist('wgt') || isempty(wgt)) ; wgt=ones(size(xref)); end
 wgt=wgt/sum(wgt(inds));
 swgt=sqrt(wgt);
% oswgt=1./swgt(:);
%
 nall=(size(xall,2)); %; number of frames
 rgyr=zeros(1,nall);
 for i=1:nall
  if (~mod(i,100))
   disp(['==> Computing Rgyr for set ',num2str(i)]);
  end
%
  com=[wgt(inds)'*xall(inds,i) wgt(inds)'*yall(inds,i) wgt(inds)'*zall(inds,i)];
%
  X=[ swgt(inds).*(xall(inds,i)-com(1)) swgt(inds).*(yall(inds,i)-com(2)) swgt(inds).*(zall(inds,i)-com(3)) ]; % COM subtracted
% compute Rgyr :
  rgyr(i)=norm(X(:),'fro'); % Frobenius norm
 end
%
end
