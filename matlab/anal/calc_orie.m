function [xall, yall, zall, comall, pcall]=calc_orie(xall,yall,zall,wgt,inds)
% compute center of mass and principal components
%
 if (~exist('inds')) ; inds=[1:size(xall,1)] ; 
 else
  disp(['Computing radius of gyration (Rgyr) using a subset of atoms : ',num2str(inds')]);
 end
%
 if (~exist('wgt') || isempty(wgt)) ; wgt=ones(size(xref)); end
 wgt=wgt/sum(wgt(inds));
 swgt=sqrt(wgt);
%
 nall=(size(xall,2)); %; number of frames
 rgyr=zeros(1,nall);
 for i=1:nall
  if (~mod(i,100))
   disp(['==> Computing orientation to laboratory frame for set ',num2str(i)]);
  end
%
  com=[wgt(inds)'*xall(inds,i) wgt(inds)'*yall(inds,i) wgt(inds)'*zall(inds,i)];
%
  X=[ swgt(inds).*(xall(inds,i)-com(1)) swgt(inds).*(yall(inds,i)-com(2)) swgt(inds).*(zall(inds,i)-com(3)) ]; % COM subtracted
  X2=transpose(X)*X;  % not quite the inertia tensor
  X2=trace(X2)*eye(3) - X2; 
% diagonalize to compute principal components
  [pc,pv]=eig(X2);
% decide whether to move actual coordinates:
  Xnew=[ xall(:,i)-com(1) yall(:,i)-com(2) zall(:,i)-com(3) ] * pc ; % apply pc^-1 which is same as pc^T
  xall(:,i) = Xnew(:,1);
  yall(:,i) = Xnew(:,2);
  zall(:,i) = Xnew(:,3);
% save
  comall(1:3,i)=com ;
  pcall(1:3,1:3,i)=pc ;
 end
%
end
