function pcc=msa2proj(aln, v, inds, qgranth);
 if (nargin<4)
  coor=aln2coor(aln);
 else
  coor=aln2coor(aln,qgranth);
 end
%
 u=v(:,inds); % extract evs
 ncor=size(u,2);
% normalize ev
 for i=1:ncor
  u(:,i)=u(:,i)/norm(u(:,i));
 end
%
 nseq=size(aln,1);
 pcc=zeros(ncor,nseq);
% compute projection for each sequence
 for j=1:nseq
  for i=1:ncor
   pcc(i,j)=coor(j,:)*u(:,i) ;
  end
 end
