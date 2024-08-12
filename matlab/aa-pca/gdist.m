% compute distance between all sequence pairs, or to a reference point
function d=gdist(coor,refcoor)
if (nargin>1) % disdtance to reference coordinate set
 d=sqrt(sum(bsxfun(@minus,coor,refcoor).^2,2)) ;
else
 nseq=size(coor,1);
 d=zeros(1,nseq*(nseq-1)/2);
 ind=1;
 for i=1:nseq-1
  di=nseq-i ;
  d(ind:ind+di-1)=sqrt(sum(bsxfun(@minus,coor(i+1:end,:),coor(i,:)).^2,2));
  ind=ind+di;
 end
end
