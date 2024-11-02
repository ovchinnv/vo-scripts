function ppi=splint(pp)
% integrate a spline in pp form segment-by-segment
 ppi=pp ; % initialize to the integrand
 ppi.order=pp.order+1 ; % increase order ;
 ppi.coefs=bsxfun(@times,pp.coefs(:,1:pp.order),1./[pp.order:-1:1] ); % integrate each polynomial ;
 ppi.coefs(:,ppi.order)=0 ; % now the poly will compute as zero at all break points (fixed below)
 for i=2:ppi.pieces % compute leading term correctly
  ppi.coefs(i,ppi.order) =(ppi.breaks(i)-ppi.breaks(i-1)).^[ppi.order-1:-1:0]  * ppi.coefs(i-1,:)'; % match from prev segment
 end
