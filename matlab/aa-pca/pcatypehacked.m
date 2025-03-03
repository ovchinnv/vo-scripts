% add projections by HA type
cs=['rgbmcy'];
ts=['osv*^'];

inds={} ;
for hatype = {'H1', 'H2', 'H5', 'H6', 'H8', 'H9', 'H11', 'H12', 'H13', 'H16', 'H17', 'H18', 'H3', 'H4', 'H7', 'H10', 'H14', 'H15'};
 inds=[ inds ; {find(strcmp({strains.hatype},char(hatype)))}] ;
 leg=[ leg hatype ];
end

xoff = 3;
yoff = 3;
zoff = 3;

i=1;
for i=1:size(inds,1);
 inds2=cell2mat(inds(i));
 pcc=msa2proj(msamat(inds2,:), v, vind, qgrantham);
 ms=10;
 mc=cs(mod(i-1,numel(cs))+1);
 mt=ts(mod(i-1,numel(ts))+1);
% compute average coords : 
 pcav=mean(pcc,2);
 if (q2d)
  scatter(pcc(1,:), pcc(2,:), ms, mc, mt) ; view(2);
  text(pcav(1)+xoff, pcav(2)+yoff, char(leg(i+1)),'fontsize',16);
 else
  scatter3(pcc(1,:), pcc(2,:), pcc(3,:), ms, mc, mt) ; view(3);
  text(pcav(1)+xoff, pcav(2)+yoff, pcav(3)+zoff, char(leg(i+ileg)),'fontsize',16);
 end
% also compute cluster center coordinates
 [~,dim]=max(size(pcav));
 d1=sum(bsxfun(@minus,pcc,pcav).^2,dim); % distances from average
 [~,ind]=min(d1) ;% index of the strain that is closest to this cluster
 ccenter(i)=inds2(ind);
%
end
legend(leg);

return
% (re) compute projections for the center strains to make sure we've got the right ones :
pcc=msa2proj(msamat(ccenter,:), v, vind, qgrantham);
ms=20;
mc='k';
mt='*';
if (q2d)
  scatter(pcc(1,:), pcc(2,:), ms, mc, mt) ; view(2);
else
  scatter3(pcc(1,:), pcc(2,:), pcc(3,:), ms, mc, mt) ; view(3);
end

pdist=gdist(pcc');
pcatree=seqlinkage(pdist,'average',leg(end-length(inds)+1:end))
ttype='radial' ;
ttype='equaldaylight' ;
%f=figure(2);
%plot(pcatree,'Type',ttype)

