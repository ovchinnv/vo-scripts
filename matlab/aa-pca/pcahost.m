% add projections by host
cs=['rgbmcy'];
ts=['osv*^'];

inds={} ;
for host = unique({strains.host})
 inds=[ inds ; {find(strcmp({strains.host},char(host)))}] ;
 leg=[ leg host ];
end

xoff = 3;
yoff = 3;
zoff = 3;

i=1;
for i=1:size(inds,1);
 pcc=msa2proj(msamat(cell2mat(inds(i)),:), v, vind, qgrantham);
 ms=10 ;
 mc=cs(mod(i-1,numel(cs))+1);
 mt=ts(mod(i-1,numel(ts))+1);
% compute average coords : 
 pcav=mean(pcc,2);
 if (q2d)
  scatter(pcc(1,:), pcc(2,:), ms, mc, mt) ; view(2);
  text(pcav(1)+xoff, pcav(2)+yoff, char(leg(i+1)),'fontsize',16)
 else
  scatter3(pcc(1,:), pcc(2,:), pcc(3,:), ms, mc, mt) ; view(3);
  text(pcav(1)+xoff, pcav(2)+yoff, pcav(3)+zoff, char(leg(i+1)),'fontsize',16)
 end
end
legend(leg);
