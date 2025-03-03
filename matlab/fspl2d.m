function fspl2d=pp2val(spl2d,kshow,zshow)
% here, we read directly the spline file, compute the forces and plot (i.e. without using any other data)
% should match results from showspl2d
% written as a function
ka=spl2d.breaks{1};
z0=spl2d.breaks{2};
nka=numel(ka);
nz0=numel(z0);
korder=spl2d.order{1};
zorder=spl2d.order{2};
%
nz=numel(zshow);
nk=numel(kshow);
f2d=zeros( nz, nk );
%
%if (qnative) % cannot use the matlab function below ! Must have the wrong forat for the spline file, and the help is unhelpful
%   f2d(iz,ik)=fnval(spl2d, { kshow(ik) zshow(iz) } );
%else
for i=1:nk
% go by brute force (see showspl2d):
 if     (kshow(i)<ka(1))     ; ik=1;
 elseif (kshow(i)>ka(nka-1)) ; ik=nka-1 ;
 else ; ik=find(ka<=kshow(i),1,'last') ;
 end
 kap=(kshow(i)-ka(ik)).^[korder-1:-1:0]; % curvature grid
 for j=1:nz
  if     (zshow(j)<z0(1))     ; iz=1 ;
  elseif (zshow(j)>z0(nz0-1)) ; iz=nz0-1 ;
  else ; iz=find(z0<=zshow(j),1,'last') ; % the last point that is smaller than z
  end
  zp=(zshow(j)-z0(iz)).^[zorder-1:-1:0]; % NOTE : computing "by hand"
  coefs=spl2d.coefs(:,:,ik,iz) ; % after permuting dimensions (see mkspl2d)
  fspl2d(j,i) = kap*coefs*zp' ; % quadratic
 end % over z(j)
end % over k(i)
