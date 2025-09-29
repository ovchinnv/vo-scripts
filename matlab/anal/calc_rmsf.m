function RMSF=calc_rmsf(x,y,z);
% compute root-mean-square flluctuations
 [natom,nframes]=size(x);
 X=[x ; y ; z; ]' ;% concatenate coordinates ; the variables must be in different columns; the observations in different rows
% variances
 VX=var(X,0);
% combine variances
 VR=zeros(1,natom);
 for i=1:3 % over components
  VR=VR+VX( (i-1)*natom + 1 : i*natom );
 end
 RMSF=sqrt(VR);
end
