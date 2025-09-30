function tmpmol=solvate(xsolu,ysolu,zsolu,solvent_pdb,buf,xovl,yovl,zovl,segname,shape,overlap,density)
%
% NOTE :xyz-ovl are optional arrays with which we compute overlap
shell_=1; % solvate in a shell
ball_ =2; % solvate in a ball
%
if ~exist('quiet') ; quiet=0 ; end
%return
nsolu=numel(xsolu); assert(nsolu==numel(ysolu)); assert(nsolu==numel(zsolu));
if (exist('xovl') || exist('yovl') || exist('zovl'))
 assert(exist('xovl')>0); assert(exist('yovl')>0); assert(exist('zovl')>0);
 novl=numel(xovl);
 assert(novl==numel(yovl));
 assert(novl==numel(zovl));
 qovl=1;
else
 qovl=0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~exist('buf') buf=9 ; end
if ~exist('overlap') overlap=2.75; end
if ~exist('density') density=1 ; end % g/mL
if ~exist('solvent_pdb') solvent_pdb='WAT.pdb' ; end % g/mL
if ~exist('segname') segname='W' ; end
if ~exist('shape') shape=shell_ ; end

Navogadro=6.022e23;
mw_h20=18.02;
qsavemem=1; % might also be faster, actually

xmin=min(xsolu)-buf ; xmax=max(xsolu)+buf;
ymin=min(ysolu)-buf ; ymax=max(ysolu)+buf;
zmin=min(zsolu)-buf ; zmax=max(zsolu)+buf;

if (shape == ball_)
 xcen=0.5*(xmin+xmax); ycen=0.5*(ymin+ymax); zcen=0.5*(zmin+zmax);
 rball=0.5*max([abs(xmax-xmin),abs(ymax-ymin),abs(zmax-zmin)])+buf; % ball radius
 r2ball=rball^2;
 xmin=xcen-rball ; xmax=xcen+rball ;
 ymin=ycen-rball ; ymax=ycen+rball ;
 zmin=zcen-rball ; zmax=zcen+rball ;
end

if (~quiet)
 fprintf('==> Solute contains %d atoms\n',nsolu)
 fprintf('==> Solute buffer is %12.5f Angstroms\n',buf)
 fprintf('==> Overlap set contains %d atoms\n',novl)
 fprintf('==> Solvation box dimensions are : (%12.5f - %12.5f)x(%12.5f - %12.5f)x( %12.5f - %12.5f )\n',xmin,xmax,ymin,ymax,zmin,zmax)
end
% compute average distance between water molecules
molecules_per_angstrom3 = density * ( 1e6 * 1e-30 ) / mw_h20 * Navogadro ;
angstrom_per_molecule = 1./ molecules_per_angstrom3^(1./3) ; %spacing
dwat=angstrom_per_molecule;

% create coords :
xwat=[xmin:dwat:xmax+0.5*dwat] ;
ywat=[ymin:dwat:ymax+0.5*dwat] ;
zwat=[zmin:dwat:zmax+0.5*dwat] ;
%
nx=numel(xwat);ny=numel(ywat);nz=numel(zwat);
%
if (~quiet)
 fprintf('==> Number of water molecules in the initial solvation box is: %d x %d x %d = %d\n', nx,ny,nz,nx*ny*nz)
end
if (~quiet); 
 if (shape==shell_) 
  fprintf('==> Deleting molecules %12.5f Ang or farther from solute, or %12.5f Ang or closer to solute\n',buf,overlap);
 end
 if (shape==ball_) 
  fprintf('==> Deleting molecules %12.5f Ang or farther from solute center, or %12.5f Ang or closer to solute\n',rball,overlap);
 end
end

overlap2=overlap.^2;
ioks=zeros(nx*ny*nz,3) ;
iok=0;

if (shape==ball_)
 for i=1:nx
  if (~quiet); fprintf('==> Checking molecule #%d of %d\n',1+(i-1)*ny*nz, nx*ny*nz); end
  dxx = abs(xwat(i) - xcen) ;
%  if ( dxx>rball ) ; continue ; end ;% this is basically impossible b/c the box is based on sphere radius
  dxx=dxx.^2;
  for j=1:ny
   dyy = abs(ywat(j) - ycen) ;
%   if ( dyy>rball ) ; continue ; end
   dyy = dxx + dyy.^2;
   if ( dyy>r2ball ) ; continue ; end
   for k=1:nz
    dzz = abs(zwat(k) - zcen) ;
%    if ( dzz>rball ) ; continue ; end
    dzz = dyy + dzz.^2;
    if ( dzz>r2ball ) ; continue ; end
% if we are still here, it is still possible that the molecule overlaps with the additional overlap set, checked next
    if (qovl) % this is only executed if the overlap set is present and we haven't been kicked out.
%     if ( min((xwat(i) - xovl).^2 + (ywat(j) - yovl).^2 + (zwat(k) - zovl).^2) < overlap2 ) ; continue ; end
     if ( any( (xwat(i) - xovl).^2 + (ywat(j) - yovl).^2 + (zwat(k) - zovl).^2 < overlap2 )) ; continue ; end
    end
    iok=iok+1;
    ioks(iok,1)=i ;    ioks(iok,2)=j ;    ioks(iok,3)=k ;
   end %k
  end %j
 end %i
elseif (shape==shell_)
 buf2=buf.^2;
 for i=1:nx
  if (~quiet); fprintf('==> Checking molecule #%d of %d\n',1+(i-1)*ny*nz, nx*ny*nz); end
  dxx = xwat(i) - xsolu ;
  mdxx=min(abs(dxx));
  if ( mdxx>buf ) ; continue ; end
  dxx=dxx.^2;
  for j=1:ny
   dyy = ywat(j) - ysolu ;
   mdyy=min(abs(dyy));
   if ( mdyy>buf ) ; continue ; end
   dyy = dxx + dyy.^2;
   mdyy=min(dyy);
   if ( mdyy>buf2 ) ; continue ; end
   for k=1:nz
    dzz = zwat(k) - zsolu ;
    mdzz=min(abs(dzz));
    if ( mdzz>buf2 ) ; continue ; end
    dzz = dyy + dzz.^2;
    mdzz=min(dzz);
    if ( mdzz>buf2 || mdzz < overlap2) ; continue ; end
% if we are still here, it is still possible that the molecule overlaps with the additional overlap set, checked next
    if (qovl) % this is only executed if the overlap set is present and we haven't been kicked out.
%     if ( min((xwat(i) - xovl).^2 + (ywat(j) - yovl).^2 + (zwat(k) - zovl).^2) < overlap2 ) ; continue ; end
     if ( any( (xwat(i) - xovl).^2 + (ywat(j) - yovl).^2 + (zwat(k) - zovl).^2 < overlap2 )) ; continue ; end
    end
    iok=iok+1;
    ioks(iok,1)=i ;    ioks(iok,2)=j ;    ioks(iok,3)=k ;
   end %k
  end %j
 end %i
end
ioks=ioks(1:iok,:); % truncate to valid entries
XWAT=xwat(ioks(:,1));
YWAT=ywat(ioks(:,2));
ZWAT=zwat(ioks(:,3));


if (~quiet) ;
 fprintf('==> There are %d molecules remaining\n',numel(XWAT)); 
 fprintf('==> Writing solvent PDB file %s\n',solvent_pdb);
end
%scatter3(XWAT, YWAT, ZWAT) ;
tmpmol=seq2watpdb(XWAT,YWAT,ZWAT,'',segname,solvent_pdb,1); % best to use short segment names so that they can be appended to
% NOTE : seq2watpdb runs very slowly in octave
