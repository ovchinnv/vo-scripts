function tmpmol=solvate(xsolu,ysolu,zsolu,solvent_pdb,buf,xovl,yovl,zovl,segname,overlap,density)
%
% NOTE :xyz-ovl are optional arrays with which we compute overlap
%
if ~exist('quiet') ; quiet=0 ; end
%return
nsolu=numel(xsolu); assert(nsolu==numel(ysolu)); assert(nsolu==numel(zsolu));
if (exist('xovl') | exist('yovl') | exist('zovl'))
 assert(exist('xovl')>0); assert(exist('yovl')>0); assert(exist('zovl')>0);
 novl=numel(xovl);
 assert(novl==numel(yovl));
 assert(novl==numel(zovl));
 qovl=1;
else
 qovl=0;
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~exist('buf') buf=9 ; end
if ~exist('overlap') overlap=2.75; end
if ~exist('density') density=1 ; end % g/mL
if ~exist('solvent_pdb') solvent_pdb='WAT.pdb' ; end % g/mL
if ~exist('segname') segname='W' ; end

Navogadro=6.022e23;
mw_h20=18.02;
qsavemem=1; % might also be faster, actually

xmin=min(xsolu)-buf ; xmax=max(xsolu)+buf;
ymin=min(ysolu)-buf ; ymax=max(ysolu)+buf;
zmin=min(zsolu)-buf ; zmax=max(zsolu)+buf;

if (~quiet)
 fprintf('==> Solute contains %d atoms\n',nsolu)
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
 fprintf('==> Number of water molecules in solvation box is: %d x %d x %d = %d\n', nx,ny,nz,nx*ny*nz)
end
if (~quiet); fprintf('==> Deleting molecules %12.5f Ang or farther from solute, or %12.5f Ang or closer to solute\n',buf,overlap) ; end

if (1) % now, only have low memory version
 buf2=buf.^2;
 overlap2=overlap.^2;
 ioks=zeros(nx*ny*nz,3) ;
 iok=0;
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
    if ( mdzz>buf2 | mdzz < overlap2) ; continue ; end
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
 ioks=ioks(1:iok,:); % truncate to valid entries
 XWAT=xwat(ioks(:,1));
 YWAT=ywat(ioks(:,2));
 ZWAT=zwat(ioks(:,3));
end
if (~quiet) ; 
 fprintf('==> There are #%d molecules remaining\n',numel(XWAT)); 
 fprintf('==> Writing solvent PDB file %s\n',solvent_pdb);
end
%scatter3(XWAT, YWAT, ZWAT) ;
tmpmol=seq2watpdb(XWAT,YWAT,ZWAT,'',segname,solvent_pdb,1); % best to use short segment names so that they can be appended to
