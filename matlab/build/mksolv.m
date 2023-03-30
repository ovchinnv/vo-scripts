% create solvation file with water oxygens
addpath('~/scripts/matlab/build')

if (~exist('pdbfile')) ; pdbfile='spy.pdb'; end

if (~exist('qpdb')) ; qpdb=0 ; end
if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 molecule=pdbread(pdbfile);
 qpdb=1; %
end
%
if ~exist('quiet') ; quiet=0 ; end
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~exist('buf') buf=9 ; end
if ~exist('overlap') overlap=2.5; end
if ~exist('selection')
 selection='[1:natom]'' & cellfun(''isempty'',strfind(segid,''WAT''))' ;
end
if ~exist('density') density=1 ; end % g/mL
if ~exist('solvent_pdb') solvent_pdb='WAT.pdb' ; end % g/mL

Navogadro=6.022e23;
mw_h20=18.02;
qsavemem=1; % might also be faster, actually

pdb=molecule.Model.Atom;

anum= [pdb.AtomSerNo]';
aname={pdb.AtomName}' ;
rname={pdb.resName}' ;
segid={pdb.segID}' ;
resid=[pdb.resSeq]' ;
insertion=char({pdb.iCode});
occu =[pdb.occupancy]' ;
bet  =[pdb.tempFactor]' ;
xpdb =[pdb.X]';
ypdb =[pdb.Y]';
zpdb =[pdb.Z]';
chainid=[pdb.chainID];
element=[pdb.element];
natom=length(pdb);

iselection=eval(selection) ;

xsolu=xpdb(iselection);
ysolu=ypdb(iselection);
zsolu=zpdb(iselection);
nsolu=numel(xsolu);

xmin=min(xsolu)-buf ; xmax=max(xsolu)+buf;
ymin=min(ysolu)-buf ; ymax=max(ysolu)+buf;
zmin=min(zsolu)-buf ; zmax=max(zsolu)+buf;

if (~quiet)
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

if (qsavemem)
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
    iok=iok+1;
    ioks(iok,1)=i ;    ioks(iok,2)=j ;    ioks(iok,3)=k ;
   end %k
  end %j
 end %i
 ioks=ioks(1:iok,:); % truncate to valid entries
 XWAT=xwat(ioks(:,1));
 YWAT=ywat(ioks(:,2));
 ZWAT=zwat(ioks(:,3));
else
 [XWAT,YWAT,ZWAT]=meshgrid(xwat,ywat,zwat); % take these to be the centers of water molecules
% flatten
 XWAT=XWAT(:);
 YWAT=YWAT(:);
 ZWAT=ZWAT(:);
 NWAT=numel(XWAT);
% compute distance to solvent
 dx=repmat(xsolu',NWAT,1)-repmat(XWAT,1,nsolu);
 dy=repmat(ysolu',NWAT,1)-repmat(YWAT,1,nsolu);
 dz=repmat(zsolu',NWAT,1)-repmat(ZWAT,1,nsolu);
 d2=dx.^2 + dy.^2 + dz.^2 ;
 d2min=min(d2,[],2) ; % the smallest squared distance to any atom in selection

 iclose=find(d2min<buf^2) ; % solvent inside shell
 ishell=find(d2min(iclose)>overlap^2) ; % too far for a clash
%
 iok=iclose(ishell) ;
%
 XWAT=XWAT(iok);
 YWAT=YWAT(iok);
 ZWAT=ZWAT(iok);
end
if (~quiet) ; 
 fprintf('==> There are #%d molecules remaining\n',numel(XWAT)); 
 fprintf('==> Writing solvent PDB file %s\n',solvent_pdb);
end
%scatter3(XWAT, YWAT, ZWAT) ;
tmpmol=seq2watpdb(XWAT,YWAT,ZWAT,'','W',solvent_pdb,1); % best to use short segment names so that they can be appended to
