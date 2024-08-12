% create solvation file with water oxygens
% iterative algorithm to minimize using large solute segments
addpath('~/scripts/matlab/build')

pdbfile='np60_m.pdb';

if (~exist('qpdb')) ; qpdb=0 ; end
if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
% molecule=pdbread(pdbfile);
 molecule=readpdb(pdbfile,1); % octave/matlab compatible ; 2nd arg for verbosity
 qpdb=1;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buf=9;

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
element=cellfun(@(s) s(1), aname); % from struc
natom=length(pdb);

% loop over multiple selections, solvating each selection separately; checking must be done on the entire protein; each sucessive water selection cannot ovelap w/ prior
%icheck=(~cellfun('isempty',strfind(segid,'NP'))) ; % all atoms in a segment
icheck= [ (~cellfun('isempty',strfind(segid,'NP'))) & (element~='H') ] ; % exclude hydrogens from checking for a little extra speed
xchk=xpdb(icheck);ychk=ypdb(icheck);zchk=zpdb(icheck);

isel=icheck ;
xsolu=xpdb(isel);
ysolu=ypdb(isel);
zsolu=zpdb(isel);
solvent_pdb=['WAT.pdb'];
watmol=solvate(xsolu,ysolu,zsolu,solvent_pdb,buf,xchk,ychk,zchk,['W'],2); % assume that a single pdb will be written
watpdb=watmol.Model.Atom ; % grab waters out of the solvent segment to use for checking in the next steps :
nallwat=length(watpdb) ;
xwat=[watpdb.X]';
ywat=[watpdb.Y]';
zwat=[watpdb.Z]';
wnam={watpdb.AtomName}';
