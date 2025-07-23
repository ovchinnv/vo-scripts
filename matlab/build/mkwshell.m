% create solvation file with water oxygens (water shell)
addpath('~/scripts/matlab/build')

if (~exist('pdbfile')); pdbfile='b12_m.pdb';end

if (~exist('qpdb')) ; qpdb=0 ; end
if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
% molecule=pdbread(pdbfile); % matlab native
 molecule=readpdb(pdbfile,1); % octave/matlab compatible ; 2nd arg for verbosity
 qpdb=1;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buf=8.5; % solvation buffer thickness

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
% could also use a single selection, which would be faster for some protein shapes
%icheck=(~cellfun('isempty',strfind(segid,'NP'))) ; % all atoms in any segment whose name contains the string "NP"
icheck= [ ~(cellfun('isempty',(strfind(segid,'HC'))) & cellfun('isempty',(strfind(segid,'LC'))) & cellfun('isempty',(strfind(segid,'S')))...
          & cellfun('isempty',(strfind(segid,'XWAT'))))...
          & (element~='H') ] ; % exclude hydrogens from checking for a little extra speed
xchk=xpdb(icheck);ychk=ypdb(icheck);zchk=zpdb(icheck);
% create selections
% (1) split into multiple :
%sels={ ismember(segid,'HC1') ; ismember(segid,'LC1') ; ismember(segid,'HC2') ; ismember(segid,'LC2') }; % could include sugars
% (2) one large selection :
sels={ ismember(segid,'HC1') | ismember(segid,'LC1') | ismember(segid,'HC2') | ismember(segid,'LC2') }; % could include sugars

nsels=numel(sels);
for isel=1:nsels
 iss=num2str(isel);
 sel=sels{isel};
 xsolu=xpdb(sel);
 ysolu=ypdb(sel);
 zsolu=zpdb(sel);
 solvent_pdb=['WAT',iss,'.pdb'];
 watmol=solvate(xsolu,ysolu,zsolu,solvent_pdb,buf,xchk,ychk,zchk,['W',iss]);
 watpdb=watmol.Model.Atom ; % grab waters out of the solvent segment to use for checking in the next steps :
 nallwat=length(watpdb) ;
 xwat=[watpdb.X]';
 ywat=[watpdb.Y]';
 zwat=[watpdb.Z]';
 wnam={watpdb.AtomName}';
%% iwat=[1:nallwat] ; % all atoms
 iwat=find(ismember(strtrim(wnam),'OH2')); % only water oxygens
 xchk=[xchk;xwat(iwat)]; ychk=[ychk;ywat(iwat)]; zchk=[zchk;zwat(iwat)]; % add new waters to overlap check for next iter
end
%solvate(xsolu,ysolu,zsolu,solvent_pdb,buf,xsolu,ysolu,zsolu)
