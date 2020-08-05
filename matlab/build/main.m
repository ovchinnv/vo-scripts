% MATLAB routines for trajectory analysis
%
addpath('~/scripts/matlab/build')
%
name='1eo8';

% segment ids for the protein chains : (optional)
chain_segid={ 'H', 'L', 'A', 'B' ; ...
              'HC  ', 'LC  ','HA1 ', 'HA2 ' };

ch2seg= @(x) char(chain_segid(1+find(ismember(chain_segid(:), x)))); % get segment ID from chainID
ch2segt= @(x) strtrim(char(chain_segid(1+find(ismember(chain_segid(:), x))))); % get segment ID from chainID; trimmed to removed leading/trailing spaces

pdbfile=[name,'.pdb']; % in lieu of a "structure" file

if (~exist('qpdb')) ; qpdb=0 ; end

if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 molecule=pdbread(pdbfile);
 qpdb=1; %
end

pdb=molecule.Model.Atom;

% manipulate PDB format and write out pdb files for charmm
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
%
% mark all chains (replacement of explicit code)
for ch = chain_segid(1,:)
 cch=char(ch)
 mcmd=['chain',cch,'=ismember(chainid,''',cch,'''); segid=ch2seg(''',cch,''');for i=find(chain',cch,'); pdb(i).segID=segid; end']
 eval(mcmd);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change charmm terminology
pdb=fixcharmm(pdb);
% create new model
mol2.Model.Atom=pdb;
% write out pdb files
if (1)
pdbout(mol2,[ch2segt('H'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainH));
pdbout(mol2,[ch2segt('L'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainL));
pdbout(mol2,[ch2segt('B'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainB));
pdbout(mol2,[ch2segt('A'),'.pdb'],xpdb,ypdb,zpdb,[],[],find(chainA));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% now process heterogeneous atom entries %%%%%%%%%%%%
if (isfield(molecule.Model,'HeterogenAtom'))
 pdbh=molecule.Model.HeterogenAtom;
 anum= [pdbh.AtomSerNo]';
 aname={pdbh.AtomName}' ;
 rname={pdbh.resName}' ;
 segid={pdbh.segID}' ;
 resid=[pdbh.resSeq]' ;
 insertion=char({pdbh.iCode});
 occu =[pdbh.occupancy]' ;
 bet  =[pdbh.tempFactor]' ;
 xpdb =[pdbh.X]';
 ypdb =[pdbh.Y]';
 zpdb =[pdbh.Z]';
 chainid=[pdbh.chainID];
 element=[pdbh.element];
%
% water
 water=ismember(rname,'HOH');
 resnum=1; % renumber residues
 for i=find(water)'
  pdbh(i).resName='TP3';
  pdbh(i).AtomNameStruct.chemSymbol='O';
  pdbh(i).AtomNameStruct.remoteInd='H';
  pdbh(i).AtomNameStruct.branch='2';
  pdbh(i).segID='XWAT';
  pdbh(i).element='';
  pdbh(i).resSeq=resnum;%sprintf('%5d',resnum);
  resnum=resnum+1; % renumber residues
 end

 mol3.Model.Atom=pdbh;
% write out pdb files
 pdbout(mol3,'xwat.pdb',xpdb,ypdb,zpdb,[],[],find(water));

% rename NAG atoms to match CHARMM carbohydrate force field
%
 nag=ismember(rname,'NAG');
 nagc8=nag & ismember(aname,'C8') ;
 nago7=nag & ismember(aname,'O7') ;
 nagc7=nag & ismember(aname,'C7') ;
 nagn2=nag & ismember(aname,'N2') ;

 for i=find(nagc8)'
  pdbh(i).AtomNameStruct.chemSymbol='C';
  pdbh(i).AtomNameStruct.remoteInd='T';
 end

 for i=find(nago7)'
  pdbh(i).AtomNameStruct.chemSymbol='O';
  pdbh(i).AtomNameStruct.remoteInd='';
 end

 for i=find(nagc7)'
  pdbh(i).AtomNameStruct.chemSymbol='C';
  pdbh(i).AtomNameStruct.remoteInd='';
 end

 for i=find(nagn2)'
  pdbh(i).AtomNameStruct.chemSymbol='N';
  pdbh(i).AtomNameStruct.remoteInd='';
 end

 mol3.Model.Atom=pdbh;

% process glycans
 glyco;
end
% create pdbs for missing loops
missing;
%
% note : the missing residues are all at the ends, so omit them
return
% write pdbs with missing coordinates included
% note: will only model one internal missing loop in HA1 and HA2
pdbs={'HA1.pdb', 'HA1-MISSING-2.pdb'}
molout=combine_pdbs(pdbs);
pdbwrite('HA1-ALL.pdb', molout); system('echo END >> HA1-ALL.pdb');
%
pdbs={'HA2.pdb', 'HA2-MISSING-2.pdb'}
molout=combine_pdbs(pdbs);
pdbwrite('HA2-ALL.pdb', molout); system('echo END >> HA2-ALL.pdb');
%
