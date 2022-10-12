% restore correct sequence nomenclature in a model made by MODELLER
% using the original incomplete file as a template
% this version is unsophisticated : 
% we simply read the two PDBs side by side and transfer residue and segment names

addpath('~/scripts/matlab/build')

pdbfile0='ha_1eo8_07_fkwgfk.pdb';
% note that it is unclear which model is the "best" one
pdbfile1='ha_1eo8_07_fkwgfk.B99990001.pdb';
%pdbfile1='G120-ALL.B99990002.pdb';
pdboutfile='temp.pdb';

disp(['==>Reading pdb structure from file ',pdbfile0,' ...']);
mol0=pdbread(pdbfile0);
pdb0=mol0.Model.Atom;

disp(['==>Reading pdb structure from file ',pdbfile1,' ...']);
mol1=pdbread(pdbfile1);
pdb1=mol1.Model.Atom;

anum     = [pdb0.AtomSerNo]';
chainid  = [pdb0.chainID];
rname    = {pdb0.resName}' ;
resid    = [pdb0.resSeq]' ;
insertion= char({pdb0.iCode});
segid    = {pdb0.segID}' ;
natom    = length(pdb0);

pdbout=pdb1 ; % copy pdb structure

i0=1;
for i1=1:length(pdb1)
 pdbout(i1).chainID=char(chainid(i0));
 pdbout(i1).resSeq=resid(i0);
 if (~isempty(insertion)) ; pdbout(i1).iCode=char(insertion(i0)); end
 pdbout(i1).segID=char(segid(i0));
% check if the next atom is in a different residue :
 if (i1<length(pdb1))
  if ~pdbeq(pdb1(i1+1),pdb1(i1))
   % cycle forward in pdb0
   while ~( i0>=length(pdb0) | ~pdbeq(pdb0(min(i0+1,length(pdb0))),pdb0(i0)) )
    i0=min(i0+1, length(pdb0));
   end
   i0=min(i0+1, length(pdb0));
  end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdbout=fixcharmm(pdbout);
mol2.Model.Atom=pdbout;
pdbwrite(pdboutfile, mol2);
system(['echo END >> ',pdboutfile]);
