% restore correct sequence nomenclature in a model made by MODELLER
% using the original incomplete file as a template
% this version is unsophisticated : 
% we simply read the two PDBs side by side and transfer residue and segment names
% 8/18 minor update : relabel all models so that we can decide on the best one more systematically later


pdbfile0='HA1-ALL.pdb';
% note that it is unclear which model is the "best" one
pdbmodelname='HA1-ALL.B9999';
pdboutname='HA1-MODEL';
imodel=1;
emodel=10;

disp(['==>Reading pdb structure from file ',pdbfile0,' ...']);
mol0=pdbread(pdbfile0);
zeropad='00000000000' ; % padding string
pdb0=mol0.Model.Atom;


for i=imodel:emodel
%
pdbfile1=[pdbmodelname,zeropad(1:3-floor(log10(i))),num2str(i),'.pdb'];
pdboutfile=[pdboutname,num2str(i),'.pdb'];

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
  if pdbgt(pdb1(i1+1),pdb1(i1))
   % cycle forward in pdb0
   while ~pdbgt(pdb0(i0+1),pdb0(i0))
    i0=i0+1;
   end
   i0=i0+1;
  end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdbout=fixcharmm(pdbout);
mol2.Model.Atom=pdbout;
pdbwrite(pdboutfile, mol2);
system(['echo END >> ',pdboutfile]);

end
