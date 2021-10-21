function molout=concat_pdbs(pdbs);
% concatenate PDB files

% read first pdb:
npdbs=length(pdbs);
if (npdbs<=1)
 warning(' Number of PDBs to combine should be greater than one');
% return; % can still proceed
end

pdbfile=char(pdbs(1));
disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
mol1=pdbread(pdbfile);
pdb1=mol1.Model.Atom;

% process remaining files

for i=2:npdbs
 pdbfile=char(pdbs(i));
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 mol2=pdbread(pdbfile);
 pdb2=mol2.Model.Atom;
% combine
% loop over records in pdb2, 
 natom2=length(pdb2);
%
 i2=length(pdb1);
 pdb1(i2+1:i2+natom2)=pdb2;

end % npdbs
% write combined pdb
% need to renumber atom numbers in the pdb because output is sorted based on the serial numbers
for i=1:length(pdb1)
 pdb1(i).AtomSerNo=i;
end
molout.Model.Atom=pdb1;

