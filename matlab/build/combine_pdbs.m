function molout=combine_pdbs(pdbs,qoctpdb);
% combine two PDB files into one, taking into account possible
% isertions; useful for including missing loops
if (nargin<2)
 qoctpdb=exist('OCTAVE_VERSION','builtin');
end
qloud=1;
%
% read first pdb:
%
npdbs=length(pdbs);
if (npdbs<=1)
 warning(' Number of PDBs to combine should be greater than one');
% return; % can still proceed
end

pdbfile=char(pdbs(1));
disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
if (qoctpdb)
 mol1=readpdb(pdbfile,qloud);
else
 mol1=pdbread(pdbfile);
end
pdb1=mol1.Model.Atom;

% process remaining files

for i=2:npdbs
 pdbfile=char(pdbs(i));
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 if (qoctpdb)
  mol2=readpdb(pdbfile,qloud);
 else
  mol2=pdbread(pdbfile);
 end
 pdb2=mol2.Model.Atom;
% combine
% loop over records in pdb2, 
 natom2=length(pdb2);
%
 for j=1:natom2
% binary search limits
% find first entry in PDB1 such that j-record in PDB2 is larger
  i1=1;
  i2=length(pdb1);
  while i2-i1>1
   imid = fix((i1+i2)/2);
   if pdbgt(pdb1(imid),pdb2(j)) % pdb1(imid) greater than pdb2(j)
    i2=imid;
   else
    i1=imid;
   end % if
  end % while
% insert record
  if pdbgt(pdb1(i1),pdb2(j))
   pdb1(i1+1:end+1)=pdb1(i1:end);
   pdb1(i1)=pdb2(j);
  elseif pdbgt(pdb1(i2),pdb2(j))
   pdb1(i2+1:end+1)=pdb1(i2:end);
   pdb1(i2)=pdb2(j);
  else % insert after i2
   pdb1(i2+2:end+1)=pdb1(i2+1:end);
   pdb1(i2+1)=pdb2(j);
  end

 end % pdb2 atoms
end % npdbs
% write combined pdb
% need to renumber atom numbers in the pdb because output is sorted based on the serial numbers
for i=1:length(pdb1)
 pdb1(i).AtomSerNo=i;
end
molout.Model.Atom=pdb1;

