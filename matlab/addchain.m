% add chain to pdn using matlab
% note : must remove invalid fields ; in this case they are REMARK xxx fields because remarks 
% must be numbered

emr=pdbread('emr-scwrl.pdb');

pdb=emr.Model.Atom;
natom=length(pdb);

for i=1:natom
 if (pdb(i).segID=='EMR1')
  pdb(i).chainID='A';
 elseif (pdb(i).segID=='EMR2')
  pdb(i).chainID='B';
 end
end

% write pdb file
emr.Model.Atom=pdb;
%
pdbwrite('emr-scwrl-chain.pdb', emr);


