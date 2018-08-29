% truncate PDB to limits provided in another pdb using alignment
addpath('~/scripts/matlab/build')

if (~exist('read'))
 read=1;
end

if (read)

pdbfile0='heavy.pdb';
pdbfile1='HC.pdb';
pdboutfile='heavy_trunc.pdb';

disp(['==>Reading pdb structure from file ',pdbfile0,' ...']);
mol0=pdbread(pdbfile0);
pdb0=mol0.Model.Atom;
[seq0,resnum0]=pdbseq(pdb0);

disp(['==>Reading pdb structure from file ',pdbfile1,' ...']);
mol1=pdbread(pdbfile1);
pdb1=mol1.Model.Atom;
[seq1,resnum1]=pdbseq(pdb1);

read=0;

end

[score,al]=swalign(seq0,seq1);
%showalignment(al)
% in case alignment is truncated, as when using swalign, find beginning residue indices in the alignment
seq=al(1,:); seq=strrep(seq,'-',''); ires(1)=strfind(seq0,seq)-1;
seq=al(3,:); seq=strrep(seq,'-',''); ires(2)=strfind(seq1,seq)-1;

nres=size(al,2);
% align resids with sequence alignment :
for i=1:nres
% first sequence
 if (al(1,i)=='-')
  alresn(1,i)=0;
 else
  ires(1)=ires(1)+1;
  alresn(1,i)=resnum0(ires(1));
 end
% second
 if (al(3,i)=='-')
  alresn(2,i)=0;
 else
  ires(2)=ires(2)+1;
  alresn(2,i)=resnum1(ires(2));
 end
end
%
% truncate pdb to desired set
ibeg=find( alresn(2,:)==max(1,min(alresn(2,:))),1);
while ((alresn(1,ibeg)==0) & (ibeg<nres))
 ibeg=ibeg+1;
end

iend=find( alresn(2,:)==min(111,max(alresn(2,:))),1);
while ((alresn(1,iend)==0) & (iend>0))
 iend=iend-1;
end

% find atom limits:
resid0=[pdb0.resSeq];
icut=find(resid0==alresn(1,ibeg),1,'first');
ecut=find(resid0==alresn(1,iend),1,'last');

pdbout=pdb0(icut:ecut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdbout=fixcharmm(pdbout);
mol2.Model.Atom=pdbout;
pdbwrite(pdboutfile, mol2);
system(['echo END >> ',pdboutfile]);
