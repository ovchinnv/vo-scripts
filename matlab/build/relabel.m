% restore correct residue numbering
% change the numbering of pdb1 to that of pdb0
%


close all;
addpath('~/scripts/matlab/build');

if (~exist('read'))
 read=1;
end

if (read)

name='bg505sosip_3bnc60';

pdbfile0='G120-ALL.pdb';
pdbfile1=['../',name,'.pdb'];
pdboutfile=[name,'.pdb'];

disp(['==>Reading pdb structure from file ',pdbfile0,' ...']);
mol0=pdbread(pdbfile0);
pdb0=mol0.Model.Atom;

anum0     = [pdb0.AtomSerNo]';
chainid0  = [pdb0.chainID];
rname0    = {pdb0.resName}' ;
resid0    = [pdb0.resSeq]' ;
%insertion0= char({pdb0.iCode});
insertion0= {pdb0.iCode}';
segid0    = {pdb0.segID}' ;
natom0    = length(pdb0);

[seq0,resnum0,resnumi0]=pdbseqi(pdb0);
resnumi0(find(cellfun('isempty',resnumi0)))={' '}; % replace empty strings with spaces

disp(['==>Reading pdb structure from file ',pdbfile1,' ...']);
mol1=pdbread(pdbfile1);
pdb1=mol1.Model.Atom;
% make sure only to consider chain A

chainid1  = [pdb1.chainID];
inda=(find(chainid1=='A'));
pdb1=pdb1(inda);

anum1     = [pdb1.AtomSerNo]';
chainid1  = [pdb1.chainID];
rname1    = {pdb1.resName}' ;
resid1    = [pdb1.resSeq]' ;
%insertion1= char({pdb1.iCode});
insertion1= {pdb1.iCode}';
segid1    = {pdb1.segID}' ;
natom1    = length(pdb1);
insertion1(find(cellfun('isempty',insertion1)))={' '}; % replace empty strings with spaces

[seq1,resnum1,resnumi1]=pdbseqi(pdb1);
resnumi1(find(cellfun('isempty',resnumi1)))={' '}; % replace empty strings with spaces

read=0;

end

%[score,al]=nwalign(seq0,seq1); % this is a global aligment, but in this case it includes missing residues !
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
  alresni(1,i)=' ';
 else
  ires(1)=ires(1)+1;
  alresn(1,i)=resnum0(ires(1));
  alresni(1,i)=resnumi0(ires(1));
 end
% second
 if (al(3,i)=='-')
  alresn(2,i)=0;
  alresni(2,i)=' ';
 else
  ires(2)=ires(2)+1;
  alresn(2,i)=resnum1(ires(2));
  alresni(2,i)=resnumi1(ires(2));
 end
end

% change residue numbering
pdbout=pdb1;

ialn=1; % index into alignment
i=1; % index into pdb to change

while ( (i<=natom1) & (ialn<=nres) )
 if ( resid1(i) < alresn(2,ialn))
  i=i+1; % skip atom
  continue
 elseif ( resid1(i) > alresn(2,ialn))
  ialn=ialn+1;
  continue
 else % they must be equal -- proceed to compare insertion codes
  if ( char(insertion1(i)) < char(alresni(2,ialn)))
   i=i+1; % skip atom
   continue
  elseif ( char(insertion1(i)) > char(alresni(2,ialn)))
   ialn=ialn+1;
   continue
  else % they must be equal, so we found the correct residue  -- replace
   pdbout(i).resSeq=alresn(1,ialn);
   pdbout(i).iCode=char(alresni(1,ialn));
   i=i+1; % next atom
  end
 end
end

pdbout=fixcharmm(pdbout);
mol2.Model.Atom=pdbout;
pdbwrite(pdboutfile, mol2);
system(['echo END >> ',pdboutfile]);
