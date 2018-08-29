function [seq,resid]=pdbseq(pdb);

% aux functions to translate aa code
aas={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HIS' 'HSD' 'HSE' 'GLY' 'TRP' 'CYS' ;
      'A'   'R'   'D'   'Q'   'L'   'T'   'E'   'I'   'F'   'K'   'S'   'V'   'M'   'N'   'P'   'Y'   'H'   'H'   'H'   'G'   'W'   'C' } ;
aa1= @(x) char(aas(find(ismember(aas(:), x))+1)); % get 1 letter code from 3 letter code
aa3= @(x) char(aas(find(ismember(aas(:), x))-1)); % get 3 letter code from 1 letter code

rnum=[pdb.resSeq];
rname={pdb.resName};
%
natom=length(rnum);
seq='';

if natom<1
return
end

ires=1;
seq(ires)=aa1(rname(1));
resnum(ires)=rnum(1);
rprev=rnum(1);

for i=2:natom
 rthis=rnum(i);
 if(rthis~=rprev)
  ires=ires+1;
  seq(ires)=aa1(rname(i));
  resnum(ires)=rnum(i);
  rprev=rthis;
 end % if
end % for

if (nargout>1)
 resid=resnum;
end

end % function
