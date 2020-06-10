
%pdbin='../HC.pdb'
%pdbout='HC'
% mutations to make (for now, a single file)
%mresid=[ 71 ];
%minsert={ ' ' } ;
%mrname={ 'ALA' };

%%%%%%%%%%%%%%%%%%%%%%%
aas={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HIS' 'HSD' 'HSE' 'GLY' 'TRP' 'CYS' ;
      'A'   'R'   'D'   'Q'   'L'   'T'   'E'   'I'   'F'   'K'   'S'   'V'   'M'   'N'   'P'   'Y'   'H'   'H'   'H'   'G'   'W'   'C' } ;
aa1= @(x) char(aas(find(ismember(aas(:), x))+1)); % get 1 letter code from 3 letter code
aa3= @(x) char(aas(find(ismember(aas(:), x))-1)); % get 3 letter code from 1 letter code


mol=pdbread(pdbin);
pdb=mol.Model.Atom;
resid=[pdb(:).resSeq];
insertion=char({pdb.iCode});
if isempty(insertion)
 insertion(:,1)=' '; % compat
end
rname={pdb(:).resName};

mstr='';

for i=1:length(mresid)
 id=mresid(i);
 ins=char(minsert(i));
 if isempty(ins)
  ins(:,1)=' '; % compat
 end
 name=char(mrname(i));
 minds=find( resid(:)==id & insertion==ins ) ;
 if (~isempty(minds))
  oldname=rname(minds(1));
  for ind=minds'
   pdb(ind).resName=name;
  end
  mstr=[mstr,'_',aa1(oldname),num2str(id),strtrim(ins),aa1(name)];
 end
end

if (~isempty(mstr))
 if (~exist('pdbout'))
  pdbout=pdbin;
 end
 ind=strfind(pdbout,'.');
 if(~isempty(ind))
  newfile=[pdbout(1:ind(end)-1),mstr,'.pdb'];
 else
  newfile=[pdbout,mstr,'.pdb'];
 end 
 mol.Model.Atom=pdb;
 pdbwrite(newfile,mol);
end

