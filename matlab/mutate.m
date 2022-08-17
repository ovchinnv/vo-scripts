% example vars :
%pdbin='3lzg-5naz-fus|hux-candidate-3.pdb'
%pdbout='mut-ss';
%mresid=[ 310 93 ] ;
%minsert={ ' ' ' ' } ;
%mrname={ 'CYS' 'CYS' } ;
%msegid={ 'HA11' 'HA31' } ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aas={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HSD' 'HIS' 'HSE' 'GLY' 'TRP' 'CYS' ;
      'A'   'R'   'D'   'Q'   'L'   'T'   'E'   'I'   'F'   'K'   'S'   'V'   'M'   'N'   'P'   'Y'   'H'   'H'   'H'   'G'   'W'   'C' } ;
aa1= @(x) char(aas(find(ismember(aas(:), x),1)+1)); % get 1 letter code from 3 letter code
aa3= @(x) char(aas(find(ismember(aas(:), x),1)-1)); % get 3 letter code from 1 letter code

nts={ 'ADE' 'THY' 'CYT' 'GUA' 'URA' 'RDV';
      'A'   'T'   'C'   'G'   'U' 'R'} ; % add remdesivir (RDV), though nonstandard !
nt1= @(x) char(nts(find(ismember(nts(:), x),1)+1));
nt3= @(x) char(nts(find(ismember(nts(:), x),1)-1));

res=[aas nts];
res1= @(x) char(res(find(ismember(res(:), x),1)+1));
res3= @(x) char(res(find(ismember(res(:), x),1)-1));

if ~exist('qreadpdb') ; qreadpdb=1 ; end % assume that have a valid name in pdbin
%
if (qreadpdb)
 mol=pdbread(pdbin);
end
%
pdb=mol.Model.Atom;
%
resid=[pdb(:).resSeq];
segid={pdb.segID}' ;
insertion=char({pdb.iCode});
if isempty(insertion)
 insertion(:,1)=' '; % compat
end
rname={pdb(:).resName};

mstr='';

qseg=exist('msegid');

for i=1:length(mresid)
 id=mresid(i);
 ins=char(minsert(i));
 if isempty(ins)
  ins(:,1)=' '; % compat
 end
 name=char(mrname(i));
% try to convert 1-letter codes to 3-letter per PDB
 if (length(name)==1) ; name=aa3(name) ; end
 if (qseg)
  seg=char(msegid(i));
  fprintf('%s%d%s%s%s...','Looking for residue ',id,ins,'in segment ',seg)
  minds=find( resid(:)==id & insertion==ins & ismember(segid,seg)) ;
 else
  fprintf('%s%d%s...','Looking for residue ',id,ins)
  minds=find( resid(:)==id & insertion==ins ) ;
 end
 if (~isempty(minds))
  oldname=rname(minds(1));
  for ind=minds'
   pdb(ind).resName=name;
  end
  newmut=[res1(oldname),num2str(id),strtrim(ins),res1(name)];
  fprintf(['Made mutation ', newmut,'\n']);
  mstr=[mstr,'_',newmut];
 else
  fprintf('not found\n')
 end
end

if (~isempty(mstr))
 if (~exist('pdbout','var'))
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
