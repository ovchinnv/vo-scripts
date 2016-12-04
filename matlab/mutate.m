
pdbin='HH_2.0_47.pdb'
pdbout='HH_47'

% mutations to make (for now, a single file)
% (1) - unfolds at Cterm
mresid=[ 20 26 28 ];
mrname={ 'ASP' 'GLN' 'GLN' };
% (2) -- further atempts to stabilize the helix neat C-term
mresid=[ 20 26 28 ];
mrname={ 'ASP' 'GLN' 'ASP' }; % unfolds

mresid=[ 20 28 ];
mrname={ 'ASP' 'ASP' }; % unstable; similar to previous case

mresid=[ 28 ];
mrname={ 'GLN' }; % original, most successful mutation

%%%%%%%%%%%%%%%%%%%%%%%
aas={ 'ALA' 'ARG' 'ASP' 'GLN' 'LEU' 'THR' 'GLU' 'ILE' 'PHE' 'LYS' 'SER' 'VAL' 'MET' 'ASN' 'PRO' 'TYR' 'HIS' 'HSD' 'HSE' 'GLY' 'TRP' 'CYS' ;
      'A'   'R'   'D'   'Q'   'L'   'T'   'E'   'I'   'F'   'K'   'S'   'V'   'M'   'N'   'P'   'Y'   'H'   'H'   'H'   'G'   'W'   'C' } ;
aa1= @(x) char(aas(find(ismember(aas(:), x))+1)); % get 1 letter code from 3 letter code
aa3= @(x) char(aas(find(ismember(aas(:), x))-1)); % get 3 letter code from 1 letter code


mol=pdbread(pdbin);
pdb=mol.Model.Atom;
resid=[pdb(:).resSeq];
rname={pdb(:).resName};

mstr='';

for i=1:length(mresid)
 id=mresid(i);
 name=char(mrname(i));
 minds=find( resid==id );
 if (~isempty(minds))
  oldname=rname(minds(1));
  for ind=minds
   pdb(ind).resName=name;
  end
  mstr=[mstr,'_',aa1(oldname),num2str(id),aa1(name)];
 end
end

if (~isempty(mstr))
 if (isempty(pdbout))
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

