% process missing residue fields in the pdb file using Remark465 fields
%
r465=molecule.Remark465 ;
% find the beginning of the sequences
for i=1:size(r465,1)
 a=r465(i,:);
 ind=strfind(a,'SSSEQI');
 if ~isempty(ind)
  break;
 end
end
%
if isempty(ind)
 disp(' Could not find missing residue record ... ');
 return
end

miss=r465(i+1:end,:);
mnum=size(miss,1);
%[mrname, mchain, mresid]=strread(miss','%s %s %s'); % OK for matlab, but not octave, so reimplement below
mrname=cell(mnum,1);
mchain=cell(mnum,1);
mresid=cell(mnum,1);
for i=1:mnum ; [mrname(i), mchain(i), mresid(i)]=strread(miss(i,:),'%s %s %s'); end
%mnum=length(mrname);
% define an insertion code (if any)
mins=cell(mnum,1);
%mins=blanks(mnum)';
for i=1:mnum
 mr=char(mresid(i));
 if(isstrprop(mr(end),'alpha'))
  mresid(i)={ mr(1:end-1) };
  mins(i)={ mr(end) };
%  mins(i)=mr(end);
 else
  mins(i)={''};
 end
end
% now mresid contains only numbres so we cast it to integer type:
mresid=str2num(char(mresid));

% extract missing sequences and write dummy pdb files with missing coordinates for modeling
% define matrix with sequence limits; these are indices into 'miss'
mseqs=zeros(0,2);
nseqs=1;
i=1;
mseqs(nseqs,1)=1; % first seauence starts at one
while i<mnum
% determine if the next residue in sequence follows the current one
 qadjacent=0;
 if (strcmp(mchain(i),mchain(i+1))) % chain IDs agree, proceed
  ins1=char(mins(i));
  ins2=char(mins(i+1));
  if (isempty(ins1))
   if (isempty(ins2)) % both insertion codes blank; rank residuse by resid
    qadjacent = ( mresid(i+1)==mresid(i)+1 ) ;
   else % insertion code present on following residue ; adjacent if resids the same and the code is A
    qadjacent = ( mresid(i)==mresid(i+1) ) & strcmp(ins2,'A') ;
   end
  else % insertion code present on current residue
   if (isempty(ins2)) % inserted sequence ended, check resid sequence
    qadjacent=( mresid(i+1)==mresid(i)+1 );
   else % both insertion codes present ; adjacent if resids the same and insertion codes sequential
    qadjacent = ( mresid(i)==mresid(i+1) ) & ((ins1+1)==(ins2+0)) ;
   end
  end
 end
% fprintf('%d %d %d %d\n', i, qadjacent, ins1, ins2)
 if ~qadjacent
  mseqs(nseqs,2)=i; % terminate this sequence
  nseqs=nseqs+1;    % start new sequence
  mseqs(nseqs,1)=i+1;
 end
 i=i+1;
end
mseqs(nseqs,2)=i; % terminate final sequence
%
% now that we have the sequence limits, write dummy pdbs
%
for iseq=1:nseqs
 mpdb=struct(pdb(1));
 mx=9999;
 my=mx;
 mz=mx;
 ibeg=mseqs(iseq,1);
 iend=mseqs(iseq,2);
 mchainid=char(mchain(ibeg));
 msegid=ch2seg(mchainid);
 for i=ibeg:iend
  ind=i-ibeg+1;
  mpdb(ind).AtomSerNo=ind;
  mpdb(ind).AtomNameStruct.chemSymbol='C';
  mpdb(ind).AtomNameStruct.remoteInd='A';
  mpdb(ind).AtomNameStruct.branch='';
  mpdb(ind).AtomName='CA';
  mpdb(ind).chainID=mchainid;
  mpdb(ind).resName=char(mrname(i));
  mpdb(ind).resSeq=mresid(i);
%  cins=char(mins);
  mpdb(ind).iCode=char(mins(i));%cins(i);
  mpdb(ind).X=mx;
  mpdb(ind).Y=my;
  mpdb(ind).Z=mz;
  mpdb(ind).occupancy=0.0;
  mpdb(ind).tempFactor=0.0;
  mpdb(ind).element='';
  mpdb(ind).charge='';
  mpdb(ind).segID=msegid;
 end
% print this pdb :
 mol4.Model.Atom=fixcharmm(mpdb);
 mind=sum(char(mchain(mseqs(1:iseq,1)))==mchainid); % current number of missing loops for this chain
 mname=[strtrim(msegid),'-MISSING-',num2str(mind),'.pdb']; % missing loop file name
 pdbwrite(mname, mol4);
end

