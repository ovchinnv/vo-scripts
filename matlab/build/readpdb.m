function mol=readpdb(fname,qloud);
% octave/matlab replacement for MATLAB's pdbread, which is very slow for large pdbs ;
% note that this isn't perfect either -- it is slower than the matlab reader for small files
% most fields are not read, the list is updated on demand
%
if (~exist('qloud')) ; qloud=0 ; end
%if (qread)
%
if (qloud) ; fprintf('Reading file "%s"\n', fname); tic ; end
fid=fopen(fname,'r');
%ncol=76 ;
%file=fscanf(fid,['%',num2str(ncol),'c']);
file=fscanf(fid,'%c'); % perhaps the fastest way to read
data=strsplit(file,{'\n'});

nlines=numel(data);
if (qloud) ; toc ; fprintf('File contains %d lines\n', nlines);end
% for each line, look at the header, then can group lines into categories for processing
% this might be optimizable without using cell functions
if (qloud) ; fprintf('Sorting records based on prefix\n');end
header=cellfun(@(l) l(1:min(6,length(l))), data, 'UniformOutput',0);
header=upper(strtrim(header)) ;

iremark=strcmp(header,'REMARK');
imodel=strcmp(header,'MODEL');
iendmdl=strcmp(header,'ENDMDL');
iatom=strcmp(header,'ATOM');
ihetatm=strcmp(header,'HETATM');
ilink=strcmp(header,'LINK');
iterminal=strcmp(header,'TER');

% initialize model
mol=struct('Model',{{}});
%=====================================================================
% process different types of record based on the index above
% add text from remark fields:
if (qloud) ; toc ; fprintf('Processing REMARK fields\n');end
dchar=char(data(find(iremark))); % grab lines that begin with REMARK in header
iremark=[];
remark=[];
if ~isempty(dchar)
 iremark=dchar(:,8:11); % remark number : VO changed from 8:10, and below from 11:80
 remark=dchar(:,12:80); % rest of remark line
end
ind=0;
for i=1:size(iremark,1);
 ir=str2num(iremark(i,:)) ;
 if(isempty(ir)) ; ir=999-ind ; ind=ind+1; % automatically generate remark number, starting from 999 (compat)
  field=['Remark',num2str(ir)];
  mol=setfield(mol(1),field,strtrim([iremark(i,:),remark(i,:)])); % restore entire line after 'REMARK'
 else
  field=['Remark',num2str(ir)];
  if (isfield(mol(1),field))
   fieldval=eval(['mol(1).',field,';']);
   fieldval = [ fieldval ; remark(i,:) ]; % update field value
  else
   fieldval = remark(i,:); % set field value
  end
  mol=setfield(mol(1),field,fieldval);
 end
end
%=====================================================================
% process LINK fields for glyco to work
if (qloud) ; toc ; fprintf('Processing LINK fields\n');end
dchar=char(data(find(ilink)));
nlink=size(dchar,1);
Link=repmat(struct(),1,nlink);
% there are usually few link fields, so we do not have to be 100% efficient
for ilink=1:nlink % NOTE: perhaps the simplest way, which might need mods in the future
 Link(ilink).remove1=dchar(ilink,13:16);
 Link(ilink).altLoc1=dchar(ilink,17);
 Link(ilink).resName1=dchar(ilink,18:20);
 Link(ilink).chainID1=dchar(ilink,22);
 Link(ilink).resSeq1=str2double(dchar(ilink,23:26));
 Link(ilink).iCode1=dchar(ilink,27);
 Link(ilink).AtomName2=dchar(ilink,43:46);
 Link(ilink).altLoc2=dchar(ilink,47);
 Link(ilink).resName2=dchar(ilink,48:50);
 Link(ilink).chainID2=dchar(ilink,52);
 Link(ilink).resSeq2=str2double(dchar(ilink,53:56));
 Link(ilink).iCode2=dchar(ilink,57);
 Link(ilink).sym1=dchar(ilink,60:65);
 Link(ilink).sym2=dchar(ilink,67:72);
end
mol(1).Link=Link;
% from matlab :
%            NumOfLink = NumOfLink+1;
%            PDB_struct.Link(NumOfLink) = ...
%                struct('remove1',{tline(13:16)},...
%                'altLoc1',{tline(17)},...
%                'resName1',{tline(18:20)},...
%                'chainID1',{tline(22)},...
%                'resSeq1',{str2double(tline(23:26))},...
%                'iCode1',{tline(27)},...
%                'AtomName2',{tline(43:46)},...
%                'altLoc2',{tline(47)},...
%                'resName2',{tline(48:50)},...
%                'chainID2',{tline(52)},...
%                'resSeq2',{str2double(tline(53:56))},...
%                'iCode2',{tline(57)},...
%                'sym1',{tline(60:65)},...
%                'sym2',{tline(67:72)});
%=====================================================================
% ATOM & HETATM :
% for now, treating as one model
%fmt='(A6,I5,1X,A4,1X,A4,1X,A5,3X,3F8.3,2F6.2,1X,3X,2X,A4)'
%  mpdb(ind).AtomSerNo=ind;
%  mpdb(ind).AtomNameStruct.chemSymbol='O';
%  mpdb(ind).AtomNameStruct.remoteInd='H';
%  mpdb(ind).AtomNameStruct.branch='2';
%  mpdb(ind).AtomName='OH2';
%  mpdb(ind).altLoc='A' ; % alternative coordinates for this atom
%  mpdb(ind).chainID=mchainid;
%  mpdb(ind).resName='TP3';
%  mpdb(ind).resSeq=ires; % assume 1 atom per residue
%  mpdb(ind).iCode=ins;
%  mpdb(ind).X=x(iatom);
%  mpdb(ind).Y=y(iatom);
%  mpdb(ind).Z=z(iatom);
%  mpdb(ind).occupancy=0.0;
%  mpdb(ind).tempFactor=0.0;
%  mpdb(ind).element='';
%  mpdb(ind).charge='';
%  mpdb(ind).altLoc='';
%  mpdb(ind).segID=msegids{ifile}(1:4);

% 5/25 : ATOM and HETATM records are identical except for the 6-char header
% so treat them together below, hopefully not too slowly
afields={'ATOM', 'HETATM'};
iatoms={ iatom, ihetatm };
strucnames={'Atom','HeterogenAtom'};

for iatype=1:2

if (qloud) ; toc ; fprintf(['Processing ',afields{iatype},' records\n']);end
if (qloud) ; toc ; fprintf('Generating character matrix\n');end
dchar=char(data(find(iatoms{iatype})));
if(numel(dchar)==0) ; continue ; end
if (qloud) ; toc ; fprintf('serial...');end
serial=dchar(:,7:11);if (~any(serial(:)=='*')) ; serial=str2num(serial) ; end
if (qloud) ; toc ; fprintf('aname...');end
aname=dchar(:,13:16); % this produces a string with spaces, which needs to be trimmed, which is a slow serial operation, so need to optimize
%anmtrim=strtrim(aname) ; [~,anamelen]=find(anmtrim==' ');j=j-1; % will not work if there are no blanks !
anmtrim=strtrim(aname) ; anamelen=sum(anmtrim~=' ',2); % convert to 0/1 sum nonzero chars (assume consecutiveness)
if (qloud) ; toc ; fprintf('altloc...');end
altloc=dchar(:,17);
if (qloud) ; toc ; fprintf('resname...');end
resname=dchar(:,18:21);
rnmtrim=strtrim(resname) ; rnamelen=sum(rnmtrim~=' ',2); % as above
if (qloud) ; toc ; fprintf('chain...');end
chain=dchar(:,22);
if (qloud) ; toc ; fprintf('resnum...');end
%resnum=str2num(dchar(:,23:26));
resnum=(dchar(:,23:26));
if (qloud) ; toc ; fprintf('ins...');end
ins=dchar(:,27);
if (qloud) ; toc ; fprintf('xx...');end
xx=str2num(dchar(:,31:38)); % what goes wrong here ?
%xx=(dchar(:,31:38));
if (qloud) ; toc ; fprintf('yy...');end
yy=str2num(dchar(:,39:46));
%yy=(dchar(:,39:46));
if (qloud) ; toc ; fprintf('zz...');end
zz=str2num(dchar(:,47:54));
%zz=(dchar(:,47:54));
if (qloud) ; toc ; fprintf('occ...');end
occ=str2num(dchar(:,55:60));
%occ=(dchar(:,55:60));
if (qloud) ; toc ; fprintf('bet...');end
bet=str2num(dchar(:,61:66));
%bet=(dchar(:,61:66));
if (qloud) ; toc ; fprintf('segid...');end
segid=dchar(:,73:76);
%sidtrim=strtrim(segid) ; sidlen=sum(sidtrim~=' ',2); % does not work is segid is blank, which is common ! So use cells
natom=size(serial,1);
% data conversions to match MATLAB pdb reader:
%if (qloud) ; toc ; fprintf('Trimming text strings\n');end
%anamec=strtrim(mat2cell(aname,ones(natom,1),4)); % conversion to cell; FASTER THAN STRTRIM, but slightly slower than explicit length
%resname=strtrim(mat2cell(resname,ones(natom,1),4));
segid=strtrim(mat2cell(segid,ones(natom,1),4));
ins=strtrim(mat2cell(ins,ones(natom,1),1));
%aname=strtrim(aname) ;
%resname=strtrim(resname) ;
%segid=strtrim(segid) ;
%chain=strtrim(chain) ; % will destroy if blank
% now, see how fast we can populate a MATLAB-style pdb structure
if (qloud) ; toc ; fprintf('Transforming data to Matlab format\n');end
% see if preallocating helps -- not really
dummyatom=struct('AtomSerNo',{},'AtomName',{},'altLoc',{},'resName',{},'chainID',{},'resSeq',{},'iCode',{},...
                'X',{},'Y',{},'Z',{},'occupancy',{},'tempFactor',{},'segID',{},'element',{},'charge',{},'AtomNameStruct',{});
%Atom=repmat(dummyatom,1,natom); % VERY SLOW
%Atom=repmat(dummyatom,natom,1); % VERY SLOW
%Atom=struct(); % generally, faster
%Atom=repmat(struct(),natom,1); % faster
Atom=repmat(struct(),1,natom); % often, fastest
AtomNameStruct=struct('chemSymbol','','remoteInd','','branch','');
for i=1:natom
% this loop is _terrible_ in octave; to the point that we cannot use it !
 if (qloud) ; if (mod(i,1000)==0) ; fprintf('%d / %d atoms processed\n',i,natom); end ; end
 Atom(i).AtomSerNo=serial(i,:);
% Atom(i).AtomName=strtrim(aname(i,:)); % I think the strtrim function is very slow, esp. when used serially
 Atom(i).AtomName=anmtrim(i,1:anamelen(i));
% Atom(i).AtomName
% Atom(i).AtomName=anamec{i};
%
% ioff=2-abs(sign(aname(i,1)-32)) ; % allow to ignore the 1st char if it is blank
% AtomNameStruct.chemSymbol=aname(i,ioff:2);
 AtomNameStruct.chemSymbol=anmtrim(i,1);
 AtomNameStruct.remoteInd=aname(i,3);
 AtomNameStruct.branch=aname(i,4);
%continue
%
 Atom(i).altLoc=altloc(i); % this ignores alternative coords
% Atom(i).resName=strtrim(resname(i,:)); % same strategy as above for aname
 Atom(i).resName=rnmtrim(i,1:rnamelen(i));
% Atom(i).resName=resname{i};
 Atom(i).chainID=chain(i);
% Atom(i).resSeq=resnum(i);
 Atom(i).resSeq=sscanf(resnum(i,:),'%d');
% Atom(i).iCode=strtrim(ins(i));
% Atom(i).iCode=strtrim(ins(i));
 Atom(i).iCode=ins{i};
% Atom(i).segID=strtrim(segid(i,:)); % slow
% Atom(i).segID=sidtrim(i,1:sidlen(i)); % does not work with blank segids
 Atom(i).segID=segid{i};
 Atom(i).element='';
 Atom(i).charge='  ';
 Atom(i).AtomNameStruct=AtomNameStruct ;
%continue
 Atom(i).X=xx(i);
 Atom(i).Y=yy(i);
 Atom(i).Z=zz(i);
 Atom(i).occupancy=occ(i);
 Atom(i).tempFactor=bet(i);
%continue
% Atom(i).X=sscanf(xx(i,:),'%f'); % this adds at least 10% extra time
% Atom(i).Y=sscanf(yy(i,:),'%f');
% Atom(i).Z=sscanf(zz(i,:),'%f');
% Atom(i).occupancy=sscanf(occ(i,:),'%f');
% Atom(i).tempFactor=sscanf(bet(i,:),'%f');
end
%
eval(['mol(1).Model.',strucnames{iatype},'=Atom;']) ; %clear Atom ;
%
end % iatype
% processs terminal records
if (qloud) ; toc ; fprintf('Processing TER records\n');end
dchar=char(data(find(iterminal)));
[numter,maxcol]=size(dchar);
if (maxcol>=11) 
 serial=dchar(:,7:11);if (~any(serial(:)=='*')) ; serial=str2num(serial) ; end
else
 serial=repmat('',numter,1);
end
if (maxcol>20) ; resname=dchar(:,18:21) ; else ; resname=repmat('',numter,1);end
if (maxcol>21) ; chain=dchar(:,22) ; else ; chain=repmat('',numter,1);end
if (maxcol>25) ; resnum=str2num(dchar(:,23:26)); else ; resnum=repmat('',numter,1);end
if (maxcol>26) ; ins=dchar(:,27); else ; ins=repmat('',numter,1);end
%
if (qloud) ; toc ; fprintf('Transforming data to match Matlab format\n');end
Terminal=struct();
for i=1:numter
 Terminal(i).SerialNo=serial(i,:);
 Terminal(i).resName=strtrim(resname(i,:));
 Terminal(i).chainID=chain(i);
 Terminal(i).resSeq=resnum(i);
 Terminal(i).iCode=strtrim(ins(i));
end
mol(1).Model.Terminal=Terminal ; %clear Atom ;

if (qloud) ; toc ; end
%qread=0;
%end % qread
