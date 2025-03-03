% writing an octave/matlab replacement for MATLAB's pdbread, which is very slow

fname='c15x3.pdb';
%fname='spy-np60_m.pdb';

qloud=1;
if (~exist('qloud')) ; qloud=0 ; end

if (~exist('qread')) ; qread=1 ; end
if (qread)

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
iterminal=strcmp(header,'TER');

% initialize model
mol=struct('Model',{{}});
% add text from remark fields:
if (qloud) ; toc ; fprintf('Processing REMARK fields\n');end
dchar=char(data(find(iremark)));
iremark=dchar(:,8:10);
remark=dchar(:,11:80);
ind=0;
for i=1:size(iremark,1);
 ir=str2num(iremark(i,:)) ;
 if(isempty(ir)) ; ir=999-ind ; ind=ind+1;
  field=['Remark',num2str(ir)];
  mol=setfield(mol(1),field,strtrim([iremark(i,:),remark(i,:)]));
 else
  field=['Remark',num2str(ir)];
  mol=setfield(mol(1),field,remark(i,:));
 end
end
% process different types of record based on the index above
% start with ATOM, which are most important
% for now, treating as one model
%fmt='(A6,I5,1X,A4,1X,A4,1X,A5,3X,3F8.3,2F6.2,1X,3X,2X,A4)'
%  mpdb(ind).AtomSerNo=ind;
%  mpdb(ind).AtomNameStruct.chemSymbol='O';
%  mpdb(ind).AtomNameStruct.remoteInd='H';
%  mpdb(ind).AtomNameStruct.branch='2';
%  mpdb(ind).AtomName='OH2';
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

if (qloud) ; toc ; fprintf('Extracting data from ATOM records\n');end
%
dchar=char(data(find(iatom)));
serial=dchar(:,7:11);if (~any(serial(:)=='*')) ; serial=str2num(serial) ; end
aname=dchar(:,13:16);
resname=dchar(:,18:21);
chain=dchar(:,22);
resnum=str2num(dchar(:,23:26));
ins=dchar(:,27);
xx=str2num(dchar(:,31:38));
yy=str2num(dchar(:,39:46));
zz=str2num(dchar(:,47:54));
occ=str2num(dchar(:,55:60));
bet=str2num(dchar(:,61:66));
segid=dchar(:,73:76);
natom=size(serial,1);
% data conversions to match MATLAB pdb reader:
%if (qloud) ; toc ; fprintf('Trimming text strings\n');end
%aname=strtrim(mat2cell(aname,ones(natom,1),4));
%resname=strtrim(mat2cell(resname,ones(natom,1),4));
%segid=strtrim(mat2cell(segid,ones(natom,1),4));
%ins=strtrim(mat2cell(ins,ones(natom,1),1));
%aname=strtrim(aname) ;
%resname=strtrim(resname) ;
%segid=strtrim(segid) ;
%chain=strtrim(chain) ; % will destroy if blank
%
% now, see how fast we can populate a MATLAB-style pdb structure
if (qloud) ; toc ; fprintf('Transforming data to match Matlab format\n');end
% see if preallocating helps -- nope
%dummyatom=struct('AtomSerNo',{},'AtomName',{},'altLoc',{},'resName',{},'chainID',{},'resSeq',{},'iCode',{},...
%                'X',{},'Y',{},'Z',{},'occupancy',{},'tempFactor',{},'segID',{},'element',{},'charge',{});
%Atom=repmat(dummyatom,natom,1);
Atom=struct();
AtomNameStruct=struct('chemSymbol','','remoteInd','','branch','');
for i=1:natom
 Atom(i).AtomSerNo=serial(i,:);
 Atom(i).AtomName=strtrim(aname(i,:));
% Atom(i).AtomName=aname{i};
%
 ioff=2-abs(sign(aname(i,1)-32)) ; % allow to ignore the 1st char if it is blank
 AtomNameStruct.chemSymbol=aname(i,ioff:2);
 AtomNameStruct.remoteInd=aname(i,3);
 AtomNameStruct.branch=aname(i,4);
%
 Atom(i).altLoc='';
 Atom(i).resName=strtrim(resname(i,:));
% Atom(i).resName=resname{i};
 Atom(i).chainID=chain(i);
 Atom(i).resSeq=resnum(i);
 Atom(i).iCode=strtrim(ins(i));
% Atom(i).iCode=ins{i};
 Atom(i).X=xx(i);
 Atom(i).Y=yy(i);
 Atom(i).Z=zz(i);
 Atom(i).occupancy=occ(i);
 Atom(i).tempFactor=bet(i);
 Atom(i).segID=strtrim(segid(i,:));
% Atom(i).segID=segid{i};
 Atom(i).element='';
 Atom(i).charge='  ';
 Atom(i).AtomNameStruct=AtomNameStruct ;
end
%
mol(1).Model.Atom=Atom ; %clear Atom ;
%
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
qread=0;
end % qread

% see if matlab will write it ...
pdbwrite(['o-',fname],mol);
