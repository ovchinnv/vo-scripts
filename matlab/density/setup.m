%
addpath '~/scripts/matlab/matdcd-1.0' ;
%
datapath='../' ;
pdbfile='dmpc.pdb'; % in lieu of a "structure" file
basename='dmpc';
flag='npt';

% npt equilibration file indices :
inds=6:12 ;

fnames={};

for i = inds
 fnames=[fnames { [datapath,basename,num2str(i),flag] } ];
end

% read pdb file

if (~exist('qpdb')) 
 qpdb=0 ; 
end

if (~qpdb) % read pdb file
 disp(['==>Reading pdb structure file ',pdbfile,'...']);
 model=pdbread([datapath,pdbfile]);
 qpdb=1; %
end

% extract some useful information from pdb (and change the data structure)
pdb=model.Model.Atom ; % this is the actual pdb file
natom=length(pdb);
%
anum=zeros(1,natom) ;
aname=cell(1,natom) ;
rname=cell(1,natom) ;
segid=cell(1,natom) ;
resid=zeros(1,natom) ;
occu=zeros(1,natom);
bet=zeros(1,natom) ;
anumber=zeros(1,natom) ; % atomic number

mass=zeros(1,natom) ; % atomic number

for i=1:natom
 atom=pdb(i);
 anum(i)=atom.AtomSerNo ;
 aname(i)={atom.AtomName} ;
 rname(i)={atom.resName} ;
 segid(i)={atom.segID} ;
 resid(i)=atom.resSeq ;
 occu(i)=atom.occupancy ;
 bet(i)=atom.tempFactor ;
end

% note : to find cell array matches, use "ismember", eg : ind = ismember(aname,'OH2') will return 1s for matches
% note : can use cell functions : e.g. to find indices with atom names beginning with O have : ind = cellfun(@(s) s(1)=='O', aname) ; % here 's' is the cell value
% in the above, the function @(s) is called with the arguments 'aname' ; it is defined inline by "s(1)==O"

% atomic element (crude) :
element=cellfun(@(s) s(1), aname);
%
% now populate atomic number array :
%
proton=find(element=='H');
carbon=find(element=='C');
nitrogen=find(element=='N');
oxygen=find(element=='O');
phosphorus=find(element=='P');

anumber(proton)=1;
anumber(carbon)=6;
anumber(nitrogen)=7;
anumber(oxygen)=8;
anumber(phosphorus)=15;

mass(proton)=1.00794;
mass(carbon)=12.0107;
mass(nitrogen)=14.00674;
mass(oxygen)=15.9994;
mass(phosphorus)=30.973761;

% make sure everything is assigned : 
find(anumber==0)
%
lipid=ismember(segid,'MEMB'); % flag all lipid atoms
dmpc=ismember(rname,'DMPC'); % dmpc
chol=ismember(rname,'CHL1'); % cholesterol

