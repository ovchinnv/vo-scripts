%
addpath '~/scripts/matlab/matdcd-1.0' ;
%
% read pdb file
% skip if already read
if (exist('qpsf')) ; if (qpsf==1) ; return ; end ; end

if (~exist('qpdb')) ; qpdb=0 ; end

if (~qpdb) % read pdb file
 global molecule;
 disp(['==>Reading pdb structure file ',pdbfile,' ...']);
 if (~exist('qoctpdb')) ; qoctpdb=0 ; end
 if(qoctpdb)
  molecule=readpdb(pdbfile,1); % custom partial octave-compatible (faster?) ; 2nd argument for verbosity ; might need to add a path
 else
  molecule=pdbread(pdbfile); % native matlab
 end
 qpdb=1; %
end

% extract some useful information from pdb (and change the data structure)
pdb=molecule.Model.Atom ; % this is the actual pdb file
natom=length(pdb);
%
anum= [pdb.AtomSerNo]';
aname={pdb.AtomName}' ;
rname={pdb.resName}' ;
segid={pdb.segID}' ;
resid=[pdb.resSeq]' ;
insertion={pdb.iCode}';
%insertion=char({pdb.iCode});
occu =[pdb.occupancy]' ;
bet  =[pdb.tempFactor]' ;
xpdb =[pdb.X]';
ypdb =[pdb.Y]';
zpdb =[pdb.Z]';

anumber=zeros(natom,1) ; % atomic number
mass=zeros(natom,1) ; % mass

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
sulphur=find(element=='S');
sulfur=sulphur;

anumber(proton)=1;
anumber(carbon)=6;
anumber(nitrogen)=7;
anumber(oxygen)=8;
anumber(phosphorus)=15;
anumber(sulfur)=16;

mass(proton)=1.00794;
mass(carbon)=12.0107;
mass(nitrogen)=14.00674;
mass(oxygen)=15.9994;
mass(phosphorus)=30.973761;
mass(sulfur)=32.065;

% make sure everything is assigned :
unkn=find(anumber==0);
if ( unkn )
 warning('The following atoms do not have an atomic number assigned:');
 unkn
 aname(unkn)
end
%
qpsf=1;
