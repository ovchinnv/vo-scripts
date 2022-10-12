% read dcd(s) average ligand position for each frame and write out new dcd
addpath('~/scripts/matlab/anal');

nrep=4 ; % number of rles replicas
name='hp35c27' ;

pdbfile = [ './',name,'_now.pdb'];

struc ;

bdcd=1; % trajectory file limits
edcd=1;

basename=[name,'_aligned_now'];
ext='dcd';
%
idcd=bdcd:edcd ;
dcdnames={};
for i = idcd
 dcdnames=[dcdnames { [basename,'.',ext] } ];
end
%
% read dcds :
traj;

% define ligand selections
ligind=[];
for irep=1:nrep
 ligind=[ligind find(ismember(segid,['REP',num2str(irep)])) ];
end

% loop over all frames and replace intantaneous ligand coordinates with averaged ones
nligatom=size(ligind,1);
for i=1:iframe
 xave=mean(reshape(xall(ligind,i),nligatom,[]),2);
 yave=mean(reshape(yall(ligind,i),nligatom,[]),2);
 zave=mean(reshape(zall(ligind,i),nligatom,[]),2);
% now put average coords into all inst. coords
 for irep=1:nrep
  xall(ligind(:,irep),i) = xave ;
  yall(ligind(:,irep),i) = yave ;
  zall(ligind(:,irep),i) = zave ;
 end
end 

% write new dcd :
writedcd([basename,'ave.dcd'], xall, yall, zall);

qdcd=0;
clear xall yall zall

